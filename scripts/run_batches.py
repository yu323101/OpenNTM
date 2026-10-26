#!/usr/bin/env python3
"""Run a Snakemake workflow in bounded sample batches.

Large cohorts can exhaust disk even when every rule is written correctly: a
``temp()`` output is only deleted once the job that consumes it has finished, so
a fast producer (trimming, single-threaded) schedules far ahead of a slow
consumer (alignment, many threads) and thousands of intermediates pile up.
Submitting a bounded set of sample targets per Snakemake call caps how many of
those intermediates can exist at once, without changing any rule.

Run it from the workflow directory so that relative paths in the configuration
keep working.
"""

import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Submit per-sample Snakemake targets in bounded batches.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--samplesheet", type=Path, required=True, help="CSV with one row per sample")
    parser.add_argument("--results-dir", type=Path, required=True, help="Workflow results directory")
    parser.add_argument(
        "--target",
        default="logs/rules/core_call_variants/{sample}.done",
        help="Per-sample target, relative to --results-dir; {sample} is substituted",
    )
    parser.add_argument(
        "--prefer-input",
        default=None,
        help=(
            "Optional per-sample input, relative to --results-dir. Samples where it already "
            "exists are submitted first, so their intermediates are consumed instead of new ones."
        ),
    )
    parser.add_argument(
        "--cleanup",
        action="append",
        default=[],
        metavar="GLOB_TEMPLATE",
        help=(
            "Optional per-sample file to delete after a successful batch, relative to "
            "--results-dir; repeatable. Only files of samples whose target now exists are removed."
        ),
    )
    parser.add_argument("--sample-column", default="sample_id", help="Samplesheet column holding sample IDs")
    parser.add_argument(
        "--configfile",
        action="append",
        default=[],
        help="Passed through as --configfile; repeatable, later files win",
    )
    parser.add_argument("--batch-size", type=int, default=9, help="Samples submitted per Snakemake call")
    parser.add_argument("--cores", default="24", help="Passed through as --cores")
    parser.add_argument(
        "--max-batches", type=int, default=0, help="Stop after this many batches (0 = all)"
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help=(
            "Keep going after a failed batch: failed samples are appended to --skip-file "
            "and skipped on later runs. Stops if every sample in a batch fails."
        ),
    )
    parser.add_argument(
        "--skip-file",
        type=Path,
        default=None,
        help="Samples to skip, one per line (default: <results-dir>/failed_samples.txt)",
    )
    parser.add_argument("--snakemake", default="snakemake", help="Snakemake executable")
    parser.add_argument("--log", type=Path, default=None, help="Log file (default: <results-dir>/logs/run_batches.log)")
    parser.add_argument("--dry-run", action="store_true", help="Dry-run the first batch and exit")
    return parser.parse_args()


def read_samples(sheet: Path, column: str) -> list:
    with sheet.open(newline="") as handle:
        return [row[column].strip() for row in csv.DictReader(handle) if row.get(column)]


def resolve(results_dir: Path, template: str, sample: str) -> Path:
    return results_dir / template.format(sample=sample)


def load_skip(path: Path) -> set:
    """One sample ID per line; missing file means nothing to skip."""
    if not path.exists():
        return set()
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def record_failures(path: Path, samples) -> None:
    with path.open("a") as handle:
        for sample in samples:
            handle.write(sample + "\n")


def order_samples(samples, results_dir: Path, target_template: str, prefer_template=None) -> list:
    """Pending samples, with those whose preferred input already exists placed first."""
    pending = [s for s in samples if not resolve(results_dir, target_template, s).exists()]
    if not prefer_template:
        return pending
    ready = [s for s in pending if resolve(results_dir, prefer_template, s).exists()]
    ready_set = set(ready)
    return ready + [s for s in pending if s not in ready_set]


def main() -> int:
    args = parse_args()
    if args.batch_size < 1:
        print("--batch-size must be at least 1", file=sys.stderr)
        return 2

    args.skip_file = args.skip_file or args.results_dir / "failed_samples.txt"
    log_path = args.log or args.results_dir / "logs" / "run_batches.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    def log(message: str) -> None:
        line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}"
        print(line, flush=True)
        with log_path.open("a") as handle:
            handle.write(line + "\n")

    samples = read_samples(args.samplesheet, args.sample_column)
    skip = load_skip(args.skip_file)
    order = order_samples(
        [s for s in samples if s not in skip], args.results_dir, args.target, args.prefer_input
    )
    total_batches = (len(order) + args.batch_size - 1) // args.batch_size
    log(
        f"pending {len(order)}/{len(samples)} samples, batch size {args.batch_size}"
        f"{', dry-run' if args.dry_run else ''}"
    )
    if not order:
        return 0

    for index, start in enumerate(range(0, len(order), args.batch_size), start=1):
        if args.max_batches and index > args.max_batches:
            log(f"stopped after {args.max_batches} batch(es)")
            return 0
        batch = order[start : start + args.batch_size]
        command = [args.snakemake]
        for config in args.configfile:
            command += ["--configfile", config]
        command += ["--cores", str(args.cores), "--rerun-incomplete", "--keep-going"]
        if args.dry_run:
            command.append("-n")
        command += [str(resolve(args.results_dir, args.target, sample)) for sample in batch]

        log(f"batch {index}/{total_batches}: {len(batch)} samples")
        code = subprocess.run(command).returncode
        if code != 0:
            if not args.continue_on_error:
                log(f"batch failed with exit code {code}; no further batches submitted")
                return code
            missing = [
                s for s in batch if not resolve(args.results_dir, args.target, s).exists()
            ]
            record_failures(args.skip_file, missing)
            log(
                f"batch {index} failed with exit code {code}; recorded {len(missing)} "
                f"failed sample(s) in {args.skip_file.name} and continuing"
            )
            if len(missing) == len(batch):
                log("every sample in the batch failed; stopping in case it is systematic")
                return code
            continue
        if args.dry_run:
            log("dry-run finished")
            return 0

        freed = 0
        for sample in batch:
            if not resolve(args.results_dir, args.target, sample).exists():
                continue
            for template in args.cleanup:
                path = resolve(args.results_dir, template, sample)
                if path.exists():
                    freed += path.stat().st_size
                    path.unlink()
        if args.cleanup:
            log(f"batch done, freed {freed / 1e9:.1f} GB of intermediates")
        else:
            log("batch done")

    log(f"all {total_batches} batches finished")
    return 0


if __name__ == "__main__":
    sys.exit(main())
