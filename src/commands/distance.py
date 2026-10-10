"""Calculate pairwise SNP distances with :program:`snp-dists`."""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import typer
from typer import Option

from src.utils.logging_config import get_logger

logger = get_logger("commands.distance")


def build_snp_dists_command(
    alignment_file: Path,
    threads: int = 1,
    executable: str = "snp-dists",
) -> List[str]:
    """Build the external command used for compatible SNP distances.

    ``snp-dists`` defaults to counting only A/C/G/T and ignoring ambiguous
    bases and gaps.  ``-k`` is also intentional: the legacy Python code did
    not uppercase input, so lower-case masked bases remain ignored here.
    ``-L -m -t`` computes only unique pairs and gives us a parseable
    three-column stream; the converter below restores the legacy
    upper-triangle output format.
    """
    if threads < 1:
        raise ValueError("threads must be >= 1")

    return [
        executable,
        "-q",
        "-k",
        "-L",
        "-m",
        "-t",
        "-j",
        str(threads),
        str(alignment_file),
    ]


def read_alignment_names(alignment_file: Path) -> List[str]:
    """Return FASTA record names in their original order."""
    names: List[str] = []
    with open(alignment_file, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(">"):
                continue
            name = line[1:].strip()
            if not name:
                raise ValueError(f"Empty FASTA record name on line {line_number}")
            if any(character.isspace() for character in name):
                raise ValueError(
                    "FASTA record names used for distance calculation must not contain "
                    f"whitespace: {name!r}"
                )
            if name in names:
                raise ValueError(f"Duplicate FASTA record name: {name!r}")
            names.append(name)

    if not names:
        raise ValueError(f"No FASTA records found in {alignment_file}")
    return names


def _parse_molten_rows(molten_output: str) -> Dict[Tuple[str, str], int]:
    """Parse ``snp-dists -L -m -t`` output into an ordered-pair lookup."""
    distances: Dict[Tuple[str, str], int] = {}
    data_seen = False

    for line_number, raw_line in enumerate(molten_output.splitlines(), start=1):
        if not raw_line.strip():
            continue

        fields = next(csv.reader([raw_line], delimiter="\t"))
        if len(fields) != 3:
            # This fallback makes the parser useful with whitespace-formatted
            # fixtures, while retaining tabs for sample names in real output.
            fields = raw_line.split()
        if len(fields) != 3:
            raise ValueError(
                f"Invalid snp-dists molten output on line {line_number}: {raw_line!r}"
            )

        left, right, distance_text = (field.strip() for field in fields)
        try:
            distance = int(distance_text)
        except ValueError:
            if not data_seen and left.lower() in {"sequence_1", "sample_1"}:
                continue
            raise ValueError(
                f"Invalid distance on snp-dists output line {line_number}: "
                f"{distance_text!r}"
            ) from None

        if not left or not right:
            raise ValueError(f"Empty sample name on snp-dists output line {line_number}")
        if distance < 0:
            raise ValueError(f"Negative distance on snp-dists output line {line_number}")

        key = (left, right)
        previous = distances.get(key)
        if previous is not None and previous != distance:
            raise ValueError(f"Conflicting distances for ordered pair {key!r}")
        distances[key] = distance
        data_seen = True

    if not distances:
        raise ValueError("snp-dists produced no distance rows")
    return distances


def molten_to_upper_triangle(
    molten_output: str,
    output: Path,
    sample_order: Optional[Sequence[str]] = None,
) -> None:
    """Convert ``snp-dists -L -m -t`` output to the legacy three-column TSV.

    The legacy command emits each unordered pair once, in FASTA order, as
    ``sample_i<TAB>sample_j<TAB>distance`` for ``i < j``.  ``snp-dists``
    ``snp-dists`` molten output contains the unique triangle (including the
    diagonal), so we validate it and emit the same upper-triangle stream for
    downstream tools.  Full molten output is also accepted for compatibility
    with fixtures and manual callers.
    """
    distances = _parse_molten_rows(molten_output)

    if sample_order is None:
        inferred: List[str] = []
        for left, right in distances:
            if left not in inferred:
                inferred.append(left)
            if right not in inferred:
                inferred.append(right)
        sample_order = inferred
    else:
        sample_order = list(sample_order)

    if not sample_order or len(sample_order) != len(set(sample_order)):
        raise ValueError("sample_order must contain unique, non-empty sample names")

    expected = set(sample_order)
    observed = {name for pair in distances for name in pair}
    unknown = observed - expected
    missing = expected - observed
    if unknown:
        raise ValueError(f"snp-dists output contains unknown samples: {sorted(unknown)!r}")
    if missing:
        raise ValueError(f"snp-dists output is missing samples: {sorted(missing)!r}")

    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for index, left in enumerate(sample_order):
            for right in sample_order[index + 1 :]:
                forward = distances.get((left, right))
                reverse = distances.get((right, left))
                if forward is None and reverse is None:
                    raise ValueError(f"snp-dists output is missing pair {left!r}, {right!r}")
                if forward is not None and reverse is not None and forward != reverse:
                    raise ValueError(f"Asymmetric distances for pair {left!r}, {right!r}")
                writer.writerow([left, right, forward if forward is not None else reverse])


def calculate_distance(
    alignment_file: Path,
    output: Path,
    threads: int = 1,
    executable: str = "snp-dists",
) -> None:
    """Calculate pairwise SNP distances from an alignment using ``snp-dists``."""
    logger.info(
        "Calculating pairwise distance from %s to %s with %s thread(s)",
        alignment_file,
        output,
        threads,
    )
    names = read_alignment_names(alignment_file)
    command = build_snp_dists_command(alignment_file, threads, executable)

    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            "snp-dists was not found; install it with the mtb-evo conda environment"
        ) from exc
    except subprocess.CalledProcessError as exc:
        details = (exc.stderr or exc.stdout or "").strip()
        message = "snp-dists failed"
        if details:
            message = f"{message}: {details}"
        raise RuntimeError(message) from exc

    molten_to_upper_triangle(result.stdout, output, names)


def distance_cmd(
    alignment: Path = Option(..., "--alignment", "-a", help="Alignment FASTA file"),
    output: Path = Option(..., "--output", "-o", help="Output distance file"),
    threads: int = Option(
        1,
        "--threads",
        "-t",
        min=1,
        help="Threads passed to snp-dists (default: 1)",
    ),
) -> None:
    """Calculate pairwise SNP distances with the configured thread count."""
    calculate_distance(alignment, output, threads)
    logger.info("Distance calculation completed: %s", output)
    typer.echo(f"Distance calculation completed: {output}")
