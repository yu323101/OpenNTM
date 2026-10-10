"""Calculate pairwise SNP distances with the compiled ``snp-dists`` tool."""

from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List

import typer
from typer import Option

from src.utils.logging_config import get_logger


logger = get_logger("commands.distance")


def _read_alignment_names(alignment_file: Path) -> List[str]:
    """Read and validate FASTA names without loading the alignment sequences."""
    names: List[str] = []
    with alignment_file.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            name = line[1:].strip()
            if not name:
                raise ValueError(f"Empty FASTA header in {alignment_file}")
            names.append(name)

    if not names:
        raise ValueError(f"No FASTA records found in {alignment_file}")
    if len(names) != len(set(names)):
        raise ValueError(f"Duplicate FASTA headers found in {alignment_file}")
    return names


def _write_upper_triangle(
    lower_triangle_file: Path,
    output: Path,
    alignment_names: List[str],
) -> int:
    """Convert ``snp-dists -L -b`` output to the legacy three-column format."""
    with lower_triangle_file.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            raise ValueError("snp-dists produced an empty distance matrix")

        matrix_names = [field.strip() for field in header[1:]]
        if matrix_names != alignment_names:
            raise ValueError(
                "snp-dists sample order does not match the FASTA alignment: "
                f"expected {alignment_names[:3]}..., got {matrix_names[:3]}..."
            )

        output.parent.mkdir(parents=True, exist_ok=True)
        pair_count = 0
        row_index = 0
        with output.open("w", encoding="utf-8") as out_handle:
            for row in reader:
                if not row or not any(field.strip() for field in row):
                    continue
                if row_index >= len(matrix_names):
                    raise ValueError("snp-dists produced more rows than sequence names")

                row_name = row[0].strip()
                expected_name = matrix_names[row_index]
                if row_name != expected_name:
                    raise ValueError(
                        f"Unexpected snp-dists row {row_index + 1}: "
                        f"expected {expected_name!r}, got {row_name!r}"
                    )

                values = row[1:]
                if len(values) < row_index + 1:
                    raise ValueError(
                        f"Incomplete lower triangle for {row_name!r}: "
                        f"expected {row_index + 1} values, got {len(values)}"
                    )

                for column_index in range(row_index):
                    try:
                        distance = int(values[column_index].strip())
                    except ValueError as exc:
                        raise ValueError(
                            f"Invalid distance for {matrix_names[column_index]!r}, "
                            f"{row_name!r}: {values[column_index]!r}"
                        ) from exc
                    if distance < 0:
                        raise ValueError(f"Negative distance for {row_name!r}")
                    out_handle.write(
                        f"{matrix_names[column_index]}\t{row_name}\t{distance}\n"
                    )
                    pair_count += 1

                row_index += 1

    if row_index != len(matrix_names):
        raise ValueError(
            f"snp-dists produced {row_index} rows for {len(matrix_names)} sequences"
        )

    expected_pairs = len(alignment_names) * (len(alignment_names) - 1) // 2
    if pair_count != expected_pairs:
        raise ValueError(
            f"Expected {expected_pairs} unique pairs, converted {pair_count}"
        )
    return pair_count


def calculate_distance(
    alignment_file: Path,
    output: Path,
    threads: int = 1,
    executable: str = "snp-dists",
) -> None:
    """Run ``snp-dists`` and preserve MTB-evo's three-column output contract.

    The default ``snp-dists`` rules count only A/C/G/T mismatches and ignore
    ambiguous bases and gaps. This matches the previous Python implementation.
    The compiled tool writes a lower-triangle matrix internally; this wrapper
    converts it to the historical upper-triangle ``sample1\tsample2\tdistance``
    format consumed by downstream commands.
    """
    if threads < 1:
        raise ValueError("threads must be >= 1")
    if not alignment_file.is_file():
        raise FileNotFoundError(f"Alignment file does not exist: {alignment_file}")

    alignment_names = _read_alignment_names(alignment_file)
    executable_path = shutil.which(executable)
    if executable_path is None:
        candidate = Path(executable)
        if candidate.is_file() and candidate.stat().st_mode & 0o111:
            executable_path = str(candidate)
    if executable_path is None:
        raise RuntimeError(
            "snp-dists executable was not found. Install the workflow environment "
            "or pass --snp-dists /path/to/snp-dists."
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        executable_path,
        "-q",
        "-b",
        "-L",
        "-j",
        str(threads),
        str(alignment_file),
    ]
    logger.info(
        "Calculating pairwise SNP distance with snp-dists threads=%d: %s -> %s",
        threads,
        alignment_file,
        output,
    )

    with tempfile.TemporaryDirectory(
        prefix=".snp-dists-",
        dir=str(output.parent),
    ) as temporary_dir:
        lower_triangle_file = Path(temporary_dir) / "lower_triangle.tsv"
        with lower_triangle_file.open("w", encoding="utf-8") as raw_output:
            completed = subprocess.run(
                command,
                stdout=raw_output,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )

        if completed.returncode != 0:
            details = (completed.stderr or "").strip()
            suffix = f": {details}" if details else ""
            raise RuntimeError(
                f"snp-dists failed with exit code {completed.returncode}{suffix}"
            )

        pair_count = _write_upper_triangle(
            lower_triangle_file,
            output,
            alignment_names,
        )

    logger.info("Wrote %d unique pairwise distances to %s", pair_count, output)


def distance_cmd(
    alignment: Path = Option(..., "--alignment", "-a", help="Alignment FASTA file"),
    output: Path = Option(..., "--output", "-o", help="Output distance file"),
    threads: int = Option(
        1,
        "--threads",
        "-t",
        help="Number of snp-dists threads (default: 1; controlled by YAML in the workflow)",
    ),
    snp_dists: str = Option(
        "snp-dists",
        "--snp-dists",
        help="snp-dists executable or absolute path",
    ),
) -> None:
    """Calculate pairwise SNP distances using the compiled snp-dists tool."""
    calculate_distance(alignment, output, threads=threads, executable=snp_dists)
    logger.info("Distance calculation completed: %s", output)
    typer.echo(f"Distance calculation completed: {output}")
