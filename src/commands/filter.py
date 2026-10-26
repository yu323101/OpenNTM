"""Filter core SNP alignments."""

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set

import typer
from typer import Option

from src.utils.logging_config import get_logger


logger = get_logger("commands.filter")


def filter_core_snp(
    wild_loci_file: Path,
    alignment_file: Path,
    threshold: int,
    output_prefix: str,
    reference_name: str = "H37Rv",
) -> None:
    """
    Filter core SNP alignment based on missing/mixed ratio and polymorphism.

    This function implements the logic from 2nd_loci_filt_fa_bak.pl:
    1. Read wild-type bases and coordinates
    2. Read all sample sequences
    3. Calculate missing/mixed ratio for each position
    4. Identify polymorphic positions
    5. Filter positions: missing/mixed < threshold AND polymorphic
    6. Output filtered alignment with the reference appended
    """
    # Step 1: Read wild-type bases
    logger.info(
        "Filtering core SNP alignment wild_loci=%s alignment=%s threshold=%s output_prefix=%s",
        wild_loci_file,
        alignment_file,
        threshold,
        output_prefix,
    )
    wild_bases: Dict[int, str] = {}  # index -> base
    loci_coords: Dict[int, int] = {}  # index -> coordinate

    with open(wild_loci_file) as f:
        for idx, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) >= 2:
                coord = int(cols[0])
                base = cols[1]
                wild_bases[idx] = base
                loci_coords[idx] = coord

    num_loci = len(wild_bases)
    logger.info("Read %d wild-type loci", num_loci)
    typer.echo(f"Read {num_loci} wild-type loci")

    # Step 2: Read sample sequences
    sequences: Dict[str, str] = {}
    sample_names: List[str] = []

    with open(alignment_file) as f:
        current_name = None
        current_seq = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_seq)
                current_name = line
                sample_names.append(current_name)
                current_seq = []
            elif line:
                current_seq.append(line)

        if current_name is not None:
            sequences[current_name] = "".join(current_seq)

    num_samples = len(sample_names)
    logger.info("Read %d sample sequences", num_samples)
    typer.echo(f"Read {num_samples} sample sequences")

    # Verify sequence lengths match
    if not sample_names:
        raise ValueError(f"No FASTA headers found in {alignment_file}")

    seq_lengths = {len(sequences.get(name, "")) for name in sample_names}
    if len(seq_lengths) > 1:
        raise ValueError(f"Sequence lengths do not match: {seq_lengths}")

    actual_length = seq_lengths.pop() if seq_lengths else 0
    if actual_length != num_loci:
        raise ValueError(
            f"Alignment length ({actual_length}) does not match "
            f"number of loci ({num_loci})"
        )

    # Step 3 & 4: Calculate missing/mixed ratio and identify polymorphic positions
    thr = threshold / 100.0
    mis_count: Dict[int, int] = defaultdict(int)
    bases_at_pos: Dict[int, Set[str]] = defaultdict(set)

    for idx in range(num_loci):
        for name in sample_names:
            base = sequences[name][idx]
            if base in "N?":
                mis_count[idx] += 1
            else:
                bases_at_pos[idx].add(base)

    # Step 5: Determine which positions to keep
    is_polymorphic: Dict[int, bool] = {
        idx: len(bases) > 1 for idx, bases in bases_at_pos.items()
    }

    keep: Dict[int, bool] = {}
    kept_positions = 0
    for idx in range(num_loci):
        gap_ratio = mis_count[idx] / num_samples
        # Keep if: gap_ratio < threshold AND is polymorphic
        if gap_ratio < thr and is_polymorphic.get(idx, False):
            keep[idx] = True
            kept_positions += 1
        else:
            keep[idx] = False

    logger.info(
        "Kept %d out of %d positions with threshold=%d%%",
        kept_positions,
        num_loci,
        threshold,
    )
    typer.echo(f"Kept {kept_positions} out of {num_loci} positions "
               f"(threshold: {threshold}%)")

    # Step 6: Output filtered alignment
    output_fa = f"{output_prefix}.fadel-InvMisF{threshold}.bak.fa"
    output_loc = f"{output_prefix}.fadel-InvMisF{threshold}.bak.loc"

    with open(output_fa, "w") as f_fa, open(output_loc, "w") as f_loc:
        # Output sample sequences
        for name in sample_names:
            f_fa.write(f"{name}\n")
            filtered_seq = "".join(
                sequences[name][idx]
                for idx in range(num_loci)
                if keep.get(idx, False)
            )
            f_fa.write(f"{filtered_seq}\n")

        # Output reference row (wild-type bases). The name must match the
        # reference the calls were made against, not a fixed H37Rv label.
        f_fa.write(f">{reference_name}\n")
        reference_seq = "".join(
            wild_bases[idx] for idx in range(num_loci) if keep.get(idx, False)
        )
        f_fa.write(f"{reference_seq}\n")

        # Output coordinates
        for idx in range(num_loci):
            if keep.get(idx, False):
                f_loc.write(f"{loci_coords[idx]}\n")

    typer.echo(f"Output files:")
    typer.echo(f"  - Alignment: {output_fa}")
    typer.echo(f"  - Coordinates: {output_loc}")


def filter_cmd(
    wild_loci: Path = Option(..., "--wild-loci", "-w", help="Wild-type loci file"),
    alignment: Path = Option(..., "--alignment", "-a", help="Alignment FASTA file"),
    threshold: int = Option(5, "--threshold", "-t", help="Minimum valid samples threshold (%)"),
    output_prefix: str = Option(
        "all_strains", "--output-prefix", "-o", help="Output file prefix"
    ),
    reference_name: str = Option(
        "H37Rv",
        "--reference-name",
        "-r",
        help="Name for the reference sequence row; must match the active reference",
    ),
) -> None:
    """Filter core SNP alignment."""
    filter_core_snp(wild_loci, alignment, threshold, output_prefix, reference_name)
    logger.info("Core SNP filtering completed: %s.fadel-InvMisF%s.bak.fa", output_prefix, threshold)
    typer.echo("Core SNP filtering completed!")
