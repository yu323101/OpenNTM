"""Filter FASTA alignments by wild-type loci and missing-data thresholds."""

from pathlib import Path
from typing import Dict, List, Set

import typer
from typer import Option


def filter_fasta(
    wild_loci_file: Path,
    fasta_file: Path,
    threshold: int,
    output_prefix: str,
    reference_name: str = "H37Rv",
) -> None:
    """
    Filter FASTA alignment based on wild-type loci and missing data.
    
    Args:
        wild_loci_file: Path to wild-type loci file (from wild-extract)
        fasta_file: Path to input FASTA file (merged sequences)
        threshold: Missing data threshold percentage (e.g., 5 for 5%)
        output_prefix: Prefix for output files
    
    Outputs:
        {output_prefix}.bak.fa: Filtered FASTA alignment
        {output_prefix}.bak.loc: Coordinates of retained positions
    """
    # Parse threshold
    thr = threshold / 100.0
    
    # Read wild-type loci
    wd = {}  # Position -> wild-type base
    loc = {}  # Position -> location info
    typer.echo(f"  Reading wild-type loci from {wild_loci_file}")
    
    with open(wild_loci_file, 'r') as f:
        for n, line in enumerate(f):
            line = line.strip()
            if line:
                cols = line.split('\t')
                if len(cols) >= 2:
                    loc[n] = cols[0]
                    wd[n] = cols[1]
    
    num_loci = len(wd)
    typer.echo(f"  Loaded {num_loci} wild-type loci")
    
    # Read FASTA sequences
    sequences = {}  # Name -> sequence
    names = []  # Order of names
    seq_length = 0
    
    typer.echo(f"  Reading FASTA sequences from {fasta_file}")
    
    with open(fasta_file, 'r') as f:
        current_name = None
        current_seq = []
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                    names.append(current_name)
                
                # Start new sequence
                current_name = line
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_name:
            sequences[current_name] = ''.join(current_seq)
            names.append(current_name)
    
    num_samples = len(sequences)
    typer.echo(f"  Loaded {num_samples} sequences")
    
    # Check sequence lengths
    lengths = set(len(seq) for seq in sequences.values())
    if len(lengths) != 1:
        typer.echo(f"❌ Error: Sequence lengths do not match: {lengths}")
        raise ValueError(f"Sequence lengths do not match: {lengths}")
    
    seq_length = lengths.pop()
    typer.echo(f"  Sequence length: {seq_length}")
    
    # Check if sequence length matches number of loci
    if seq_length != num_loci:
        typer.echo(f"⚠️  Warning: Sequence length ({seq_length}) != number of loci ({num_loci})")
    
    # Count missing/ambiguous bases at each position
    mis = {}  # Position -> count of missing/ambiguous
    typen = set()  # Positions with variation
    
    for name, seq in sequences.items():
        type_at_pos = {}  # Position -> base at this position
        
        for i in range(min(seq_length, num_loci)):
            nuc = seq[i]
            
            # Check if ambiguous (N or ?)
            if nuc in ['N', '?', 'n']:
                mis[i] = mis.get(i, 0) + 1
            else:
                # Check for variation
                if i in type_at_pos and type_at_pos[i] != nuc:
                    typen.add(i)
                else:
                    type_at_pos[i] = nuc
    
    # Determine which positions to keep
    out = {}  # Position -> 1 (keep) or 0 (filter)
    
    for i in range(min(seq_length, num_loci)):
        if i in mis:
            gap = mis[i] / num_samples
            if gap >= thr:
                out[i] = 0  # Filter out (too much missing data)
            else:
                out[i] = 1  # Keep
        else:
            out[i] = 1  # Keep (no missing data)
    
    # Count retained positions
    retained = sum(1 for v in out.values() if v == 1 and v in typen)
    typer.echo(f"  Retained {retained} positions after filtering")
    
    # Write filtered FASTA
    output_fa = Path(f"{output_prefix}.bak.fa")
    typer.echo(f"  Writing filtered FASTA to {output_fa}")
    
    with open(output_fa, 'w') as f:
        # Write each sample sequence
        for name in names:
            f.write(f"{name}\n")
            seq = sequences[name]
            
            for k in range(min(seq_length, num_loci)):
                if out.get(k, 0) == 1 and k in typen:
                    base = seq[k]
                    f.write(base)
            
            f.write('\n')
        
        # Write the reference sequence (wild-type bases), labelled with the
        # name of the reference actually used for variant calling.
        f.write(f">{reference_name}\n")
        for k in range(min(seq_length, num_loci)):
            if out.get(k, 0) == 1 and k in typen:
                if k in wd:
                    f.write(wd[k])
        f.write('\n')
    
    # Write coordinates
    output_loc = Path(f"{output_prefix}.bak.loc")
    typer.echo(f"  Writing coordinates to {output_loc}")
    
    with open(output_loc, 'w') as f:
        for k in range(min(seq_length, num_loci)):
            if out.get(k, 0) == 1 and k in typen:
                if k in loc:
                    f.write(f"{loc[k]}\n")


def filter_fasta_cmd(
    wild_loci: Path = Option(..., "--wild-loci", "-w", help="Wild-type loci file"),
    fasta: Path = Option(..., "--fasta", "-f", help="Input FASTA file"),
    threshold: int = Option(5, "--threshold", "-t", help="Missing data threshold (%)"),
    output_prefix: str = Option("filtered", "--output-prefix", "-o", help="Output file prefix"),
    reference_name: str = Option(
        "H37Rv",
        "--reference-name",
        "-r",
        help="Name for the reference sequence row; must match the active reference",
    ),
) -> None:
    """Filter FASTA alignment based on wild-type loci (replaces 2nd_loci_filt_fa_bak.pl)."""
    filter_fasta(wild_loci, fasta, threshold, output_prefix, reference_name)
    typer.echo(f"FASTA filtering completed!")
