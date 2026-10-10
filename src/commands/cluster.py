"""Identify connected genomic proximity clusters from pairwise distances."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import typer
from typer import Option

from src.utils.logging_config import get_logger

logger = get_logger("commands.cluster")


class _UnionFind:
    def __init__(self, items: Iterable[str]) -> None:
        self.parent = {item: item for item in items}

    def find(self, item: str) -> str:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            # Stable roots make results reproducible across Python versions.
            if left_root < right_root:
                self.parent[right_root] = left_root
            else:
                self.parent[left_root] = right_root


def read_distance_pairs(distance_matrix: Path) -> Tuple[List[str], List[Tuple[str, str, int]]]:
    """Read the three-column upper-triangle distance format emitted by distance.py."""
    samples = set()
    pairs: List[Tuple[str, str, int]] = []
    with open(distance_matrix, "r", encoding="utf-8", errors="ignore") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or not fields[0].strip() or not fields[1].strip():
                continue
            left, right = fields[0].strip(), fields[1].strip()
            try:
                distance = int(fields[2])
            except ValueError:
                if line_number == 1 and fields[0].lower() in {"sample", "sample_id", "name"}:
                    continue
                raise ValueError(f"Invalid distance on line {line_number}: {fields[2]!r}")
            if distance < 0:
                raise ValueError(f"Negative distance on line {line_number}: {distance}")
            samples.update((left, right))
            pairs.append((left, right, distance))

    return sorted(samples), pairs


def assign_clusters(
    distance_matrix: Path,
    output: Path,
    threshold: int = 12,
) -> Path:
    """Write connected components of the distance graph as cluster assignments.

    A cluster is a connected component: each edge has distance <= threshold.
    This is intentionally not interpreted as proof of transmission; samples
    should be checked against collection dates, locations and epidemiological
    records before an outbreak response is initiated.
    """
    if threshold < 0:
        raise ValueError("threshold must be >= 0")

    samples, pairs = read_distance_pairs(distance_matrix)
    uf = _UnionFind(samples)
    for left, right, distance in pairs:
        if distance <= threshold:
            uf.union(left, right)

    groups: Dict[str, List[str]] = {}
    for sample in samples:
        groups.setdefault(uf.find(sample), []).append(sample)

    # Rename roots by the lexicographically smallest member so cluster IDs are
    # stable when the input row order changes.
    ordered_groups = sorted((sorted(members) for members in groups.values()), key=lambda x: x[0])
    cluster_by_sample: Dict[str, Tuple[str, int]] = {}
    for index, members in enumerate(ordered_groups, start=1):
        cluster_id = f"C{index:04d}"
        for sample in members:
            cluster_by_sample[sample] = (cluster_id, len(members))

    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "cluster_id", "cluster_size", "snp_threshold"])
        for sample in sorted(samples):
            cluster_id, size = cluster_by_sample[sample]
            writer.writerow([sample, cluster_id, size, threshold])

    logger.info("Wrote %d samples in %d clusters to %s", len(samples), len(ordered_groups), output)
    return output


def cluster_cmd(
    distance_matrix: Path = Option(..., "--distance-matrix", help="Three-column pairwise distance TSV"),
    output: Path = Option(..., "--output", "-o", help="Cluster assignment TSV"),
    threshold: int = Option(12, "--threshold", help="Maximum SNP distance for a graph edge"),
) -> None:
    """Group samples into genomic proximity clusters."""
    typer.echo(f"Cluster assignments written to {assign_clusters(distance_matrix, output, threshold)}")
