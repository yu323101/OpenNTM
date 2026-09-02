"""Screen annotated variants for known TB drug-resistance markers.

This module reports a genotype-based resistance risk screen.  Absence of a
marker is deliberately reported as ``NO_KNOWN_MARKER`` rather than susceptible:
the supplied database is not a complete substitute for phenotypic DST or a
validated targeted NGS assay.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import typer
from typer import Option

from src.core.report_inputs import ANNOTATED_VARIANTS_DIRNAME, LEGACY_ANNOTATED_VARIANTS_DIRNAME
from src.utils.logging_config import get_logger

logger = get_logger("commands.drug_resistance")
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DRUG_DB = PROJECT_ROOT / "data" / "drug_db" / "Total_drug_resistance_mutation.txt"

DRUG_ORDER = [
    "RIFAMPICIN",
    "ISONIAZID",
    "FLUOROQUINOLONES",
    "BEDAQUILINE",
    "LINEZOLID",
    "PYRAZINAMIDE",
    "ETHAMBUTOL",
    "STREPTOMYCIN",
    "ETHIONAMIDE",
    "AMIKACIN",
    "CAPREOMYCIN",
    "KANAMYCIN",
    "PARA-AMINOSALISYLIC_ACID",
    "CLOFAZIMINE",
]
DRUG_ALIASES = {
    "FLUOROQUINOLONE": "FLUOROQUINOLONES",
    "FLUOROQUINOLONES": "FLUOROQUINOLONES",
    "FQ": "FLUOROQUINOLONES",
    "RIF": "RIFAMPICIN",
    "INH": "ISONIAZID",
    "PZA": "PYRAZINAMIDE",
    "EMB": "ETHAMBUTOL",
    "AMK": "AMIKACIN",
    "CAP": "CAPREOMYCIN",
    "KM": "KANAMYCIN",
    "PAS": "PARA-AMINOSALISYLIC_ACID",
}
FQ_DRUGS = {"FLUOROQUINOLONES"}
GROUP_A_DRUGS = {"BEDAQUILINE", "LINEZOLID"}
RESULT_SCOPE = "MTB_TB_ONLY"


@dataclass(frozen=True)
class ResistanceMarker:
    drug: str
    locus: str
    position: str
    nucleotide_change: str
    source: str

    @property
    def marker_id(self) -> str:
        return f"{self.drug}:{self.position}:{self.nucleotide_change}"


def normalize_drug(value: str) -> str:
    cleaned = value.strip().upper().replace(" ", "_")
    return DRUG_ALIASES.get(cleaned, cleaned)


def _position_tokens(value: str) -> List[str]:
    return re.findall(r"\d+", value or "")


def _change_parts(value: str) -> Optional[Tuple[str, str]]:
    if "/" not in value:
        return None
    ref, alt = value.strip().upper().split("/", 1)
    if not ref or not alt or "-" in {ref, alt}:
        return None
    return ref, alt


def load_resistance_markers(drug_db: Path) -> Dict[str, List[ResistanceMarker]]:
    """Index database markers by genomic coordinate."""
    by_position: Dict[str, List[ResistanceMarker]] = {}
    with open(drug_db, "r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if line_number == 1 or len(row) < 5:
                continue
            drug = normalize_drug(row[0])
            positions = _position_tokens(row[2])
            if not positions or not _change_parts(row[4]):
                continue
            marker = ResistanceMarker(
                drug=drug,
                locus=row[1].strip(),
                position=row[2].strip(),
                nucleotide_change=row[4].strip().upper(),
                source=row[9].strip() if len(row) > 9 else "",
            )
            for position in positions:
                by_position.setdefault(position, []).append(marker)
    return by_position


def _variant_alleles(ref: str, alt: str) -> Set[Tuple[str, str]]:
    """Return normalized allele forms for SNP/indel rows from the annotation step."""
    ref = ref.strip().upper()
    alt = alt.strip().upper()
    if not ref or not alt:
        return set()
    if len(alt) == 1 and alt in "ACGT":
        return {(ref, alt)}
    if alt.startswith("+") and len(alt) > 1:
        inserted = alt[1:]
        return {(ref, ref + inserted), ("", inserted)}
    if alt.startswith("-") and len(alt) > 1:
        deleted = alt[1:]
        return {(ref, ref.replace(deleted, "", 1)), (deleted, "")}
    return {(ref, alt)}


def _matched_markers(annotated_file: Path, markers_by_position: Dict[str, List[ResistanceMarker]]) -> List[ResistanceMarker]:
    matched: Dict[str, ResistanceMarker] = {}
    with open(annotated_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) < 3 or not _position_tokens(row[0]):
                continue
            alleles = _variant_alleles(row[1], row[2])
            if not alleles:
                continue
            for position in _position_tokens(row[0]):
                for marker in markers_by_position.get(position, []):
                    marker_alleles = _change_parts(marker.nucleotide_change)
                    if marker_alleles and marker_alleles in alleles:
                        matched[marker.marker_id] = marker
    return sorted(matched.values(), key=lambda item: (DRUG_ORDER.index(item.drug) if item.drug in DRUG_ORDER else 999, item.marker_id))


def classify_resistance(drugs: Iterable[str]) -> str:
    resistant = {normalize_drug(drug) for drug in drugs}
    if {"ISONIAZID", "RIFAMPICIN"}.issubset(resistant):
        if FQ_DRUGS.intersection(resistant) and GROUP_A_DRUGS.intersection(resistant):
            return "XDR"
        if FQ_DRUGS.intersection(resistant):
            return "pre-XDR"
        return "MDR"
    if "RIFAMPICIN" in resistant:
        return "RR"
    if "ISONIAZID" in resistant:
        return "ISONIAZID_RESISTANT"
    if not resistant:
        return "NO_KNOWN_MARKER"
    if len(resistant) == 1:
        return "SINGLE_DRUG_RESISTANCE"
    return "POLYRESISTANCE"


def _interpretation(classification: str) -> str:
    return {
        "NO_KNOWN_MARKER": "未检出数据库收录的耐药突变；不能据此判定敏感",
        "SINGLE_DRUG_RESISTANCE": "检出单一药物相关耐药标记，需结合药敏确认",
        "ISONIAZID_RESISTANT": "异烟肼耐药风险标记，需结合药敏确认",
        "RR": "利福平耐药（RR-TB）风险标记，需结合药敏确认",
        "MDR": "MDR-TB 风险标记（异烟肼+利福平），需结合药敏确认",
        "pre-XDR": "pre-XDR-TB 风险标记，需结合药敏确认",
        "XDR": "XDR-TB 风险标记，需结合药敏确认",
        "POLYRESISTANCE": "检出多种药物相关耐药标记，需结合药敏确认",
    }.get(classification, "需结合药敏确认")


def screen_resistance(annotated_file: Path, drug_db: Path, output: Path, sample_id: str = "") -> Path:
    markers_by_position = load_resistance_markers(drug_db)
    markers = _matched_markers(annotated_file, markers_by_position)
    drugs = sorted({marker.drug for marker in markers}, key=lambda drug: DRUG_ORDER.index(drug) if drug in DRUG_ORDER else 999)
    classification = classify_resistance(drugs)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "sample_id", "classification", "resistant_drugs", "resistance_mutations",
            "marker_count", "evidence", "interpretation", "organism_scope",
        ])
        writer.writerow([
            sample_id or annotated_file.name.removesuffix("_annotated.txt"),
            classification,
            ";".join(drugs) or "-",
            ";".join(marker.marker_id for marker in markers) or "-",
            len(markers),
            "known_genotypic_marker" if markers else "no_known_marker_detected",
            _interpretation(classification),
            RESULT_SCOPE,
        ])
    logger.info("Resistance screen %s: %s (%d markers)", annotated_file, classification, len(markers))
    return output


def _annotated_variant_inputs(samples_dir: Path) -> List[Path]:
    """Find canonical and legacy annotated-variant locations once each."""
    inputs = set()
    for sample_dir in samples_dir.glob("*"):
        if not sample_dir.is_dir():
            continue
        for dirname in (ANNOTATED_VARIANTS_DIRNAME, LEGACY_ANNOTATED_VARIANTS_DIRNAME):
            inputs.update((sample_dir / "report_inputs" / dirname).glob("*_annotated.txt"))
    return sorted(inputs)


def screen_samples(samples_dir: Path, drug_db: Path, output_dir: Path, manifest_output: Path) -> Path:
    """Run the MTB/TB-only screen for standardized and legacy inputs."""
    inputs = _annotated_variant_inputs(samples_dir)
    sample_output_root = output_dir if output_dir.name == "samples" else output_dir / "samples"
    rows: List[Dict[str, str]] = []
    for annotated in inputs:
        sample_id = annotated.name.removesuffix("_annotated.txt")
        sample_dir = annotated.parent.parent.parent
        output = sample_output_root / sample_dir.name / "table" / f"{sample_id}_drug_resistance.tsv"
        screen_resistance(annotated, drug_db, output, sample_id=sample_id)
        with open(output, "r", encoding="utf-8", newline="") as handle:
            rows.append(next(csv.DictReader(handle, delimiter="\t")))

    manifest_output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id", "classification", "resistant_drugs", "resistance_mutations",
        "marker_count", "evidence", "interpretation", "organism_scope",
    ]
    with open(manifest_output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if not inputs:
        logger.warning("No annotated variant files found under %s; wrote an empty manifest", samples_dir)
    return manifest_output


def drug_resistance_cmd(
    annotated: Optional[Path] = Option(None, "--annotated", help="One annotated variant TSV"),
    drug_db: Path = Option(DEFAULT_DRUG_DB, "--drug-db", help="Resistance marker database TSV"),
    output: Optional[Path] = Option(None, "--output", "-o", help="One-sample output TSV"),
    samples_dir: Optional[Path] = Option(None, "--samples-dir", help="Results/samples directory for batch mode"),
    output_dir: Optional[Path] = Option(None, "--output-dir", help="Results root for batch mode"),
    manifest_output: Optional[Path] = Option(None, "--manifest-output", help="Batch manifest TSV"),
) -> None:
    """Screen one sample or all standardized samples for known resistance markers."""
    if annotated:
        if not output:
            raise typer.BadParameter("--output is required with --annotated")
        result = screen_resistance(annotated, drug_db, output)
    elif samples_dir:
        root = output_dir or samples_dir.parent
        manifest = manifest_output or root / "logs" / "drug_resistance_manifest.tsv"
        result = screen_samples(samples_dir, drug_db, root, manifest)
    else:
        raise typer.BadParameter("provide --annotated or --samples-dir")
    typer.echo(f"Drug resistance results written to {result}")
