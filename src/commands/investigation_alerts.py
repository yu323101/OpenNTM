"""Generate explainable public-health follow-up prompts from pipeline outputs."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List

import typer
from typer import Option

from src.commands.drug_resistance import RESULT_SCOPE
from src.utils.logging_config import get_logger

logger = get_logger("commands.investigation_alerts")


def _read_tsv(path: Path) -> List[Dict[str, str]]:
    with open(path, "r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        return [{key: str(value or "") for key, value in row.items()} for row in csv.DictReader(handle, delimiter="\t")]


def generate_alerts(cluster_assignments: Path, drug_resistance_dir: Path, output: Path) -> Path:
    """Create action prompts; these are not transmission or clinical diagnoses."""
    cluster_rows = _read_tsv(cluster_assignments)
    resistance_files = sorted(drug_resistance_dir.rglob("*_drug_resistance.tsv"))
    resistance_rows: Dict[str, Dict[str, str]] = {}
    for path in resistance_files:
        for row in _read_tsv(path):
            if row.get("sample_id"):
                resistance_rows[row["sample_id"]] = row

    alerts: List[Dict[str, str]] = []
    for row in cluster_rows:
        sample_id = row.get("sample_id", "")
        cluster_size = int(row.get("cluster_size", "0") or 0)
        cluster_id = row.get("cluster_id", "")
        if cluster_size >= 2:
            alerts.append({
                "alert_type": "GENOMIC_PROXIMITY_CLUSTER",
                "priority": "MEDIUM",
                "sample_id": sample_id,
                "cluster_id": cluster_id,
                "resistant_drugs": "",
                "classification": "",
                "message": "该样本与同簇样本基因组距离达到阈值；建议结合采样时间、地点和接触史开展流调核查，不代表已证实传播。",
            })

        resistance = resistance_rows.get(sample_id)
        classification = resistance.get("classification", "") if resistance else ""
        organism_scope = resistance.get("organism_scope", "") if resistance else ""
        # Older per-sample files have no scope column; retain compatibility for
        # those files. If a scope is present, only the bundled MTB/TB screen is
        # eligible for TB-specific alert wording.
        if (
            resistance
            and classification not in {"", "NO_KNOWN_MARKER"}
            and organism_scope in {"", RESULT_SCOPE}
        ):
            priority = (
                "CRITICAL"
                if classification == "XDR"
                else "HIGH"
                if classification in {"RR", "MDR", "pre-XDR"}
                else "MEDIUM"
            )
            alerts.append({
                "alert_type": "DRUG_RESISTANCE_RISK",
                "priority": priority,
                "sample_id": sample_id,
                "cluster_id": cluster_id,
                "resistant_drugs": resistance.get("resistant_drugs", ""),
                "classification": classification,
                "organism_scope": organism_scope or RESULT_SCOPE,
                "message": "检出 MTB/TB 耐药相关基因组标记；请由结核病专业人员结合表型药敏、既往用药史和临床资料复核，不应由本提示单独调整治疗。",
            })

    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "alert_type", "priority", "sample_id", "cluster_id", "resistant_drugs",
        "classification", "organism_scope", "message",
    ]
    with open(output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(alerts)
    logger.info("Wrote %d investigation alerts to %s", len(alerts), output)
    return output


def investigation_alerts_cmd(
    cluster_assignments: Path = Option(..., "--cluster-assignments", help="Cluster assignment TSV"),
    drug_resistance_dir: Path = Option(..., "--drug-resistance-dir", help="Results/samples directory"),
    output: Path = Option(..., "--output", "-o", help="Investigation alert TSV"),
) -> None:
    """Generate explainable epidemiological follow-up prompts."""
    typer.echo(f"Investigation alerts written to {generate_alerts(cluster_assignments, drug_resistance_dir, output)}")
