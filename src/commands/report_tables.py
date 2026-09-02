"""Generate per-sample report tables."""

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

import typer
from typer import Option

from src.core.report_inputs import annotated_variant_candidates
from src.reporting.variant_classification import (
    summarize_variant_classes,
    write_normalized_annotated_from_cns,
)
from src.utils.logging_config import get_logger

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LEGACY_BASE = PROJECT_ROOT / "data" / "report_assets"
DRUG_DB_PATH = PROJECT_ROOT / "data" / "drug_db" / "Total_drug_resistance_mutation.txt"
LEGACY_ALIAS = "N24_1"
logger = get_logger("commands.report_tables")


def _find_first(paths: List[Path]) -> Optional[Path]:
    for p in paths:
        if p.exists():
            return p
    return None


def _collect_samples(input_root: Path) -> List[str]:
    if input_root.name == "samples":
        candidates = sorted(input_root.glob("*/variant_analysis/*.cns"))
    else:
        candidates = sorted((input_root / "samples").glob("*/variant_analysis/*.cns"))

    if candidates:
        return sorted({p.stem for p in candidates})

    if input_root.name == "samples":
        snps = sorted(input_root.glob("*/variant_analysis/*.snp"))
    else:
        snps = sorted((input_root / "samples").glob("*/variant_analysis/*.snp"))
    if snps:
        return sorted({p.stem for p in snps})

    return []

def _sample_dir_name(sample_id: str) -> str:
    return sample_id if sample_id.endswith(".cleaned") else f"{sample_id}.cleaned"


def _sample_input_roots(input_root: Path, sample_id: str) -> List[Path]:
    names = [sample_id]
    canonical = _sample_dir_name(sample_id)
    if canonical != sample_id:
        names.append(canonical)

    roots: List[Path] = []
    for name in names:
        if input_root.name == "samples":
            roots.append(input_root / name)
        else:
            roots.append(input_root / "samples" / name)
    return roots


def _sample_output_dir(output_root: Path, sample_id: str) -> Path:
    name = _sample_dir_name(sample_id)
    return output_root / name if output_root.name == "samples" else output_root / "samples" / name


def _global_logs_dir(input_root: Path, output_root: Path) -> Path:
    if (input_root / "logs").exists():
        return input_root / "logs"
    if output_root.name == "samples":
        return output_root.parent / "logs"
    return output_root / "logs"


def _lineage_summary_candidates(input_root: Path, run_id: str) -> List[Path]:
    results_root = input_root.parent if input_root.name == "samples" else input_root
    candidates: List[Path] = []
    if run_id:
        candidates.append(results_root / "logs" / "runs" / run_id / "lineage_summary.tsv")
    candidates.append(results_root / "logs" / "lineage_summary.tsv")
    candidates.append(results_root / "logs" / "runs" / "latest" / "lineage_summary.tsv")
    return candidates


def _results_root(input_root: Path) -> Path:
    return input_root.parent if input_root.name == "samples" else input_root


def _extension_record(input_root: Path, sample_id: str, sample_roots: List[Path]) -> Dict[str, str]:
    """Load cluster and resistance fields without making them hard dependencies."""
    results_root = _results_root(input_root)
    record = {
        "cluster_id": "UNKNOWN",
        "cluster_size": "",
        "drug_resistance_class": "NOT_RUN",
        "resistant_drugs": "-",
        "organism_scope": "UNSPECIFIED",
    }

    cluster_file = results_root / "core" / "cluster_assignments.tsv"
    if cluster_file.exists():
        with open(cluster_file, "r", encoding="utf-8", errors="ignore") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if str(row.get("sample_id", "")) == sample_id:
                    record["cluster_id"] = str(row.get("cluster_id", "UNKNOWN"))
                    record["cluster_size"] = str(row.get("cluster_size", ""))
                    break

    resistance_candidates = [
        root / "table" / f"{sample_id}_drug_resistance.tsv" for root in sample_roots
    ]
    resistance_file = _find_first(resistance_candidates)
    if resistance_file and resistance_file.exists():
        with open(resistance_file, "r", encoding="utf-8", errors="ignore") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"), None)
        if row:
            record["drug_resistance_class"] = str(row.get("classification", record["drug_resistance_class"]))
            record["resistant_drugs"] = str(row.get("resistant_drugs", record["resistant_drugs"]))
            record["organism_scope"] = str(row.get("organism_scope", record["organism_scope"]))
    return record



def _resolve_legacy_base(user_base: Optional[Path] = None) -> Path:
    candidates: List[Path] = []
    if user_base is not None:
        candidates.append(user_base)

    env_base = os.environ.get("MTB_EVO_REPORT_TEMPLATE_DIR")
    if env_base:
        candidates.append(Path(env_base))

    candidates.append(DEFAULT_LEGACY_BASE)

    for base in candidates:
        if (base / "scripts" / "generate_clinical_variant_report.py").exists() and (
            base / "scripts" / "generate_variant_table_final.py"
        ).exists():
            return base

    raise FileNotFoundError(
        f"Legacy table scripts not found under {DEFAULT_LEGACY_BASE}. "
        "Set --legacy-base-dir or MTB_EVO_REPORT_TEMPLATE_DIR to override."
    )


def _safe_link_or_copy(src: Path, dst: Path) -> None:
    src = src.resolve()
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        os.symlink(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def _build_clinical_variant_table_from_cns(cns_file: Path, annotated_file: Path, out_csv: Path) -> None:
    summary = summarize_variant_classes(cns_file, annotated_file)
    total = summary.total
    snp_base = summary.snp if summary.snp > 0 else 1
    total_base = total if total > 0 else 1
    rows = [
        ["总变异数", total, "100%", "-", "/"],
        ["同义突变", summary.syn, f"{summary.syn / snp_base * 100:.1f}%", "不改变氨基酸 通常无临床意义", "低"],
        ["非同义突变", summary.nonsyn, f"{summary.nonsyn / snp_base * 100:.1f}%", "改变氨基酸 可能影响蛋白功能", "高"],
        ["基因间区变异", summary.intergenic, f"{summary.intergenic / snp_base * 100:.1f}%", "非编码区 通常无临床意义", "低"],
        ["点突变(SNP)", summary.snp, f"{summary.snp / total_base * 100:.1f}%", "常见变异类型", "中"],
        ["插入缺失(Indel)", summary.indel, f"{summary.indel / total_base * 100:.1f}%", "可能影响较大", "高"],
    ]

    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["变异类型", "数量", "占比", "临床意义", "建议关注程度"])
        writer.writerows(rows)


def _build_lineage_table(
    input_root: Path,
    sample_id: str,
    out_csv: Path,
    run_id: str = "",
) -> Dict[str, str]:
    summary_lineage = _find_first(_lineage_summary_candidates(input_root, run_id))

    record = {
        "sample_id": sample_id,
        "lineage_major": "UNKNOWN",
        "lineage_sub": "",
        "method": "",
        "marker_count_total": "",
        "marker_count_hit": "",
        "status": "NO_LINEAGE_FILE",
        "message": "Lineage output not found",
    }
    record.update(_extension_record(input_root, sample_id, _sample_input_roots(input_root, sample_id)))

    def _load_first_row(tsv_path: Path) -> Dict[str, str] | None:
        with open(tsv_path, "r", encoding="utf-8", errors="ignore") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                return {k: str(v) for k, v in row.items()}
        return None

    if summary_lineage and summary_lineage.exists():
        with open(summary_lineage, "r", encoding="utf-8", errors="ignore") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if str(row.get("sample_id", "")) == sample_id:
                    record.update({k: str(v) for k, v in row.items()})
                    break

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["项目", "结果", "说明"])
        writer.writerow(["样本ID", record.get("sample_id", sample_id), ""])
        writer.writerow(["大谱系", record.get("lineage_major", "UNKNOWN"), "主要谱系判定"])
        writer.writerow(["亚谱系", record.get("lineage_sub", ""), "细分亚谱系判定"])
        writer.writerow(["基因组近缘簇", record.get("cluster_id", "UNKNOWN"), "基于SNP阈值的连通分量；不等同于确证传播"])
        writer.writerow(["簇内样本数", record.get("cluster_size", ""), "需结合流调资料解释"])
        writer.writerow(["耐药风险分级", record.get("drug_resistance_class", "NOT_RUN"), "基因组标记筛查，不替代表型药敏"])
        writer.writerow(["耐药相关药物", record.get("resistant_drugs", "-"), "数据库命中药物"])
        writer.writerow(["耐药筛查适用范围", record.get("organism_scope", "UNSPECIFIED"), "TB/MTB 数据专用；不能直接解释 NTM"])
        writer.writerow(["方法", record.get("method", ""), "谱系规则版本"])
        writer.writerow(
            [
                "命中标记数",
                f"{record.get('marker_count_hit', '')}/{record.get('marker_count_total', '')}",
                "命中marker/总marker",
            ]
        )
        writer.writerow(["状态", record.get("status", ""), "OK 表示判定成功"])
        writer.writerow(["备注", record.get("message", ""), "附加解释信息"])
    return record


def _patch_base_dir(src_script: Path, dst_script: Path, base_dir: Path) -> None:
    txt = src_script.read_text(encoding="utf-8")
    txt = re.sub(r'^BASE_DIR\s*=\s*.*$', f'BASE_DIR = r"{base_dir}"', txt, flags=re.M)
    txt = re.sub(r'^DRUG_DB\s*=\s*.*$', f'DRUG_DB = r"{DRUG_DB_PATH}"', txt, flags=re.M)
    dst_script.write_text(txt, encoding="utf-8")


def _parse_percent(s: str) -> float:
    s = s.strip().replace("%", "")
    try:
        return float(s)
    except ValueError:
        return 0.0


def _judge(value: float, rule: str) -> str:
    if rule.startswith(">"):
        target = float(rule[1:])
        return "合格" if value > target else "不合格"
    if rule.startswith("<"):
        target = float(rule[1:])
        return "优秀" if value < target else "不合格"
    return "/"


def _format_int(n: int) -> str:
    return f"{n:,}"


def _parse_leading_count(line: str) -> int:
    parts = line.strip().split()
    if not parts:
        return 0
    try:
        return int(parts[0])
    except ValueError:
        return 0


def _run_samtools_flagstat(bam_file: Path) -> Dict[str, int]:
    out: Dict[str, int] = {
        "total_reads": 0,
        "mapped_reads": 0,
        "secondary_reads": 0,
        "supplementary_reads": 0,
        "duplicate_reads": 0,
    }
    if not bam_file.exists() or shutil.which("samtools") is None:
        return out
    try:
        p = subprocess.run(
            ["samtools", "flagstat", str(bam_file)],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return out

    for line in p.stdout.splitlines():
        if "in total" in line:
            out["total_reads"] = _parse_leading_count(line)
        elif " secondary" in line:
            out["secondary_reads"] = _parse_leading_count(line)
        elif " supplementary" in line:
            out["supplementary_reads"] = _parse_leading_count(line)
        elif "duplicates" in line and "primary duplicates" not in line:
            out["duplicate_reads"] = _parse_leading_count(line)
        elif " mapped (" in line and "primary mapped" not in line:
            out["mapped_reads"] = _parse_leading_count(line)
    return out


def _run_samtools_error_rate(bam_file: Path) -> float:
    if not bam_file.exists() or shutil.which("samtools") is None:
        return 0.0
    try:
        p = subprocess.run(
            ["samtools", "stats", str(bam_file)],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return 0.0

    for line in p.stdout.splitlines():
        if line.startswith("SN") and "error rate:" in line:
            fields = line.split("\t")
            if len(fields) >= 3:
                try:
                    return float(fields[2]) * 100.0
                except ValueError:
                    return 0.0
    return 0.0


def _build_quality_table(input_root: Path, sample_id: str, out_csv: Path) -> None:
    sample_roots = _sample_input_roots(input_root, sample_id)
    bam_stat = _find_first([
        *[root / "report_inputs" / "alignment_qc" / "bam_statistic.csv" for root in sample_roots],
    ])
    depth_file = _find_first([
        *[root / "report_inputs" / "alignment_qc" / f"{sample_id}.depth" for root in sample_roots],
    ])
    bam_stat_out = _find_first([
        *[root / "alignment_qc" / f"{sample_id}_bam_stat.out" for root in sample_roots],
    ])
    dedup_metrics = _find_first([
        *[root / "alignment_qc" / "dedup_metrics.txt" for root in sample_roots],
    ])
    samtools_stats = _find_first([
        *[root / "alignment_qc" / "samtools_stats.txt" for root in sample_roots],
    ])
    bam_file = _find_first([
        *[root / "alignment_qc" / f"{sample_id}.sort.bam" for root in sample_roots],
    ])

    total_reads = mapped_reads = secondary_reads = supplementary_reads = 0
    if bam_stat_out and bam_stat_out.exists():
        for line in bam_stat_out.read_text(encoding="utf-8", errors="ignore").splitlines():
            if "in total" in line:
                total_reads = _parse_leading_count(line)
            elif " secondary" in line:
                secondary_reads = _parse_leading_count(line)
            elif " supplementary" in line:
                supplementary_reads = _parse_leading_count(line)
            elif " mapped (" in line and "primary mapped" not in line:
                mapped_reads = _parse_leading_count(line)

    # Fallback: compute from BAM when *_bam_stat.out is absent or parsing was incomplete.
    duplicate_reads = 0
    if bam_file and (total_reads == 0 or mapped_reads == 0):
        stat = _run_samtools_flagstat(bam_file)
        if total_reads == 0:
            total_reads = stat["total_reads"]
        if mapped_reads == 0:
            mapped_reads = stat["mapped_reads"]
        if secondary_reads == 0:
            secondary_reads = stat["secondary_reads"]
        if supplementary_reads == 0:
            supplementary_reads = stat["supplementary_reads"]
        duplicate_reads = stat["duplicate_reads"]

    unique_mapped = max(0, mapped_reads - secondary_reads - supplementary_reads)
    map_rate = (mapped_reads / total_reads * 100.0) if total_reads else 0.0
    unique_rate = (unique_mapped / total_reads * 100.0) if total_reads else 0.0

    mismatch_rate = 0.0
    if samtools_stats and samtools_stats.exists():
        for line in samtools_stats.read_text(encoding="utf-8", errors="ignore").splitlines():
            if line.startswith("SN") and "error rate:" in line:
                fields = line.split("\t")
                if len(fields) >= 3:
                    try:
                        mismatch_rate = float(fields[2]) * 100.0
                    except ValueError:
                        mismatch_rate = 0.0
    elif bam_file:
        mismatch_rate = _run_samtools_error_rate(bam_file)

    dup_rate = 0.0
    if dedup_metrics and dedup_metrics.exists():
        lines = dedup_metrics.read_text(encoding="utf-8", errors="ignore").splitlines()
        if len(lines) >= 2:
            fields = lines[1].split("\t")
            if len(fields) >= 8:
                try:
                    dup_rate = float(fields[7]) * 100.0
                except ValueError:
                    dup_rate = 0.0
    elif total_reads > 0 and duplicate_reads > 0:
        dup_rate = duplicate_reads / total_reads * 100.0

    mean_depth = 0.0
    cov1 = cov10 = cov30 = 0.0
    if depth_file and depth_file.exists():
        depth_vals: List[int] = []
        with open(depth_file, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) < 3:
                    continue
                try:
                    d = int(cols[2])
                except ValueError:
                    continue
                depth_vals.append(d)
        if depth_vals:
            total_pos = len(depth_vals)
            mean_depth = sum(depth_vals) / total_pos
            cov1 = sum(1 for d in depth_vals if d >= 1) / total_pos * 100.0
            cov10 = sum(1 for d in depth_vals if d >= 10) / total_pos * 100.0
            cov30 = sum(1 for d in depth_vals if d >= 30) / total_pos * 100.0

    if mean_depth == 0.0 and bam_stat and bam_stat.exists():
        with open(bam_stat, "r", encoding="utf-8", errors="ignore") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("Sample") == sample_id:
                    try:
                        mean_depth = float(str(row.get("Averaged_Depth", "0")).replace("X", ""))
                    except ValueError:
                        mean_depth = 0.0
                    if not cov1:
                        cov1 = _parse_percent(str(row.get("X1x", "0")))
                    if not cov10:
                        cov10 = _parse_percent(str(row.get("X10x", "0")))
                    if not cov30:
                        cov30 = _parse_percent(str(row.get("X30x", "0")))
                    break

    rows = [
        ["总测序读数", _format_int(total_reads), "-", "/", "测序数据总量"],
        ["成功比对读数", f"{_format_int(mapped_reads)} ({map_rate:.2f}%)", ">95%", _judge(map_rate, ">95"), "与参考基因组匹配的比例"],
        ["比对率", f"{map_rate:.2f}%", ">95%", _judge(map_rate, ">95"), "数据质量良好"],
        ["唯一比对读数", f"{_format_int(unique_mapped)} ({unique_rate:.2f}%)", ">95%", _judge(unique_rate, ">95"), "非重复的有效读数"],
        ["唯一比对率", f"{unique_rate:.2f}%", ">95%", _judge(unique_rate, ">95"), "数据可靠性高"],
        ["碱基错配率", f"{mismatch_rate:.2f}%", "<1%", _judge(mismatch_rate, "<1"), "测序准确性高"],
        ["重复读数率", f"{dup_rate:.2f}%", "<20%", _judge(dup_rate, "<20"), "无PCR重复污染"],
        ["平均测序深度", f"{mean_depth:.2f}X", ">100X", _judge(mean_depth, ">100"), "覆盖度充足"],
        ["1X覆盖度", f"{cov1:.2f}%", ">99%", _judge(cov1, ">99"), "基因组覆盖完整"],
        ["10X覆盖度", f"{cov10:.2f}%", ">95%", _judge(cov10, ">95"), "高质量覆盖区域"],
        ["30X覆盖度", f"{cov30:.2f}%", ">95%", _judge(cov30, ">95"), "可用于变异检测"],
    ]

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["指标名称", "数值", "参考标准", "结果判定", "临床意义"])
        writer.writerows(rows)


def run_report_tables(
    input_root: Path,
    output_root: Path,
    status_output: Path,
    pipeline_version: str = "v1",
    legacy_base_dir: Optional[Path] = None,
    run_id: str = "",
    stage: str = "report_tables",
) -> Path:
    _ = pipeline_version  # keep signature stable for pipeline calls

    logger.info("Running report tables with input_root=%s output_root=%s", input_root, output_root)
    samples = _collect_samples(input_root)
    if not samples:
        raise ValueError(f"No sample inferred from {input_root} (*.cns/*.snp not found)")

    legacy_base = _resolve_legacy_base(legacy_base_dir)
    script_clin = legacy_base / "scripts" / "generate_clinical_variant_report.py"
    script_detail = legacy_base / "scripts" / "generate_variant_table_final.py"

    status_rows: List[Dict[str, str]] = []
    logger.info("Discovered %d samples for table generation", len(samples))

    for sample_id in samples:
        sample_out = _sample_output_dir(output_root, sample_id)
        table_dir = sample_out / "table"
        log_dir = _global_logs_dir(input_root, output_root)
        table_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)
        sample_roots = _sample_input_roots(input_root, sample_id)

        cns = _find_first(
            [root / "variant_analysis" / f"{sample_id}.cns" for root in sample_roots]
        )
        annotated = _find_first(annotated_variant_candidates(sample_roots, sample_id))

        logger.info("Table sample start: %s", sample_id)
        issues: List[str] = []
        if not cns:
            issues.append("MISSING_CNS")
        if not annotated:
            issues.append("MISSING_ANNOTATED")

        run_logs: List[str] = []

        if not issues:
            with tempfile.TemporaryDirectory(prefix=f"mtb-evo-tab-{sample_id}-") as tmpdir:
                work = Path(tmpdir)
                for d in ["scripts", "table", "variant_analysis"]:
                    (work / d).mkdir(parents=True, exist_ok=True)

                _safe_link_or_copy(cns, work / "variant_analysis" / f"{LEGACY_ALIAS}.cns")
                write_normalized_annotated_from_cns(
                    cns, annotated, work / "variant_analysis" / f"{LEGACY_ALIAS}_annotated.txt"
                )

                clin_patched = work / "scripts" / "generate_clinical_variant_report.py"
                detail_patched = work / "scripts" / "generate_variant_table_final.py"
                _patch_base_dir(script_clin, clin_patched, work)
                _patch_base_dir(script_detail, detail_patched, work)

                for cmd, name in [
                    (["python3", str(clin_patched)], "generate_clinical_variant_report.py"),
                    (["python3", str(detail_patched)], "generate_variant_table_final.py"),
                ]:
                    try:
                        p = subprocess.run(cmd, cwd=work, check=True, capture_output=True, text=True)
                        run_logs.append(f"[{name}]\n{(p.stdout or '') + (p.stderr or '')}")
                    except subprocess.CalledProcessError as e:
                        issues.append(f"SCRIPT_FAILED:{name}")
                        run_logs.append(f"[{name}]\n{(e.stdout or '') + (e.stderr or '')}")

                # copy legacy-script tables with sample-id rename
                for src_name, dst_name in [
                    (f"临床变异检测报告_{LEGACY_ALIAS}.csv", f"临床变异检测报告_{sample_id}.csv"),
                    (f"样本级变异信息表_{LEGACY_ALIAS}.csv", f"样本级变异信息表_{sample_id}.csv"),
                ]:
                    src = work / "table" / src_name
                    if src.exists():
                        shutil.copy2(src, table_dir / dst_name)
                    else:
                        issues.append(f"MISSING_OUTPUT:{dst_name}")

        quality_csv = table_dir / f"临床测序质量报告_{sample_id}.csv"
        try:
            _build_quality_table(input_root, sample_id, quality_csv)
        except Exception as e:
            issues.append(f"QUALITY_REPORT_FAILED:{e}")

        clinical_csv = table_dir / f"临床变异检测报告_{sample_id}.csv"
        try:
            _build_clinical_variant_table_from_cns(cns, annotated, clinical_csv)
        except Exception as e:
            issues.append(f"CLINICAL_REPORT_FAILED:{e}")

        lineage_csv = table_dir / f"谱系鉴定结果_{sample_id}.csv"
        try:
            lineage_record = _build_lineage_table(input_root, sample_id, lineage_csv, run_id=run_id)
            if lineage_record.get("status") != "OK":
                issues.append(f"LINEAGE_STATUS:{lineage_record.get('status', 'UNKNOWN')}")
        except Exception as e:
            issues.append(f"LINEAGE_TABLE_FAILED:{e}")

        realized_outputs = {
            f"临床变异检测报告_{sample_id}.csv": clinical_csv.exists(),
            f"样本级变异信息表_{sample_id}.csv": (table_dir / f"样本级变异信息表_{sample_id}.csv").exists(),
            f"临床测序质量报告_{sample_id}.csv": quality_csv.exists(),
            f"谱系鉴定结果_{sample_id}.csv": lineage_csv.exists(),
        }
        filtered_issues: List[str] = []
        for issue in issues:
            if issue.startswith("MISSING_OUTPUT:"):
                name = issue.split(":", 1)[1]
                if realized_outputs.get(name, False):
                    continue
            filtered_issues.append(issue)
        issues = filtered_issues

        status = "OK" if not issues else "PARTIAL_OK"
        with open(log_dir / f"report_tables_{sample_id}.log", "w", encoding="utf-8") as f:
            f.write(f"sample_id={sample_id}\n")
            f.write(f"status={status}\n")
            f.write(f"issues={';'.join(issues)}\n")
            f.write(f"legacy_base={legacy_base}\n")
            for chunk in run_logs:
                f.write(chunk)
                f.write("\n")

        logger.info("Table sample done: %s status=%s issues=%s", sample_id, status, ";".join(issues))
        status_rows.append(
            {
                "sample_id": sample_id,
                "status": status,
                "issues": ";".join(issues),
                "table_dir": str(table_dir),
                "run_id": run_id,
                "stage": stage,
            }
        )

    status_output.parent.mkdir(parents=True, exist_ok=True)
    with open(status_output, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["sample_id", "status", "issues", "table_dir", "run_id", "stage"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(status_rows)

    logger.info("Table status summary written to %s", status_output)
    return status_output


def report_tables_cmd(
    input_root: Path = Option(Path("results"), "--input-root", help="Input root with pipeline outputs"),
    output_root: Path = Option(Path("results/samples"), "--output-root", help="Output root for per-sample reports"),
    status_output: Path = Option(
        Path("results/logs/report_tables_status.tsv"), "--status-output", help="Status TSV output"
    ),
    pipeline_version: str = Option("v1", "--pipeline-version", help="Pipeline version in output metadata"),
    legacy_base_dir: Optional[Path] = Option(
        None,
        "--legacy-base-dir",
        help="Legacy report template directory (contains scripts/generate_* files)",
    ),
    run_id: str = Option("", "--run-id", help="Run identifier stored in status metadata"),
) -> None:
    """Generate table reports by calling original legacy scripts."""
    out = run_report_tables(
        input_root,
        output_root,
        status_output,
        pipeline_version=pipeline_version,
        legacy_base_dir=legacy_base_dir,
        run_id=run_id,
        stage="report_tables",
    )
    typer.echo(f"Table status summary generated: {out}")
