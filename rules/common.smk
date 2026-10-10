import csv
from pathlib import Path

from snakemake.exceptions import WorkflowError


RESULTS_DIR = config.get("results_dir", "results")
SAMPLESHEET = Path(config.get("samplesheet", "config/samplesheet.csv"))
PYTHON_BIN = config.get("python_bin", "python")
CLI_MODULE = config.get("cli_module", "src.cli")
THREADS_CFG = config.get("threads", {})
REPORT_CFG = config.get("report", {})
WORKFLOW_CFG = config.get("workflow", {})
REFERENCE_CFG = config.get("references", {})
CLUSTERING_CFG = config.get("clustering", {})
DRUG_RESISTANCE_CFG = config.get("drug_resistance", {})

ANALYSIS_ROUTE = str(config.get("analysis_route", "TB_bowtie2"))
if ANALYSIS_ROUTE not in {"TB_bowtie2", "NTM_bwamem"}:
    raise WorkflowError("analysis_route must be TB_bowtie2 or NTM_bwamem")
MAB_ROUTE = ANALYSIS_ROUTE == "NTM_bwamem"
BWA_THREADS = int(THREADS_CFG.get("bwa", 8))
ALIGN_THREADS = BWA_THREADS if MAB_ROUTE else int(THREADS_CFG.get("bowtie2", 4))
TRIM_MIN_LENGTH_ARG = "-l 35" if MAB_ROUTE else ""

TRIM_READS_THREADS = int(THREADS_CFG.get("trim_reads", 1))
BOWTIE2_THREADS = int(THREADS_CFG.get("bowtie2", 4))
SAMTOOLS_SORT_THREADS = int(THREADS_CFG.get("samtools_sort", 2))
DEPTH_METRICS_THREADS = int(THREADS_CFG.get("depth_metrics", 1))
CALL_VARIANTS_THREADS = int(THREADS_CFG.get("call_variants", 1))
DIFF_LOCI_THREADS = int(THREADS_CFG.get("diff_loci", 1))
RECALL_THREADS = int(THREADS_CFG.get("recall", 1))
MERGE_THREADS = int(THREADS_CFG.get("merge", 1))
WILD_EXTRACT_THREADS = int(THREADS_CFG.get("wild_extract", 1))
FILTER_CORE_THREADS = int(THREADS_CFG.get("filter_core", 1))
DISTANCE_THREADS = int(THREADS_CFG.get("distance", 1))
CLUSTER_THREADS = int(THREADS_CFG.get("cluster", 1))
DRUG_RESISTANCE_THREADS = int(THREADS_CFG.get("drug_resistance", 1))
INVESTIGATION_ALERTS_THREADS = int(THREADS_CFG.get("investigation_alerts", 1))
CLUSTER_THRESHOLD = int(CLUSTERING_CFG.get("snp_threshold", 12))
DRUG_DB = str(Path(DRUG_RESISTANCE_CFG.get("database", "data/drug_db/Total_drug_resistance_mutation.txt")))
VERBOSE = bool(WORKFLOW_CFG.get("verbose", False))
RUN_ID = str(config.get("run_id", WORKFLOW_CFG.get("run_id", "latest"))).strip()
if not RUN_ID:
    RUN_ID = "latest"
STRICT = bool(REPORT_CFG.get("strict", False))
SKIP_LINEAGE = bool(REPORT_CFG.get("skip_lineage", False))
SKIP_FIGURES = bool(REPORT_CFG.get("skip_figures", False))
SKIP_TABLES = bool(REPORT_CFG.get("skip_tables", False))
SKIP_EXTENSIONS = bool(REPORT_CFG.get("skip_extensions", False))

REFERENCE_PROFILE = str(REFERENCE_CFG.get("active", "tb_ancestor"))
REFERENCE_PROFILES = REFERENCE_CFG.get("profiles", {})
ACTIVE_REFERENCE_CFG = (
    REFERENCE_PROFILES.get(REFERENCE_PROFILE, {})
    if isinstance(REFERENCE_PROFILES, dict)
    else {}
)

# New profile-style config (references.active + references.profiles.*) is preferred.
# Keep backward compatibility with old keys (h37rv/h37rv_fai/bowtie2_index/ancestor_fasta).
REFERENCE_FASTA = str(
    Path(
        ACTIVE_REFERENCE_CFG.get(
            "fasta",
            REFERENCE_CFG.get("h37rv", "data/tb.ancestor.fasta"),
        )
    )
)
REFERENCE_FASTA_FAI = str(
    Path(
        ACTIVE_REFERENCE_CFG.get(
            "fai",
            REFERENCE_CFG.get("h37rv_fai", f"{REFERENCE_FASTA}.fai"),
        )
    )
)
BOWTIE2_INDEX_PREFIX = str(
    Path(
        ACTIVE_REFERENCE_CFG.get(
            "bowtie2_index",
            REFERENCE_CFG.get("bowtie2_index", "data/bowtie2_index/tb.ancestor.fasta"),
        )
    )
)
ANCESTOR_FASTA = str(
    Path(
        ACTIVE_REFERENCE_CFG.get(
            "wildtype_fasta",
            REFERENCE_CFG.get("ancestor_fasta", REFERENCE_FASTA),
        )
    )
)
PPE_LIST = str(Path(REFERENCE_CFG.get("ppe_list", "data/PPE_INS_loci_Rv.list")))
BWA_INDEX_PREFIX = str(ACTIVE_REFERENCE_CFG.get("bwa_index", REFERENCE_FASTA))
if MAB_ROUTE and not ACTIVE_REFERENCE_CFG.get("ppe_list"):
    raise WorkflowError("NTM_bwamem requires references.profiles.<active>.ppe_list matched to that reference")
PPE_LIST = str(ACTIVE_REFERENCE_CFG.get("ppe_list", PPE_LIST))
VARSCAN_JAR = str(Path(REFERENCE_CFG.get("varscan_jar", "src/scripts/VarScan.v2.3.9.jar")))
GENOME_LENGTH = int(ACTIVE_REFERENCE_CFG.get("genome_length", REFERENCE_CFG.get("genome_length", 4411532)))

if not SAMPLESHEET.exists():
    raise WorkflowError(f"Samplesheet not found: {SAMPLESHEET}")


def load_samplesheet(path: Path):
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"sample_id", "r1", "r2"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise WorkflowError(
                f"Samplesheet must contain columns {sorted(required)}: {path}"
            )
        for row in reader:
            sample_id = str(row.get("sample_id", "")).strip()
            r1 = str(row.get("r1", "")).strip()
            r2 = str(row.get("r2", "")).strip()
            if not sample_id or not r1 or not r2:
                raise WorkflowError(f"Invalid row in samplesheet {path}: {row}")
            rows.append({"sample_id": sample_id, "r1": r1, "r2": r2})
    if not rows:
        raise WorkflowError(f"Samplesheet is empty: {path}")
    return rows


SAMPLE_ROWS = load_samplesheet(SAMPLESHEET)
SAMPLES = [row["sample_id"] for row in SAMPLE_ROWS]
SAMPLE_ROW_MAP = {row["sample_id"]: row for row in SAMPLE_ROWS}


if len(SAMPLES) != len(set(SAMPLES)):
    raise WorkflowError(f"Duplicate sample_id detected in {SAMPLESHEET}: {SAMPLES}")


def result_path(*parts):
    return str(Path(RESULTS_DIR).joinpath(*parts))


def sample_r1(sample_id: str) -> str:
    return SAMPLE_ROW_MAP[sample_id]["r1"]


def sample_r2(sample_id: str) -> str:
    return SAMPLE_ROW_MAP[sample_id]["r2"]


def per_sample(path_parts_fn):
    return [path_parts_fn(sample) for sample in SAMPLES]


def run_log_path(*parts):
    return str(Path(RESULTS_DIR).joinpath("logs", "runs", RUN_ID, *parts))


FOUNDATION_FASTP_JSONS = [
    result_path("samples", sample, "report_inputs", "fastp_qc", f"{sample}_{mate}_fastp.json")
    for sample in SAMPLES
    for mate in (1, 2)
]

ALIGNMENT_OUTPUTS = [
    result_path("samples", sample, "alignment_qc", f"{sample}.sort.bam")
    for sample in SAMPLES
]

VARIANT_OUTPUTS = [
    result_path("samples", sample, "variant_analysis", f"{sample}.{suffix}")
    for sample in SAMPLES
    for suffix in ("cns", "snp", "vars")
]

CORE_OUTPUTS = [
    result_path("core", name)
    for name in (
        "diff_loci.txt",
        "merged.fasta",
        "wildtype.fasta",
        "core_snps.fadel-InvMisF5.bak.fa",
        "core_snps.fadel-InvMisF5.bak.loc",
        "distance_matrix.txt",
    )
]

DOWNSTREAM_INPUTS = [
    result_path("samples", sample, "report_inputs", "alignment_qc", f"{sample}.depth")
    for sample in SAMPLES
] + [
    result_path("samples", sample, "report_inputs", "alignment_qc", "insert_sizes.txt")
    for sample in SAMPLES
] + [
    result_path("samples", sample, "report_inputs", "annotated_variants", f"{sample}_annotated.txt")
    for sample in SAMPLES
]

FIGURE_STATUS_OUTPUTS = [] if SKIP_FIGURES else [
    run_log_path("report_figures_status.tsv"),
    result_path("logs", "report_figures_status.tsv"),
]

TABLE_STATUS_OUTPUTS = [] if SKIP_TABLES else [
    run_log_path("report_tables_status.tsv"),
    result_path("logs", "report_tables_status.tsv"),
]

LINEAGE_SUMMARY_OUTPUTS = [] if SKIP_LINEAGE else [
    run_log_path("lineage_summary.tsv"),
    result_path("logs", "lineage_summary.tsv"),
]

FOUNDATION_STATUS_RUN = run_log_path("prepare_foundation_outputs_status.tsv")
FOUNDATION_STATUS_LATEST = result_path("logs", "prepare_foundation_outputs_status.tsv")
DOWNSTREAM_STATUS_RUN = run_log_path("prepare_downstream_inputs_status.tsv")
DOWNSTREAM_STATUS_LATEST = result_path("logs", "prepare_downstream_inputs_status.tsv")

LINEAGE_SUMMARY_RUN = run_log_path("lineage_summary.tsv")
LINEAGE_SUMMARY_LATEST = result_path("logs", "lineage_summary.tsv")
FIGURE_STATUS_RUN = run_log_path("report_figures_status.tsv")
FIGURE_STATUS_LATEST = result_path("logs", "report_figures_status.tsv")
TABLE_STATUS_RUN = run_log_path("report_tables_status.tsv")
TABLE_STATUS_LATEST = result_path("logs", "report_tables_status.tsv")

TABLE_OUTPUTS = [] if SKIP_TABLES else [
    result_path("samples", sample, "table", f"临床测序质量报告_{sample}.csv")
    for sample in SAMPLES
] + [
    result_path("samples", sample, "table", f"临床变异检测报告_{sample}.csv")
    for sample in SAMPLES
] + [
    result_path("samples", sample, "table", f"样本级变异信息表_{sample}.csv")
    for sample in SAMPLES
] + [
    result_path("samples", sample, "table", f"谱系鉴定结果_{sample}.csv")
    for sample in SAMPLES
]

WORKFLOW_STATE_DIR = result_path("logs", "workflow")
RULE_LOG_DIR = result_path("logs", "rules")


VERBOSE_FLAG = "--verbose" if VERBOSE else ""
STRICT_FLAG = "--strict" if STRICT else ""
SKIP_LINEAGE_FLAG = "--skip-lineage" if SKIP_LINEAGE else ""
SKIP_FIGURES_FLAG = "--skip-figures" if SKIP_FIGURES else ""
SKIP_TABLES_FLAG = "--skip-tables" if SKIP_TABLES else ""
