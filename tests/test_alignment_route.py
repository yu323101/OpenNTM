from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_all_routes_align_with_bwa_mem() -> None:
    core_rule = (PROJECT_ROOT / "rules/core.smk").read_text()
    common_rule = (PROJECT_ROOT / "rules/common.smk").read_text()

    assert "NTM_bwamem.sh align" in core_rule
    assert "bowtie2" not in core_rule
    assert "ALIGN_THREADS = BWA_THREADS" in common_rule
    assert '"TB_bwamem", "TB_bowtie2", "NTM_bwamem"' in common_rule


def test_reference_profiles_declare_bwa_index() -> None:
    config = (PROJECT_ROOT / "config/config.yaml").read_text()

    assert "bwa_index: data/bwa_index/tb.ancestor" in config
    assert "bwa_index: data/bwa_index/tb_h37rv" in config
    assert "bowtie2_index" not in config


def test_bwa_is_a_required_tool() -> None:
    tools = (PROJECT_ROOT / "src/utils/tools.py").read_text()

    assert '"sickle", "bwa", "samtools", "java"' in tools
