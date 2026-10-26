from pathlib import Path

from src.commands.filter import filter_core_snp

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_reference_row_uses_active_reference_name(tmp_path: Path) -> None:
    wild = tmp_path / "wild.fasta"
    wild.write_text("100\tA\n101\tC\n")
    alignment = tmp_path / "merged.fasta"
    alignment.write_text(">s1\nAC\n>s2\nAG\n")

    filter_core_snp(wild, alignment, 5, str(tmp_path / "core"), reference_name="tb.ancestor")

    output = tmp_path / "core.fadel-InvMisF5.bak.fa"
    headers = [line.strip() for line in output.read_text().splitlines() if line.startswith(">")]
    assert headers == [">s1", ">s2", ">tb.ancestor"]
    assert ">H37Rv" not in output.read_text()


def test_filter_rule_passes_reference_label() -> None:
    core_rule = (PROJECT_ROOT / "rules/core.smk").read_text()
    common_rule = (PROJECT_ROOT / "rules/common.smk").read_text()

    assert 'REFERENCE_LABEL = str(ACTIVE_REFERENCE_CFG.get("label", REFERENCE_PROFILE))' in common_rule
    assert "--reference-name {params.reference_name}" in core_rule
