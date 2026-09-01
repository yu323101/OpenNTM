from pathlib import Path
from types import SimpleNamespace

import pytest

import src.commands.distance as distance_module
from src.commands.distance import (
    build_snp_dists_command,
    calculate_distance,
    molten_to_upper_triangle,
)


def test_snp_dists_command_preserves_legacy_base_rules_and_default_thread() -> None:
    command = build_snp_dists_command(Path("alignment.fa"))

    assert command == [
        "snp-dists",
        "-q",
        "-k",
        "-L",
        "-m",
        "-t",
        "-j",
        "1",
        "alignment.fa",
    ]
    assert "-a" not in command


def test_snp_dists_command_accepts_configured_threads() -> None:
    command = build_snp_dists_command(Path("alignment.fa"), threads=8)

    assert command[command.index("-j") + 1] == "8"


def test_molten_output_is_converted_to_legacy_upper_triangle(tmp_path: Path) -> None:
    molten = (
        "sequence_1\tsequence_2\tdistance\n"
        "S1\tS1\t0\n"
        "S1\tS2\t2\n"
        "S1\tS3\t5\n"
        "S2\tS1\t2\n"
        "S2\tS2\t0\n"
        "S2\tS3\t3\n"
        "S3\tS1\t5\n"
        "S3\tS2\t3\n"
        "S3\tS3\t0\n"
    )
    output = tmp_path / "distance_matrix.txt"

    molten_to_upper_triangle(molten, output, ["S1", "S2", "S3"])

    assert output.read_text(encoding="utf-8") == "S1\tS2\t2\nS1\tS3\t5\nS2\tS3\t3\n"


def test_molten_output_accepts_reverse_only_pairs(tmp_path: Path) -> None:
    output = tmp_path / "distance.tsv"
    molten_to_upper_triangle(
        "sequence_1 sequence_2 distance\nS2 S1 7\n",
        output,
        ["S1", "S2"],
    )

    assert output.read_text(encoding="utf-8") == "S1\tS2\t7\n"


def test_molten_output_rejects_asymmetric_matrix(tmp_path: Path) -> None:
    output = tmp_path / "distance.tsv"
    molten = "S1\tS2\t2\nS2\tS1\t3\n"

    with pytest.raises(ValueError, match="Asymmetric distances"):
        molten_to_upper_triangle(molten, output, ["S1", "S2"])


def test_calculate_distance_invokes_snp_dists_and_writes_legacy_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    alignment = tmp_path / "alignment.fa"
    alignment.write_text(">S1\nACGT\n>S2\nAGGT\n", encoding="utf-8")
    output = tmp_path / "nested" / "distance.tsv"
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(stdout="sequence_1\tsequence_2\tdistance\nS1\tS2\t1\n")

    monkeypatch.setattr(distance_module.subprocess, "run", fake_run)

    calculate_distance(alignment, output, threads=4)

    assert calls[0][0] == build_snp_dists_command(alignment, threads=4)
    assert calls[0][1]["check"] is True
    assert output.read_text(encoding="utf-8") == "S1\tS2\t1\n"


def test_distance_threads_are_configured_to_one_by_default() -> None:
    config_text = Path("config/config.yaml").read_text(encoding="utf-8")

    assert "  distance: 1" in config_text
    assert "--threads {threads}" in Path("rules/core.smk").read_text(encoding="utf-8")
