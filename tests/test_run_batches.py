import importlib.util
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "run_batches", PROJECT_ROOT / "scripts/run_batches.py"
)
run_batches = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(run_batches)


def test_resolve_substitutes_sample(tmp_path: Path) -> None:
    target = run_batches.resolve(tmp_path, "logs/rules/core_call_variants/{sample}.done", "S1")
    assert target == tmp_path / "logs/rules/core_call_variants/S1.done"


def test_order_samples_skips_completed_targets(tmp_path: Path) -> None:
    template = "done/{sample}.done"
    (tmp_path / "done").mkdir()
    (tmp_path / "done/S1.done").write_text("")

    order = run_batches.order_samples(["S1", "S2", "S3"], tmp_path, template)

    assert order == ["S2", "S3"]


def test_order_samples_prefers_samples_with_existing_input(tmp_path: Path) -> None:
    (tmp_path / "trim").mkdir()
    (tmp_path / "trim/S3_1.fastq").write_text("")

    order = run_batches.order_samples(
        ["S1", "S2", "S3", "S4"],
        tmp_path,
        "done/{sample}.done",
        "trim/{sample}_1.fastq",
    )

    assert order == ["S3", "S1", "S2", "S4"]
