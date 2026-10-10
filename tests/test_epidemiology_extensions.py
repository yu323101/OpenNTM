from pathlib import Path

from src.commands.cluster import assign_clusters
from src.commands.drug_resistance import classify_resistance, screen_resistance, screen_samples
from src.commands.investigation_alerts import generate_alerts


def test_assign_clusters_is_stable_and_uses_transitive_components(tmp_path: Path) -> None:
    distances = tmp_path / "distance.tsv"
    distances.write_text("a\tb\t2\na\tc\t8\nb\tc\t2\n", encoding="utf-8")
    output = tmp_path / "clusters.tsv"

    assign_clusters(distances, output, threshold=5)

    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "sample_id\tcluster_id\tcluster_size\tsnp_threshold"
    assert lines[1:] == ["a\tC0001\t3\t5", "b\tC0001\t3\t5", "c\tC0001\t3\t5"]


def test_resistance_screen_uses_current_who_pre_xdr_and_xdr_rules(tmp_path: Path) -> None:
    database = tmp_path / "drug.tsv"
    database.write_text(
        "Drug\tLocus\tChromosome Coordinate\tGene Coordinate\tNucleotide Change\tCodon\tCodon Change\tAA\tTag\tSources\n"
        "ISONIAZID\tkatG\t10\t-\tA/G\t-\t-\t-\t-\ttest\n"
        "RIFAMPICIN\trpoB\t20\t-\tC/T\t-\t-\t-\t-\ttest\n"
        "FLUOROQUINOLONES\tgyrA\t30\t-\tG/A\t-\t-\t-\t-\ttest\n"
        "BEDAQUILINE\tRv0678\t40\t-\tT/C\t-\t-\t-\t-\ttest\n",
        encoding="utf-8",
    )
    annotated = tmp_path / "S1_annotated.txt"
    annotated.write_text("10\tA\tG\t-\t-\t-\t-\t-\t-\t-\n"
                         "20\tC\tT\t-\t-\t-\t-\t-\t-\t-\n"
                         "30\tG\tA\t-\t-\t-\t-\t-\t-\t-\n"
                         "40\tT\tC\t-\t-\t-\t-\t-\t-\t-\n", encoding="utf-8")
    output = tmp_path / "resistance.tsv"

    screen_resistance(annotated, database, output, sample_id="S1")

    row = output.read_text(encoding="utf-8").splitlines()[1].split("\t")
    assert row[1] == "XDR"
    assert classify_resistance(["ISONIAZID", "RIFAMPICIN", "FLUOROQUINOLONES"]) == "pre-XDR"


def test_alerts_are_actionable_but_do_not_claim_transmission(tmp_path: Path) -> None:
    clusters = tmp_path / "clusters.tsv"
    clusters.write_text("sample_id\tcluster_id\tcluster_size\tsnp_threshold\nS1\tC0001\t2\t5\n", encoding="utf-8")
    resistance_dir = tmp_path / "samples" / "S1" / "table"
    resistance_dir.mkdir(parents=True)
    (resistance_dir / "S1_drug_resistance.tsv").write_text(
        "sample_id\tclassification\tresistant_drugs\nS1\tMDR\tISONIAZID;RIFAMPICIN\n",
        encoding="utf-8",
    )
    output = tmp_path / "alerts.tsv"

    generate_alerts(clusters, tmp_path / "samples", output)

    text = output.read_text(encoding="utf-8")
    assert "GENOMIC_PROXIMITY_CLUSTER" in text
    assert "DRUG_RESISTANCE_RISK" in text
    assert "不代表已证实传播" in text


def test_batch_resistance_screen_supports_legacy_inputs_and_standard_output_layout(
    tmp_path: Path,
) -> None:
    database = tmp_path / "drug.tsv"
    database.write_text(
        "Drug\tLocus\tChromosome Coordinate\tGene Coordinate\tNucleotide Change\n"
        "RIFAMPICIN\trpoB\t20\t-\tC/T\n",
        encoding="utf-8",
    )
    annotated_dir = tmp_path / "samples" / "S1" / "report_inputs" / "variant_analysis"
    annotated_dir.mkdir(parents=True)
    (annotated_dir / "S1_annotated.txt").write_text("20\tC\tT\n", encoding="utf-8")

    manifest = tmp_path / "results" / "logs" / "drug_resistance_manifest.tsv"
    screen_samples(tmp_path / "samples", database, tmp_path / "results", manifest)

    sample_result = tmp_path / "results" / "samples" / "S1" / "table" / "S1_drug_resistance.tsv"
    assert sample_result.exists()
    assert "MTB_TB_ONLY" in sample_result.read_text(encoding="utf-8")
    assert "S1\tRR\tRIFAMPICIN" in manifest.read_text(encoding="utf-8")
