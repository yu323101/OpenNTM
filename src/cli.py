"""Command-line interface for MTB-Evo."""

import typer

from src.commands.build_reports import build_reports_cmd
from src.commands.cluster import cluster_cmd
from src.commands.compare_core import compare_core_cmd
from src.commands.diff_loci import diff_loci_cmd
from src.commands.distance import distance_cmd
from src.commands.filter import filter_cmd
from src.commands.filter_fasta import filter_fasta_cmd
from src.commands.format_trans import format_trans_cmd
from src.commands.investigation_alerts import investigation_alerts_cmd
from src.commands.drug_resistance import drug_resistance_cmd
from src.commands.lineage import lineage_cmd
from src.commands.merge import merge_cmd
from src.commands.ppe_filter import ppe_filter_cmd
from src.commands.prepare_downstream_inputs import prepare_downstream_inputs_cmd
from src.commands.recall import recall_cmd
from src.commands.report_figures import report_figures_cmd
from src.commands.report_tables import report_tables_cmd
from src.commands.wild_extract import wild_extract_cmd

app = typer.Typer(
    name="mtb-evo",
    help="Mycobacterium tuberculosis evolutionary analysis pipeline",
    no_args_is_help=True,
)

# Core commands
app.command("prepare-downstream-inputs")(prepare_downstream_inputs_cmd)
app.command("build-reports")(build_reports_cmd)
app.command("compare-core")(compare_core_cmd)
app.command("diff-loci")(diff_loci_cmd)
app.command("recall")(recall_cmd)
app.command("merge")(merge_cmd)
app.command("wild-extract")(wild_extract_cmd)
app.command("filter")(filter_cmd)
app.command("distance")(distance_cmd)
app.command("cluster")(cluster_cmd)
app.command("drug-resistance")(drug_resistance_cmd)
app.command("investigation-alerts")(investigation_alerts_cmd)
app.command("format-trans")(format_trans_cmd)
app.command("ppe-filter")(ppe_filter_cmd)
app.command("filter-fasta")(filter_fasta_cmd)

# Extension commands
app.command("lineage")(lineage_cmd)
app.command("report-figures")(report_figures_cmd)
app.command("report-tables")(report_tables_cmd)


def main() -> None:
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()
