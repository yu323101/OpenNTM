configfile: "config/config.yaml"

include: "rules/common.smk"
include: "rules/foundation.smk"
include: "rules/alignment_variant.smk"
include: "rules/core.smk"
include: "rules/downstream_inputs.smk"
include: "rules/lineage.smk"
include: "rules/extensions.smk"
include: "rules/reports.smk"

localrules: all, foundation_only, core_only, downstream_only, reports_only

rule all:
    input:
        rules.build_reports.output.done,

rule foundation_only:
    input:
        rules.foundation_prepare_outputs.output.done,

rule core_only:
    input:
        rules.run_core.output.done,

rule downstream_only:
    input:
        rules.prepare_downstream_inputs.output.done,

rule reports_only:
    input:
        rules.build_reports.output.done,
