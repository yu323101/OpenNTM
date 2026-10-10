"""Epidemiology-supporting extensions: proximity clusters, resistance screen and alerts."""

rule cluster_samples:
    input:
        core_done=rules.run_core.output.done,
        distance=result_path("core", "distance_matrix.txt"),
    output:
        assignments=result_path("core", "cluster_assignments.tsv"),
        done=result_path("logs", "workflow", "cluster_samples.done"),
    threads:
        CLUSTER_THREADS
    log:
        result_path("logs", "rules", "cluster_samples.log"),
    params:
        threshold=CLUSTER_THRESHOLD,
        python_bin=PYTHON_BIN,
        cli_module=CLI_MODULE,
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.assignments})" "$(dirname {output.done})"
        {params.python_bin} -m {params.cli_module} cluster \
            --distance-matrix {input.distance} \
            --output {output.assignments} \
            --threshold {params.threshold} > {log} 2>&1
        touch {output.done}
        '''


rule classify_drug_resistance:
    input:
        downstream_done=rules.prepare_downstream_inputs.output.done,
        drug_db=DRUG_DB,
    output:
        manifest=result_path("logs", "drug_resistance_manifest.tsv"),
        done=result_path("logs", "workflow", "classify_drug_resistance.done"),
    threads:
        DRUG_RESISTANCE_THREADS
    log:
        result_path("logs", "rules", "classify_drug_resistance.log"),
    params:
        samples_dir=result_path("samples"),
        output_dir=RESULTS_DIR,
        python_bin=PYTHON_BIN,
        cli_module=CLI_MODULE,
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.manifest})" "$(dirname {output.done})"
        {params.python_bin} -m {params.cli_module} drug-resistance \
            --samples-dir {params.samples_dir} \
            --output-dir {params.output_dir} \
            --drug-db {input.drug_db} \
            --manifest-output {output.manifest} > {log} 2>&1
        touch {output.done}
        '''


rule generate_investigation_alerts:
    input:
        cluster_done=rules.cluster_samples.output.done,
        cluster_assignments=rules.cluster_samples.output.assignments,
        resistance_done=rules.classify_drug_resistance.output.done,
        resistance_manifest=rules.classify_drug_resistance.output.manifest,
    output:
        alerts=result_path("logs", "investigation_alerts.tsv"),
        done=result_path("logs", "workflow", "generate_investigation_alerts.done"),
    threads:
        INVESTIGATION_ALERTS_THREADS
    log:
        result_path("logs", "rules", "generate_investigation_alerts.log"),
    params:
        drug_resistance_dir=result_path("samples"),
        python_bin=PYTHON_BIN,
        cli_module=CLI_MODULE,
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.alerts})" "$(dirname {output.done})"
        {params.python_bin} -m {params.cli_module} investigation-alerts \
            --cluster-assignments {input.cluster_assignments} \
            --drug-resistance-dir {params.drug_resistance_dir} \
            --output {output.alerts} > {log} 2>&1
        touch {output.done}
        '''
