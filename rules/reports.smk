rule build_reports:
    input:
        samplesheet=str(SAMPLESHEET),
        downstream_done=rules.prepare_downstream_inputs.output.done,
        epidemiology_done=([] if SKIP_EXTENSIONS else rules.generate_investigation_alerts.output.done),
    output:
        done=result_path("logs", "workflow", "build_reports.done"),
        lineage_summary=LINEAGE_SUMMARY_OUTPUTS,
        figure_status=FIGURE_STATUS_OUTPUTS,
        table_status=TABLE_STATUS_OUTPUTS,
    log:
        result_path("logs", "rules", "build_reports.log"),
    params:
        results_dir=RESULTS_DIR,
        python_bin=PYTHON_BIN,
        cli_module=CLI_MODULE,
        verbose_flag=VERBOSE_FLAG,
        strict_flag=STRICT_FLAG,
        skip_lineage_flag=SKIP_LINEAGE_FLAG,
        skip_figures_flag=SKIP_FIGURES_FLAG,
        skip_tables_flag=SKIP_TABLES_FLAG,
        # Extension rules, when enabled, already produce these files before
        # build-reports; avoid running the same steps a second time here.
        skip_extensions_flag="--skip-extensions",
        run_id=RUN_ID,
        lineage_run=LINEAGE_SUMMARY_RUN,
        lineage_latest=LINEAGE_SUMMARY_LATEST,
        figure_run=FIGURE_STATUS_RUN,
        figure_latest=FIGURE_STATUS_LATEST,
        table_run=TABLE_STATUS_RUN,
        table_latest=TABLE_STATUS_LATEST,
    shell:
        r'''
        mkdir -p "$(dirname {log})" "$(dirname {output.done})"
        {params.python_bin} -m {params.cli_module} build-reports \
            --samplesheet {input.samplesheet} \
            --output-dir {params.results_dir} \
            --lineage-summary-output {params.lineage_run} \
            --figure-status-output {params.figure_run} \
            --table-status-output {params.table_run} \
            --run-id {params.run_id} \
            {params.verbose_flag} \
            {params.strict_flag} \
            {params.skip_lineage_flag} \
            {params.skip_figures_flag} \
            {params.skip_tables_flag} \
            {params.skip_extensions_flag} \
            &> {log}
        if [ -f "{params.lineage_run}" ]; then cp {params.lineage_run} {params.lineage_latest}; fi
        if [ -f "{params.figure_run}" ]; then cp {params.figure_run} {params.figure_latest}; fi
        if [ -f "{params.table_run}" ]; then cp {params.table_run} {params.table_latest}; fi
        touch {output.done}
        '''
