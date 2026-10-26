rule core_trim_reads:
    input:
        r1=lambda wc: sample_r1(wc.sample),
        r2=lambda wc: sample_r2(wc.sample),
    output:
        r1_trim=temp(result_path("samples", "{sample}", "alignment_qc", "{sample}_1.fastq")),
        r2_trim=temp(result_path("samples", "{sample}", "alignment_qc", "{sample}_2.fastq")),
        s_trim=temp(result_path("samples", "{sample}", "alignment_qc", "{sample}_s.fastq")),
    threads:
        TRIM_READS_THREADS
    resources:
        trim_io=1
    params:
        min_length=TRIM_MIN_LENGTH_ARG,
    log:
        result_path("logs", "rules", "core_trim_reads", "{sample}.log"),
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.r1_trim})"
        sickle pe -t sanger {params.min_length} \
            -f "{input.r1}" \
            -r "{input.r2}" \
            -o "{output.r1_trim}" \
            -p "{output.r2_trim}" \
            -s "{output.s_trim}" \
            > {log} 2>&1
        '''


rule core_align_and_sort:
    input:
        r1_trim=rules.core_trim_reads.output.r1_trim,
        r2_trim=rules.core_trim_reads.output.r2_trim,
        s_trim=rules.core_trim_reads.output.s_trim,
    output:
        bam=result_path("samples", "{sample}", "alignment_qc", "{sample}.sort.bam"),
    threads:
        ALIGN_THREADS + SAMTOOLS_SORT_THREADS
    log:
        result_path("logs", "rules", "core_align_and_sort", "{sample}.log"),
    params:
        ref_fasta=REFERENCE_FASTA,
        ref_fai=REFERENCE_FASTA_FAI,
        bwa_index=BWA_INDEX_PREFIX,
        align_threads=ALIGN_THREADS,
        sort_threads=SAMTOOLS_SORT_THREADS,
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.bam})"

        if [ ! -f "{params.ref_fai}" ]; then
            samtools faidx "{params.ref_fasta}" >> {log} 2>&1
        fi

        bash src/scripts/NTM_bwamem.sh align "{params.ref_fasta}" "{params.bwa_index}" \
            "{input.r1_trim}" "{input.r2_trim}" "{input.s_trim}" "{output.bam}" \
            "{wildcards.sample}" {params.align_threads} {params.sort_threads} >> {log} 2>&1
        '''


rule core_depth_metrics:
    input:
        bam=rules.core_align_and_sort.output.bam,
    output:
        depth=result_path("samples", "{sample}", "report_inputs", "alignment_qc", "{sample}.depth"),
        status=result_path("samples", "{sample}", "alignment_qc", "core_gate_status.tsv"),
        mean_depth=result_path("samples", "{sample}", "report_inputs", "alignment_qc", "{sample}.mean_depth"),
    threads:
        DEPTH_METRICS_THREADS
    resources:
        depth_io=1
    log:
        result_path("logs", "rules", "core_depth_metrics", "{sample}.log"),
    params:
        genome_length=GENOME_LENGTH,
        min_mean_depth=0 if MAB_ROUTE else int(config.get("core_gate", {}).get("min_mean_depth", 10)),
        min_coverage_fraction=0 if MAB_ROUTE else float(config.get("core_gate", {}).get("min_coverage_fraction", 0.95)),
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.depth})" "$(dirname {output.status})"

        samtools depth "{input.bam}" > "{output.depth}" 2>> {log}

        # Single-pass awk computes both mean depth and coverage fraction.
        read -r mean_depth coverage_fraction <<< "$(awk -v g={params.genome_length} '{{sum+=$3; n+=1}} END {{if (n>0) printf "%.6f %.6f", sum/n, (g>0)?n/g:0; else printf "0 0"}}' "{output.depth}")"
        mean_round=$(printf '%.0f' "$mean_depth")

        passed=0
        if [ "$mean_round" -ge {params.min_mean_depth} ] && awk -v c="$coverage_fraction" 'BEGIN {{exit (c>={params.min_coverage_fraction})?0:1}}'; then
            passed=1
            msg="PASS"
        else
            msg="LOW_COVERAGE_OR_DEPTH"
        fi

        printf 'sample_id\tmean_depth\tcoverage_fraction\tpassed_core_gate\tmessage\n' > "{output.status}"
        printf '{wildcards.sample}\t%s\t%s\t%s\t%s\n' "$mean_depth" "$coverage_fraction" "$passed" "$msg" >> "{output.status}"
        printf '%.0f\n' "$mean_depth" > "{output.mean_depth}"
        '''


rule core_call_variants:
    input:
        bam=rules.core_align_and_sort.output.bam,
        depth=rules.core_depth_metrics.output.depth,
        status=rules.core_depth_metrics.output.status,
    output:
        done=result_path("logs", "rules", "core_call_variants", "{sample}.done"),
    threads:
        CALL_VARIANTS_THREADS
    resources:
        variant_io=1
    log:
        result_path("logs", "rules", "core_call_variants", "{sample}.log"),
    params:
        sample_dir=result_path("samples", "{sample}"),
        sample_id="{sample}",
        ref_fasta=REFERENCE_FASTA,
        varscan_jar=VARSCAN_JAR,
        ppe_list=PPE_LIST,
        mab_route=int(MAB_ROUTE),
        python_bin=PYTHON_BIN,
        cli_module=CLI_MODULE,
        discard_log=result_path("logs", "discard.txt"),
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.done})" "$(dirname {params.discard_log})" "{params.sample_dir}/variant_analysis"

        sample_id={params.sample_id}
        variant_dir="{params.sample_dir}/variant_analysis"
        varscan="$variant_dir/${{sample_id}}.varscan"
        cns="$variant_dir/${{sample_id}}.cns"
        vars="$variant_dir/${{sample_id}}.vars"
        var_ppe="$variant_dir/${{sample_id}}.var.ppe"
        var_for="$variant_dir/${{sample_id}}.var.for"
        snp="$variant_dir/${{sample_id}}.snp"

        passed=$(tail -n 1 "{input.status}" | cut -f4)
        mean_depth=$(tail -n 1 "{input.status}" | cut -f2)
        coverage_fraction=$(tail -n 1 "{input.status}" | cut -f3)

        if [ "{params.mab_route}" = "1" ]; then
            bash src/scripts/NTM_bwamem.sh call "{params.ref_fasta}" "{input.bam}" \
                "{params.varscan_jar}" "{params.ppe_list}" "$variant_dir" "$sample_id" \
                "{params.python_bin}" "{params.cli_module}" >> {log} 2>&1
        elif [ "$passed" = "1" ]; then
            min_cov=$(printf '%.0f' "$mean_depth")
            min_cov=$(( min_cov / 10 ))

            # Stream the mpileup into both VarScan consumers at once
            # (process substitution + tee) so the multi-GB pileup text
            # never touches disk. "-" tells VarScan to read stdin.
            samtools mpileup -q 30 -Q 30 -Bf "{params.ref_fasta}" "{input.bam}" \
                2>> {log} \
                | tee >(
                    java -jar "{params.varscan_jar}" mpileup2snp - \
                        --min-coverage "$min_cov" \
                        --min-reads2 2 --min-avg-qual 30 --min-var-freq 0.75 --p-value 99e-02 \
                        > "$varscan" 2>> {log}
                ) \
                | java -jar "{params.varscan_jar}" mpileup2cns - \
                    --min-coverage 3 --min-avg-qual 20 --min-var-freq 0.75 --strand-filter 0 --min-reads2 2 \
                    > "$cns" 2>> {log}

            # tee does not propagate the process-substitution exit code;
            # VarScan always writes a header, so an empty file means failure.
            if [ ! -s "$varscan" ]; then
                echo "VarScan mpileup2snp produced no output" >> {log}
                exit 1
            fi
            if [ ! -s "$cns" ]; then
                echo "VarScan mpileup2cns produced no output" >> {log}
                exit 1
            fi

            awk -F '[:]' '{{if($9==0 || $10==0)$0="";else print $0}}' "$varscan" > "$vars"
            {params.python_bin} -m {params.cli_module} ppe-filter --ppe-list "{params.ppe_list}" --input "$vars" --output "$var_ppe" >> {log} 2>&1
            {params.python_bin} -m {params.cli_module} format-trans --input "$var_ppe" --output "$var_for" >> {log} 2>&1
            cut -f2,3,4 "$var_for" > "$snp"

            rm -f "$varscan" "$var_ppe" "$var_for"
            echo "${{sample_id}} completed successfully" >> {log}
        else
            echo "${{sample_id}} do not meet criteria: ${{mean_depth}} ${{coverage_fraction}}" >> {params.discard_log}
            rm -f "$varscan" "$cns" "$vars" "$var_ppe" "$var_for" "$snp"
            echo "${{sample_id}} discarded (low coverage)" >> {log}
        fi

        touch {output.done}
        '''


rule core_collect_diff_loci:
    input:
        expand(result_path("logs", "rules", "core_call_variants", "{sample}.done"), sample=SAMPLES),
    output:
        diff_loci=result_path(".core_work", "diff_loci.txt"),
    threads:
        DIFF_LOCI_THREADS
    log:
        result_path("logs", "rules", "core_collect_diff_loci.log"),
    params:
        python_bin=PYTHON_BIN,
        cli_module=CLI_MODULE,
        samples_dir=result_path("samples"),
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.diff_loci})"
        {params.python_bin} -m {params.cli_module} diff-loci --snp-dir {params.samples_dir} --output {output.diff_loci} > {log} 2>&1
        '''


rule core_recall_consensus:
    input:
        diff_loci=rules.core_collect_diff_loci.output.diff_loci,
        call_done=rules.core_call_variants.output.done,
        mean_depth=rules.core_depth_metrics.output.mean_depth,
    output:
        done=result_path("logs", "rules", "core_recall_consensus", "{sample}.done"),
    threads:
        RECALL_THREADS
    log:
        result_path("logs", "rules", "core_recall_consensus", "{sample}.log"),
    params:
        python_bin=PYTHON_BIN,
        cli_module=CLI_MODULE,
        sample_id="{sample}",
        sample_dir=result_path("samples", "{sample}", "variant_analysis"),
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.done})" "{params.sample_dir}"
        cns="{params.sample_dir}/{params.sample_id}.cns"
        recall="{params.sample_dir}/{params.sample_id}.recall.fasta"

        if [ -s "$cns" ]; then
            {params.python_bin} -m {params.cli_module} recall \
                --loci {input.diff_loci} \
                --cns "$cns" \
                --depth {input.mean_depth} \
                --output "$recall" \
                > {log} 2>&1
        else
            rm -f "$recall"
            echo "skip recall: missing cns for {wildcards.sample}" > {log}
        fi

        touch {output.done}
        '''


rule core_merge_recall_fasta:
    input:
        expand(result_path("logs", "rules", "core_recall_consensus", "{sample}.done"), sample=SAMPLES),
    output:
        merged=result_path(".core_work", "merged.fasta"),
    threads:
        MERGE_THREADS
    log:
        result_path("logs", "rules", "core_merge_recall_fasta.log"),
    params:
        python_bin=PYTHON_BIN,
        cli_module=CLI_MODULE,
        samples_dir=result_path("samples"),
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.merged})"
        {params.python_bin} -m {params.cli_module} merge --fas-dir {params.samples_dir} --output {output.merged} > {log} 2>&1
        '''


rule core_extract_wildtype_loci:
    input:
        diff_loci=rules.core_collect_diff_loci.output.diff_loci,
    output:
        wildtype=result_path(".core_work", "wildtype.fasta"),
    threads:
        WILD_EXTRACT_THREADS
    log:
        result_path("logs", "rules", "core_extract_wildtype_loci.log"),
    params:
        python_bin=PYTHON_BIN,
        cli_module=CLI_MODULE,
        ancestor=ANCESTOR_FASTA,
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.wildtype})"
        {params.python_bin} -m {params.cli_module} wild-extract \
            --loci {input.diff_loci} \
            --ancestor {params.ancestor} \
            --output {output.wildtype} \
            > {log} 2>&1
        '''


rule core_filter_alignment:
    input:
        wild_loci=rules.core_extract_wildtype_loci.output.wildtype,
        merged=rules.core_merge_recall_fasta.output.merged,
    output:
        core_fa=result_path(".core_work", "core_snps.fadel-InvMisF5.bak.fa"),
        core_loc=result_path(".core_work", "core_snps.fadel-InvMisF5.bak.loc"),
    threads:
        FILTER_CORE_THREADS
    log:
        result_path("logs", "rules", "core_filter_alignment.log"),
    params:
        python_bin=PYTHON_BIN,
        cli_module=CLI_MODULE,
        output_prefix=result_path(".core_work", "core_snps"),
        reference_name=REFERENCE_LABEL,
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.core_fa})"
        {params.python_bin} -m {params.cli_module} filter \
            --wild-loci {input.wild_loci} \
            --alignment {input.merged} \
            --threshold 5 \
            --output-prefix {params.output_prefix} \
            --reference-name {params.reference_name} \
            > {log} 2>&1
        '''


rule core_compute_distance:
    input:
        core_fa=rules.core_filter_alignment.output.core_fa,
    output:
        distance=result_path(".core_work", "distance_matrix.txt"),
    threads:
        DISTANCE_THREADS
    log:
        result_path("logs", "rules", "core_compute_distance.log"),
    params:
        python_bin=PYTHON_BIN,
        cli_module=CLI_MODULE,
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.distance})"
        {params.python_bin} -m {params.cli_module} distance \
            --alignment "{input.core_fa}" \
            --output "{output.distance}" \
            --threads {threads} \
            > {log} 2>&1
        '''


rule core_publish_outputs:
    input:
        diff_loci=rules.core_collect_diff_loci.output.diff_loci,
        merged=rules.core_merge_recall_fasta.output.merged,
        wildtype=rules.core_extract_wildtype_loci.output.wildtype,
        core_fa=rules.core_filter_alignment.output.core_fa,
        core_loc=rules.core_filter_alignment.output.core_loc,
        distance=rules.core_compute_distance.output.distance,
    output:
        publish_done=result_path("logs", "workflow", "core_publish.done"),
        core_outputs=CORE_OUTPUTS,
    threads:
        1
    log:
        result_path("logs", "rules", "core_publish_outputs.log"),
    params:
        results_dir=RESULTS_DIR,
    shell:
        r'''
        set -euo pipefail
        core_dir="{params.results_dir}/core"
        mkdir -p "$(dirname {log})" "$(dirname {output.publish_done})" "$core_dir"

        cp {input.diff_loci} "$core_dir/diff_loci.txt"
        cp {input.merged} "$core_dir/merged.fasta"
        cp {input.wildtype} "$core_dir/wildtype.fasta"
        cp {input.core_fa} "$core_dir/core_snps.fadel-InvMisF5.bak.fa"
        cp {input.core_loc} "$core_dir/core_snps.fadel-InvMisF5.bak.loc"
        cp {input.distance} "$core_dir/distance_matrix.txt"

        echo "published shared core outputs to $core_dir" > {log}
        touch {output.publish_done}
        '''


rule run_core:
    input:
        publish_done=rules.core_publish_outputs.output.publish_done,
    output:
        done=result_path("logs", "workflow", "run_core.done"),
    threads:
        1
    log:
        result_path("logs", "rules", "run_core.log"),
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})" "$(dirname {output.done})"
        echo "core sample-level pipeline completed" > {log}
        touch {output.done}
        '''
