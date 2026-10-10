#!/usr/bin/env bash
# MAB route: paired and retained singleton reads, followed by legacy MAB calling.
set -euo pipefail
stage=$1; shift
if [[ "$stage" == align ]]; then
    ref=$1; prefix=$2; r1=$3; r2=$4; singles=$5; output=$6; sample=$7; threads=$8; sort_threads=$9
    mkdir -p "$(dirname "$prefix")" "$(dirname "$output")"
    (
        flock 9
        if [[ ! -f "$prefix.bwt" || ! -f "$prefix.sa" ]]; then bwa index -p "$prefix" "$ref"; fi
    ) 9>"$prefix.index.lock"
    tmp=$(mktemp -d "$(dirname "$output")/.bwa.XXXXXX")
    trap 'rm -rf "$tmp"' EXIT
    rg="@RG\\tID:${sample}\\tSM:${sample}\\tPL:illumina"
    bwa mem -t "$threads" -c 100 -R "$rg" -M "$prefix" "$r1" "$r2" | samtools view -b -o "$tmp/paired.bam" -
    bwa mem -t "$threads" -c 100 -R "$rg" -M "$prefix" "$singles" | samtools view -b -o "$tmp/single.bam" -
    samtools merge -u "$tmp/merged.bam" "$tmp/paired.bam" "$tmp/single.bam"
    samtools sort -@ "$sort_threads" -o "$output" "$tmp/merged.bam"
    samtools index "$output"
    samtools flagstat "$output"
else
    [[ "$stage" == call ]] || exit 2
    ref=$1; bam=$2; jar=$3; mask=$4; directory=$5; sample=$6; python_bin=$7; cli=$8
    mkdir -p "$directory"
    tmp=$(mktemp -d "$directory/.mab-call.XXXXXX")
    trap 'rm -rf "$tmp"' EXIT
    samtools mpileup -q 30 -Q 20 -B -O -f "$ref" "$bam" > "$tmp/pileup"
    java -jar "$jar" mpileup2snp "$tmp/pileup" --min-coverage 3 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.01 --min-freq-for-hom 0.9 --p-value 99e-02 --strand-filter 0 > "$directory/$sample.vars"
    java -jar "$jar" mpileup2cns "$tmp/pileup" --min-coverage 3 --min-avg-qual 20 --min-var-freq 0.75 --min-reads2 2 --strand-filter 0 > "$directory/$sample.cns"
    "$python_bin" -m "$cli" ppe-filter --ppe-list "$mask" --input "$directory/$sample.vars" --output "$tmp/ppe"
    "$python_bin" -m "$cli" format-trans --input "$tmp/ppe" --output "$tmp/formatted"
    perl src/scripts/mab_fix_extract.pl "$tmp/formatted" | cut -f2-4 > "$directory/$sample.snp"
    perl src/scripts/mab_mix_pileup_merge.pl "$tmp/formatted" "$tmp/pileup" > "$directory/$sample.forup"
fi
