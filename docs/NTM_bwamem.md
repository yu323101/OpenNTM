# NTM BWA-MEM route

Set `analysis_route: NTM_bwamem` in a configuration overlay. The default remains `bowtie2_legacy`.
Use a separate `results_dir` for each route to avoid reusing outputs from another route.

Configure an active reference profile with `fasta`, `fai`, `bwa_index`, `wildtype_fasta`,
`genome_length`, and **a reference-matched `ppe_list`**. No MAB mask or reference is
silently selected. The existing lineage/report assets are TB-specific: this route
supports NTM alignment and core SNP analysis; those assets must be replaced before
interpreting NTM reports.

Example overlay (replace paths with validated local files):

```yaml
analysis_route: NTM_bwamem
results_dir: results_ntm_bwamem
references:
  active: ntm
  profiles:
    ntm:
      fasta: /path/to/reference.fna
      fai: /path/to/reference.fna.fai
      bwa_index: /path/to/bwa_index/reference
      wildtype_fasta: /path/to/reference.fna
      genome_length: 5402402
      ppe_list: /path/to/reference_matched_mask.txt
```

Run `snakemake --configfile config/config.yaml config/ntm.yaml --cores 16 core_only`.

## Parameters and outputs

Sickle `pe -l 35 -t sanger`; BWA-MEM `-c 100 -M` with sample read groups,
paired and retained singleton reads aligned separately. `threads.bwa` controls
threads instead of the source script's fixed one thread. SAM/BAM conversion follows
the valid source scripts, correcting the bolletii script's reversed input/output.
Index building is serialized with `flock` (util-linux required).

mpileup uses `-q 30 -Q 20 -B -O`. A temporary pileup is retained through calling
and `.forup` generation, then removed. VarScan 2.3.9 SNP calling uses minimum
coverage 3, reads2 2, average quality 20, variant frequency 0.01, hom frequency 0.9,
p-value 0.99, strand-filter 0. CNS uses coverage 3, reads2 2, average quality 20,
variant frequency 0.75, strand-filter 0. No additional two-strand filter is applied.
The source MAB route has no 10x/95% sample gate; this route records coverage but
bypasses that legacy gate.

`.vars` contains raw candidate calls (including low-frequency calls); `.cns` contains
consensus calls. Reference-mask filtering and existing format conversion precede
`.snp`, which retains calls at frequency >=75%, matching `2_fix_extract.pl`.
`.forup` includes masked candidate records, read bases, read-position lists from
`-O`, and base qualities, matching `3.1_mix_pileup_merge.pl`.
The two small Perl helpers are copied unchanged from the user-specified workflow.

Low-frequency candidate calls are not automatically validated environmental detections.
This implementation adds the analysis route; it does not establish detection sensitivity
or specificity for environmental samples.
