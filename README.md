# mtb-evo

`mtb-evo` 是一个基于 Snakemake 的结核分枝杆菌分析流程，输入样本 `FASTQ`，输出样本级图表与表格，以及全样本共享 `core` 结果。

## 1. 当前流程（唯一主流程）

主入口：`Snakemake`

```bash
snakemake --snakefile Snakefile --configfile config/config.yaml -j 4 all
```

可选 target：

- `all`：全流程（Layer A + B + C + D）
- `foundation_only`：基础输入层（Layer A）
- `core_only`：核心分析层（Layer B）
- `downstream_only`：仅桥接输入（Layer C）
- `reports_only`：仅报告层（Layer D）

## 2. 输入格式

推荐使用 `config/samplesheet.csv`：

```csv
sample_id,r1,r2
MD001.cleaned,/path/MD001.cleaned_1.fastq.gz,/path/MD001.cleaned_2.fastq.gz
MD002.cleaned,/path/MD002.cleaned_1.fastq.gz,/path/MD002.cleaned_2.fastq.gz
```

运行前请先改写 `config/samplesheet.csv`，或在命令行覆盖 `samplesheet`：

```bash
snakemake --snakefile Snakefile --configfile config/config.yaml -n all --config samplesheet=/abs/path/to/samplesheet.csv
```

## 3. 分层设计（简版）

- Layer A：FASTQ 派生基础输入（`fastp_qc` JSON）
- Layer B：样本主分析（BAM / CNS / SNP / vars）+ 共享 `core`
- Layer C：下游桥接输入（`depth` / `insert_sizes` / `*_annotated.txt`）
- Layer D：最终报告（`figure/`、`table/`、`lineage_summary`）

## 4. 结果目录结构

```text
results/
├── core/
│   ├── diff_loci.txt
│   ├── merged.fasta
│   ├── wildtype.fasta
│   ├── core_snps.fadel-InvMisF5.bak.fa
│   ├── core_snps.fadel-InvMisF5.bak.loc
│   └── distance_matrix.txt
├── samples/
│   └── <sample>/
│       ├── alignment_qc/
│       ├── variant_analysis/
│       ├── report_inputs/
│       │   ├── fastp_qc/
│       │   ├── alignment_qc/
│       │   └── annotated_variants/
│       ├── figure/
│       └── table/
└── logs/
    ├── rules/
    ├── workflow/
    ├── runs/<run_id>/
    ├── prepare_foundation_outputs_status.tsv
    ├── prepare_downstream_inputs_status.tsv
    ├── lineage_summary.tsv
    ├── report_figures_status.tsv
    └── report_tables_status.tsv
```

说明：

- `figure/` 是正式目录名（单数）。
- `results/core/` 为全样本共享结果，样本目录不再存重复 `core` 副本。
- `.core_work/` 仅是流程内部中转目录，不作为正式交付结果。

## 5. 状态文件与 run_id

状态文件按 run_id 分批输出到：

- `results/logs/runs/<run_id>/*.tsv`

同时覆盖一份固定“latest 视图”：

- `results/logs/*.tsv`

默认 `run_id=latest`。可手动指定批次：

```bash
snakemake --snakefile Snakefile --configfile config/config.yaml -j 4 all --config run_id=run_20260418_a
```

状态 TSV 核心字段：

- `sample_id`
- `status`
- `issues`
- `run_id`
- `stage`

## 6. 参考序列配置

在 `config/config.yaml` 使用 profile 方式切换：

```yaml
references:
  active: tb_ancestor   # 或 tb_h37rv
```

当前默认：`tb_ancestor`（`data/tb.ancestor.fasta`）。

## 7. 图表与表格口径

- Figure 输出 `01-09`，同时生成 `PNG + PDF`。
- Figure/Table 使用同一套变异分类逻辑（`.cns` 主表 + 注释展开回填），避免口径分叉。
- 谱系结果保持独立表输出：`谱系鉴定结果_<sample>.csv`，不并入 `临床变异检测报告_<sample>.csv`。

## 8. 常用命令

全流程：

```bash
snakemake --snakefile Snakefile --configfile config/config.yaml -j 4 all
```

只跑基础输入层（A）：

```bash
snakemake --snakefile Snakefile --configfile config/config.yaml -j 4 foundation_only
```

只跑核心分析层（B）：

```bash
snakemake --snakefile Snakefile --configfile config/config.yaml -j 4 core_only
```

只补桥接输入（C）：

```bash
snakemake --snakefile Snakefile --configfile config/config.yaml -j 4 downstream_only
```

只跑报告层（D）：

```bash
snakemake --snakefile Snakefile --configfile config/config.yaml -j 4 reports_only
```

dry-run：

```bash
snakemake --snakefile Snakefile --configfile config/config.yaml -n all
```

或显式覆盖样本表：

```bash
snakemake --snakefile Snakefile --configfile config/config.yaml -n all --config samplesheet=/abs/path/to/samplesheet.csv
```

## 9. 当前边界与后续

当前已稳定：

- Snakemake 单轨编排
- `results/core/` 单份共享产物
- 样本级 `figure/`、`table/` 输出
- `run_id` 分批状态追踪

后续建议优先项：

1. 补充回归测试（图表/表格口径一致性）。
2. 继续细化 core 规则并增强失败诊断日志。

### NTM BWA-MEM analysis route

The optional `analysis_route: NTM_bwamem` follows the user-provided MAB workflow,
including mpileup `-O`, low-frequency candidate calls and >=75% final SNP calls.
See [configuration and output details](docs/NTM_bwamem.md).
