# mtb-evo

`mtb-evo` 是一个基于 Snakemake 的分枝杆菌 WGS/NGS 分析流程，默认配置面向
MTB/TB；输入样本 `FASTQ`，输出样本级图表与表格，以及全样本共享 `core` 结果。
研究 NTM 时必须使用匹配的参考序列、注释资源和阈值；仓库内置的 TB 专用资源
不能直接解释 NTM。

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
│           └── <sample>_drug_resistance.tsv
└── logs/
    ├── rules/
    ├── workflow/
    ├── runs/<run_id>/
    ├── prepare_foundation_outputs_status.tsv
    ├── prepare_downstream_inputs_status.tsv
    ├── lineage_summary.tsv
    ├── report_figures_status.tsv
    ├── report_tables_status.tsv
    ├── drug_resistance_manifest.tsv
    └── investigation_alerts.tsv
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

## 10. 流行病防控支持扩展

`all` 在核心距离和样本注释准备完成后，会额外生成：

- `results/core/cluster_assignments.tsv`：按 `clustering.snp_threshold` 将样本构建为基因组近缘性连通簇。
- `results/logs/drug_resistance_manifest.tsv`：将已知耐药突变标记映射到样本，并按 WHO MDR/pre-XDR/XDR 规则分级。
- `results/logs/investigation_alerts.tsv`：生成近缘簇和耐药风险的可解释流调跟进提示。

这些模块用于疾病预防控制中的筛查、风险分层和人工流调分流。近缘簇是距离图的
连通分量，不是已证实的传播链；必须结合采样时间、地点、接触史和其他流调资料。
耐药模块使用仓库内置的 MTB/TB 专用数据库，只输出基因型标记筛查结果；它不适用
于 NTM，也不能替代表型药敏（DST）、临床资料或治疗决策。NTM 样本可以在参考序列
和变异定义合适时使用通用聚类模块，但仍需使用 NTM 适用的数据库/阈值并结合流调。

也可以单独使用：

```bash
python -m src.cli cluster --distance-matrix results/core/distance_matrix.txt \
  --output results/core/cluster_assignments.tsv --threshold 12
python -m src.cli drug-resistance --samples-dir results/samples \
  --output-dir results --drug-db data/drug_db/Total_drug_resistance_mutation.txt
python -m src.cli investigation-alerts \
  --cluster-assignments results/core/cluster_assignments.tsv \
  --drug-resistance-dir results/samples --output results/logs/investigation_alerts.tsv
```

这些输出是防控筛查和人工复核的辅助证据：SNP 近缘性不等于已证实传播；未检出数据库标记也不等于药敏；耐药分级不能替代表型药敏、临床资料和专业人员判断。

Snakemake 默认启用这些扩展；在 `config/config.yaml` 中设置
`report.skip_extensions: true` 可以把扩展规则从 DAG 中排除。直接运行
`build-reports` 时使用 `--skip-extensions` 达到同样效果。扩展规则本身使用单线程
调度槽（`threads.cluster`、`threads.drug_resistance`、
`threads.investigation_alerts`）；它们当前是轻量串行步骤，主要用于资源调度。

## 11. 两两 SNP 距离实现更新

核心 alignment 的两两距离现在由 `snp-dists` 计算，以替代原来的 Python
逐字符双重循环。`snp-dists` 使用 C/OpenMP 实现，适合常规二代测序分析中对较多
样本进行距离计算；距离口径保持不变：只比较 `A/C/G/T`，忽略 `N`、`?`、gap
以及其他非标准字符，不使用 `-a`，也不做进化模型校正。

默认使用单线程，避免在共享计算节点上意外占用过多资源。线程数在
`config/config.yaml` 的 `threads.distance` 中修改，例如：

```yaml
threads:
  # snp-dists threads; default 1
  distance: 4
```

该参数是 `snp-dists` 的线程数；Snakemake 命令行的 `-j` 仍然控制可并行运行的
任务数。流程内部调用 `snp-dists -q -k -L -m -t -j N`，只计算唯一的样本对，
并将其结果转换回旧的
上三角三列格式：

```text
sample_1<TAB>sample_2<TAB>distance
```

因此现有 `results/core/distance_matrix.txt` 和读取该文件的下游接口保持兼容。
前置 alignment 构建、PPE/高变异位点屏蔽和 missing 过滤逻辑没有改变。

这里的距离结果用于常规 NGS 数据分析和潜在基因组近缘簇的筛查/预防支持，不是
临床诊断或传播认定；本次更新不引入 TB 专用距离阈值，也不产生临床结论。任何
聚类阈值都应结合研究对象、采样背景和流行病学资料单独解释。
