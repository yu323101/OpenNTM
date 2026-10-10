# OpenNTM

OpenNTM 是基于 Snakemake 的分枝杆菌全基因组测序分析项目，支持 **TB 和 NTM 两条技术路线**。输入双端 FASTQ，输出比对 BAM、共识 CNS、候选变异、筛选后 SNP、核心 SNP alignment 及样本间 SNP 距离。NTM 路线面向非结核分枝杆菌，包括自然环境样本的检测研究、基因组比较及预防分析。

## 两条技术路线

### TB：`TB_bowtie2`

沿用原 `pair_fixed_nostrandbias.py` 的主要分析逻辑，是当前默认路线。

```text
FASTQ → Sickle → Bowtie2 默认 end-to-end 比对 → samtools 排序
      → 深度/覆盖率门槛 → mpileup → VarScan
      → 双链支持筛选 → 参考屏蔽 → CNS / SNP → 核心 SNP / 距离
```

### NTM：`NTM_bwamem`

参考原 MAB BWA-MEM 脚本。双端 reads 和修剪后保留的单端 reads 分别比对，随后合并。

```text
FASTQ → Sickle（最短 35 bp）→ BWA-MEM → BAM 合并、排序、索引
      → mpileup（含 -O）→ VarScan 候选变异 / CNS
      → 参考屏蔽 → ≥75% 最终 SNP + .forup → 核心 SNP / 距离
```

| 设置 | TB 路线 | NTM 路线 |
|---|---|---|
| `analysis_route` | `TB_bowtie2` | `NTM_bwamem` |
| 比对 | Bowtie2 默认 sensitive、end-to-end | BWA-MEM `-c 100 -M`，带样本 read group |
| Sickle | 默认长度参数 | `-l 35` |
| 样本准入 | 取整平均深度 ≥10×、覆盖 ≥95%，可配置 | 记录覆盖情况，不应用 TB 准入门槛 |
| mpileup | `-q 30 -Q 30 -B` | `-q 30 -Q 20 -B -O` |
| SNP 候选最低深度 | 取整平均深度整除 10 | 固定 3× |
| SNP 候选最低平均质量 | 30 | 20 |
| SNP 候选最低频率 | 75% | 1% |
| SNP 候选链过滤 | VarScan 默认过滤，并要求双链支持 | `--strand-filter 0`，无额外双链筛选 |
| 最终 `.snp` | 按上述规则及屏蔽文件筛选 | 屏蔽后保留频率 ≥75% 的记录 |
| 共识 CNS | 深度 ≥3、平均质量 ≥20、频率 ≥75%，关闭链过滤 | 相同 |
| `.forup` | 不生成 | 保留候选变异及 read 内位置、碱基和质量 |

两条路线的 SNP 候选均要求至少 2 条变异支持 reads，p-value 参数为 0.99。NTM 的 `--min-freq-for-hom` 为 0.9。`-O` 增加 read 内碱基位置列，不改变比对或质量阈值。低频候选与最终 SNP 分开输出，便于环境样本研究中的复核。

## 安装与入口

依赖包括 Sickle、Bowtie2、BWA、samtools、VarScan、Perl、snp-dists、Snakemake 和 Python。BWA 索引构建使用 `flock`，需要 util-linux。使用仓库中的 Conda 环境定义：

```bash
conda env create -f environment.yml
conda activate mtb-evo
pip install -e .
openntm --help
```

安装后 CLI 名称为 `openntm`；旧的 `mtb-evo` 命令保留兼容。现有 Conda 环境名称仍为 `mtb-evo`。主要工作流入口为仓库根目录的 `Snakefile`。

## 输入与参考配置

样本表采用 CSV 格式：

```csv
sample_id,r1,r2
sample01,/abs/path/sample01_1.fastq.gz,/abs/path/sample01_2.fastq.gz
sample02,/abs/path/sample02_1.fastq.gz,/abs/path/sample02_2.fastq.gz
```

TB 路线默认使用 `tb_ancestor`，也可以选择 `tb_h37rv`。NTM 必须配置研究对象对应的参考序列、参考长度和屏蔽坐标，不能直接复用 TB 或其他参考的坐标文件。

NTM 配置覆盖文件示例，保存为 `config/ntm.yaml`，将示例路径及长度替换为实际资源：

```yaml
analysis_route: NTM_bwamem
samplesheet: /abs/path/ntm_samples.csv
results_dir: results_ntm_bwamem
references:
  active: ntm
  profiles:
    ntm:
      fasta: /abs/path/reference.fna
      fai: /abs/path/reference.fna.fai
      bwa_index: /abs/path/bwa_index/reference
      wildtype_fasta: /abs/path/reference.fna
      genome_length: 5402402
      ppe_list: /abs/path/reference_matched_mask.txt
```

每条路线、每个参考和不同参数的分析使用独立 `results_dir`，避免复用旧结果。屏蔽文件是每行一个参考坐标的位置列表。NTM 路线要求在活动 reference profile 中显式提供 `ppe_list`。

## 运行

TB 核心分析：

```bash
snakemake --configfile config/config.yaml --cores 16 core_only \
  --config samplesheet=/abs/path/tb_samples.csv results_dir=results_tb
```

NTM 核心分析：

```bash
snakemake --configfile config/config.yaml config/ntm.yaml --cores 16 core_only
```

正式运行前可加 `-n` 检查任务。`--cores` 是整个流程可用 CPU 总额；`threads.bowtie2`、`threads.bwa`、`threads.samtools_sort` 控制相应任务线程。`threads.distance` 控制 snp-dists 线程，默认 1。

| target | 内容 |
|---|---|
| `foundation_only` | FASTQ 基础 QC 输入 |
| `core_only` | 比对、变异调用、核心 SNP 和距离 |
| `downstream_only` | 深度、插入片段及注释桥接输入 |
| `reports_only` | 已配置资源对应的报告 |
| `all` | 基础输入、核心分析、桥接、报告及启用的扩展 |

NTM 目前可使用上述核心路线；仓库内置的 TB 谱系、耐药标记和报告资源尚未完成通用 NTM 适配。NTM 配置示例因此使用 `core_only`。生成 NTM 注释或报告前需准备匹配资源。

## 输出

```text
<results_dir>/
├── core/
│   ├── diff_loci.txt
│   ├── merged.fasta
│   ├── wildtype.fasta
│   ├── core_snps.fadel-InvMisF5.bak.fa
│   ├── core_snps.fadel-InvMisF5.bak.loc
│   └── distance_matrix.txt
├── samples/<sample>/
│   ├── alignment_qc/          # sorted BAM；NTM 同时生成索引
│   ├── variant_analysis/      # .cns、.vars、.snp；NTM 额外生成 .forup
│   ├── report_inputs/
│   ├── figure/                # 启用报告时生成
│   └── table/                 # 启用报告时生成
└── logs/
    ├── rules/
    ├── workflow/
    └── runs/<run_id>/
```

`.vars` 在 TB 路线经过双链支持筛选，在 NTM 路线保留未经参考屏蔽的候选调用，口径不同。`.cns` 是共识调用表；无覆盖或不可判读位置不能解释为无突变。共享 core 结果仅保留一份，`.core_work/` 是内部中转目录。

SNP 距离由 `snp-dists -q -k -L -m -t -j N` 计算，再转换为三列上三角格式 `sample_1<TAB>sample_2<TAB>distance`。仅比较大写 A/C/G/T，忽略模糊碱基、gap 和小写屏蔽碱基。

可用 `--config run_id=batch01` 为状态文件指定批次；`logs/` 同时保留 latest 视图。

## 近缘性与预防分析扩展

仓库支持距离聚类、耐药标记筛查和调查提示：`cluster`、`drug-resistance`、`investigation-alerts`。聚类阈值由 `clustering.snp_threshold` 配置，应按研究菌种与样本背景确定；当前默认值不是通用 NTM 标准。内置耐药标记和分级资源面向 TB。

在兼容报告流程中，可用 `report.skip_extensions: true` 排除扩展任务。环境样本的低频候选、近缘簇和调查提示用于复核及预防分析；本项目尚未验证环境检测的灵敏度和特异度。

## 详细说明

- [NTM BWA-MEM 路线与输出定义](docs/NTM_bwamem.md)
- [基础配置](config/config.yaml)
- [测试](tests/)

项目名称已改为 **OpenNTM**。现有仓库 URL、历史脚本名称和旧命令别名保留兼容；GitHub 仓库仍位于 `yu323101/MTB-evo`。
