# EBV Type 1/2 Identification Pipeline — Design Report

## 1. 概述

本报告描述了一套从 RNA-seq 数据鉴定 EBV 感染型别（I 型或 II 型）的完整生物信息学 Pipeline 的设计原理、方法选择依据及合理性论证。

**输入**：EBV 感染人类细胞的 Illumina 短读长 RNA-seq FASTQ 数据
**输出**：EBV 型别判定（Type 1 / Type 2）+ 置信度 + 潜伏类型 + 重组检测 + EBV 转录组概览
**交付形式**：Snakemake 工作流 + Markdown 报告 + JSON 结果摘要

---

## 2. 背景：EBV I 型与 II 型的分子差异

EBV 根据EBNA2基因序列的差异被分为 I 型（EBV-1，又称 A 型）和 II 型（EBV-2，又称 B 型）[1]。两型 EBNA2 核苷酸序列同一性仅约 70%，蛋白质序列同一性仅约 54%，是 EBV 基因组中差异最大的区域 [1, 2]。

除 EBNA2 外，EBNA3 家族基因（EBNA3A、EBNA3B、EBNA3C）也存在连锁的型间差异，但差异程度较 EBNA2 小 [3]。近年全基因组研究还发现 BZLF1（启动子 Zp-V3 变体与 II 型严格关联）、BZLF2、BLLF1 等基因也存在型特异性变异位点 [4, 5]。

**关键生物学差异**：I 型 EBV 体外转化 B 细胞的效率远高于 II 型；II 型感染的 B 细胞溶裂期基因表达更高 [6]。

### 参考基因组

| 型别 | 参考株 | NCBI RefSeq | 基因组大小 |
|------|--------|-------------|-----------|
| EBV-1 | B95-8 | NC_007605.1 | 171,823 bp |
| EBV-2 | AG876 | NC_009334.1 | 172,764 bp |

这两个 RefSeq 参考序列是所有主要 EBV 分型研究使用的标准参考 [1, 4, 5, 7]。

---

## 3. Pipeline 总体架构

```
FASTQ → QC/trimming → 人类reads过滤 → 双参考比对 → 基因定量+潜伏评估 → 型别分类 → 重组检测 → 转录组概览 → 最终报告
```

### 数据流图

```
                    ┌─────────┐
  FASTQ ──► fastp ──┤ trimmed │
                    └────┬────┘
                         │
              ┌──────────▼──────────┐
              │  HISAT2 → GRCh38    │
              │  (人类reads过滤)     │
              └──────────┬──────────┘
                         │ unmapped reads
              ┌──────────┴──────────┐
              │                     │
     ┌────────▼────────┐  ┌────────▼────────┐
     │ HISAT2 → EBV-1  │  │ HISAT2 → EBV-2  │
     │ (NC_007605.1)   │  │ (NC_009334.1)   │
     └────────┬────────┘  └────────┬────────┘
              │                     │
     ┌────────┴─────────────────────┴────────┐
     │              分支处理                   │
     ├──────────────┬──────────────┬──────────┤
     │              │              │          │
     ▼              ▼              ▼          ▼
  featureCounts  bcftools      bedtools    pysam
  (基因定量)     (SNP calling) (覆盖度)    (基因组覆盖)
     │              │              │          │
     ▼              ▼              ▼          ▼
  潜伏评估       SNP一致性      覆盖度评分   转录组可视化
     │              │              │          │
     └──────────────┴──────┬───────┴──────────┘
                           │
                    ┌──────▼──────┐
                    │  型别分类    │
                    │  (综合评分)  │
                    └──────┬──────┘
                           │
                    ┌──────▼──────┐
                    │  重组检测    │
                    └──────┬──────┘
                           │
                    ┌──────▼──────┐
                    │  最终报告    │
                    └─────────────┘
```

---

## 4. 各步骤方法与合理性

### 4.1 QC 与 trimming（fastp）

**方法**：使用 fastp 对原始 FASTQ 进行质量修剪，参数为 Q20 质量截断、最短 30 nt、自动接头检测。

**合理性**：fastp 是 Illumina RNA-seq 数据 QC 的事实标准工具，速度快（多线程），集 QC 报告、接头去除、质量修剪于一体 [8]。低质量碱基和接头序列会干扰比对，尤其是 EBV reads 占比极低（通常 <1%），任何噪声都会放大。

### 4.2 人类 reads 过滤（HISAT2 → GRCh38）

**方法**：将修剪后的 reads 比对至人类参考基因组 GRCh38，提取未比对上的 reads（`samtools view -f 4`）用于后续 EBV 比对。

**合理性**：
- EBV 感染的人类细胞 RNA-seq 中，>99% 的 reads 来源于人类基因组。先过滤人类 reads 可大幅减少 EBV 比对时的噪声和计算量。
- 选择 HISAT2 而非 STAR 的原因是内存效率：HISAT2 的人类基因组索引约 8 GB RAM，而 STAR 需要约 30 GB [9]。对于在标准服务器上运行的 Pipeline，这是关键考量。
- 使用 `--very-sensitive` 参数确保最大比对灵敏度，避免将 EBV reads 误比对至人类基因组。

### 4.3 双参考 EBV 比对（HISAT2 → EBV-1 + EBV-2）

**方法**：将人类未比对上的 reads 分别比对至 EBV-1 参考基因组（NC_007605.1）和 EBV-2 参考基因组（NC_009334.1），生成两个 BAM 文件。

**合理性**：
- 这是本 Pipeline 的核心策略。由于 EBNA2 在两型之间有约 30% 的序列差异 [1, 2]，来自 I 型感染的 reads 在 EBV-1 参考上的比对率会显著高于 EBV-2 参考，反之亦然。
- 双参考比对是已发表的 EBV 基因分型 Pipeline 中使用的方法 [4, 5]。Blazquez 等 [4] 在对 278 个 EBV 基因组的分析中，正是通过将 reads 比对至两型参考基因组并比较覆盖度和变异数量来进行分型。
- 使用 `--no-spliced-alignment`：EBV 是 DNA 病毒，基因组仅 175 kb，不需要剪接比对。虽然 EBV 转录本存在剪接（如 EBNA2、EBNA3），但对于分型目的，我们关注的是外显子区域的 SNP 位点，非剪接比对已足够。
- 使用 `--very-sensitive` 确保在高度变异的 EBNA2 区域也能获得最佳比对。

### 4.4 基因定量与潜伏评估

**方法**：使用 featureCounts 对两个 EBV 比对分别进行基因水平定量，计算 TPM，然后根据潜伏基因的表达模式判断潜伏类型。

**潜伏分类逻辑**：
- **Latency III**：EBNA2 + EBNA3A/3B/3C + LMP1 + LMP2A 均表达（TPM > 10）
- **Latency II**：LMP1 和/或 LMP2A 表达，EBNA2 缺失（TPM < 1）
- **Latency I**：仅 EBNA1 + EBER1/2 表达
- **Lytic-active**：溶裂期基因占总 EBV 表达 > 50%

**合理性**：
- EBV 潜伏类型决定了哪些基因在 RNA-seq 中有 reads，从而决定了分型策略的可行性 [10, 11]。
- Latency III（如 LCL）中 EBNA2 高表达，EBNA2-based 分型最可靠。
- Latency I（如 Burkitt 淋巴瘤）中 EBNA2 不表达，需依赖 EBNA3 或全基因组 SNP 一致性。
- Pipeline 设计为自适应：先评估潜伏类型，再决定使用哪种分型策略，这解决了用户"潜伏类型不确定"的需求。

### 4.5 型别分类：覆盖度评分 + SNP 一致性

#### 4.5.1 覆盖度评分

**方法**：对每个型别鉴别区域（EBNA2、EBNA3A/3B/3C、BZLF1、BZLF2、BLLF1），比较 EBV-1 比对和 EBV-2 比对中的 read 计数，计算 log2 比值。正比值 → I 型；负比值 → II 型。

**合理性**：
- 型别鉴别区域是两型 EBV 基因组差异最大的区域 [1, 3, 4, 5]。在这些区域，来自 I 型感染的 reads 会更有效地比对至 EBV-1 参考，产生更高的 read 计数。
- Blazquez 等 [4] 的全基因组分析确认了 EBNA2、EBNA3 家族、BZLF1、BZLF2、BLLF1 是型间差异最显著的基因。
- 使用 log2 比值而非绝对计数，可归一化不同区域的表达水平差异。

#### 4.5.2 SNP 一致性评分

**方法**：
1. 在 EBV-1 比对上调用 SNPs（bcftools mpileup + call）
2. 对每个型特异性 SNP 位点（EBV-1 与 EBV-2 参考基因组差异的位置），检查调用的等位基因是否匹配 EBV-1 参考碱基或 EBV-2 参考碱基
3. 计算 I 型一致性比例和 II 型一致性比例
4. 一致性更高的型别即为感染型别

**型特异性 SNP 表的生成**：
- 使用 minimap2（asm10 预设）将 EBV-1 参考基因组（NC_007605.1）与 EBV-2 参考基因组（NC_009334.1）进行全基因组比对
- 从比对中提取所有 SNP 位点（共 1,538 个）
- 注释每个 SNP 所属的基因区域
- 其中型别鉴别基因内的 SNP：EBNA2 = 251 个，EBNA3A = 636 个，EBNA3BC = 910 个，BZLF1 = 32 个，BZLF2 = 10 个，BLLF1 = 45 个

**合理性**：
- SNP 一致性是最直接的型别证据：如果观察到的等位基因与 EBV-1 参考匹配，则病毒为 I 型 [4, 5]。
- 使用 minimap2 asm10 预设（允许 10% 序列差异）是必要的，因为 EBNA2 在两型间有约 30% 差异，更严格的预设（如 asm5）会在该区域中断比对。
- 1,538 个 SNP 远超 Pipeline 要求的最小 20 个可调用 SNP，确保即使在低覆盖度情况下也有足够的分型标记。
- bcftools 是 SNP calling 的标准工具，设置最小深度 3、最小碱基质量 20、最小比对质量 20 以过滤低质量变异调用。

#### 4.5.3 综合分类

**方法**：整合覆盖度评分和 SNP 一致性评分，给出最终型别判定和置信度。

**置信度分级**：
- **High**：覆盖度和 SNP 一致性均指向同一型别，SNP 一致性 ≥ 85%，可调用 SNP ≥ 20 个
- **Medium**：两种方法一致但数据量较低，或仅有一种方法可用
- **Low**：两种方法不一致（可能重组），或数据严重不足
- **Insufficient**：总 EBV reads < 1000，建议 DNA-based PCR

**合理性**：双证据整合提供了鲁棒性。覆盖度和 SNP 一致性是独立的证据线——覆盖度反映整体比对偏好，SNP 一致性反映特定位点的等位基因。两者一致时置信度最高；不一致时提示重组或混合感染。

### 4.6 重组/混合感染检测

**方法**：
1. 独立对 EBNA2 区域进行分型（仅使用 EBNA2 内的 SNPs 和覆盖度）
2. 独立对 EBNA3 区域进行分型（仅使用 EBNA3A/3B/3C 内的 SNPs 和覆盖度）
3. 检查 BZLF1 区域的型别模式（Zp-V3 变体是 II 型的标志 [4, 5]）
4. 若 EBNA2 型别 ≠ EBNA3 型别 → 标记为型间重组株
5. 若同一区域内同时检测到两种型别信号 → 标记为混合感染

**合理性**：
- 型间重组株在 EBV 中有充分记录 [4, 12, 13]。Midgley 等 [12] 在中国人群中发现了天然存在的型间重组株；Yao 等 [13] 在免疫抑制个体中分离到重组株。
- Blazquez 等 [4] 在 278 个 EBV 基因组中鉴定了 5 个潜在的型间重组株，其中 4 个为 EBV1-EBV2 重组（I 型 EBNA2 + II 型 EBNA3），1 个为 EBV2-EBV1 重组。
- BZLF1 启动子 Zp-V3 变体与 EBV-2 严格关联（100% 的 EBV-2 基因组携带 Zp-V3，而 EBV-1 主要携带 Zp-P）[4, 5]，提供了独立的型别标记。

### 4.7 EBV 转录组概览

**方法**：生成以下可视化和表格：
- EBV 基因 TPM 柱状图（按表达量排序，按潜伏/溶裂着色）
- 潜伏 vs 溶裂基因表达比例饼图
- EBV 基因组覆盖度图
- 型别鉴别区域覆盖度比较图（EBV-1 vs EBV-2 比对并排）
- 每基因表达量表（count、TPM、类别）

**合理性**：转录组概览帮助用户理解 EBV 在细胞系中的表达状态，验证潜伏类型判定，并直观检查型别鉴别区域的覆盖度差异。

### 4.8 最终报告

**方法**：将所有结果整合为一份 Markdown 报告，包含执行摘要、方法概述、比对统计、覆盖度评分、SNP 一致性、重组分析、潜伏评估、转录组概览（含图）、局限性和建议、参考文献。同时输出机器可读的 JSON 摘要。

---

## 5. 工具选择依据

| 工具 | 用途 | 选择理由 |
|------|------|----------|
| fastp | QC/trimming | Illumina RNA-seq QC 事实标准，速度快，集成度高 [8] |
| HISAT2 | 比对 | 内存效率高（8 GB vs STAR 的 30 GB），灵敏度好 [9] |
| samtools | BAM 处理 | 标准 BAM/BAM 操作工具 |
| bcftools | SNP calling | 标准 variant calling 工具 |
| bedtools | 覆盖度计算 | 标准基因组区间运算工具 |
| featureCounts | 基因定量 | 快速、准确的 RNA-seq 基因定量工具 [14] |
| minimap2 (mappy) | 参考基因组比对 | 快速全基因组比对，asm10 预设适合近缘基因组 |
| pysam | BAM 读取 | Python BAM 接口，用于覆盖度可视化 |
| Snakemake | 工作流 | 可复现、可扩展的生物信息学工作流管理系统 |

---

## 6. 计算资源估算

| 步骤 | 内存 | 时间（每样本） | 磁盘 |
|------|------|---------------|------|
| fastp QC | <1 GB | ~5 min | ~与输入相当 |
| HISAT2 → GRCh38 | ~8 GB | ~30 min（50M PE reads） | ~5 GB BAM |
| HISAT2 → EBV | <1 GB | <1 min | <100 MB |
| bcftools SNP calling | <1 GB | <5 min | <10 MB |
| bedtools coverage | <1 GB | <1 min | <1 MB |
| featureCounts | <1 GB | <1 min | <1 MB |
| 可视化 | <2 GB | <5 min | ~5 MB 图 |
| **总计** | **~8 GB** | **~40 min** | **~5 GB** |

**执行环境**：标准服务器（8 CPU / 32 GB RAM）可处理 1-3 个样本；多样本批量处理建议使用更高配置的工作节点。

---

## 7. 局限性

1. **RNA-seq 覆盖度依赖**：Pipeline 依赖 EBV 来源的 RNA-seq reads。若细胞系处于 Latency I（EBNA2 不表达），EBNA2-based 分型不可行，Pipeline 会回退至 EBNA3 或全基因组 SNP 一致性（较低置信度）。

2. **Poly-A 捕获偏差**：若建库使用 poly-A 捕获，非 poly-A 转录本（EBER1/2、部分 BART RNAs）将代表性不足。这不影响 EBNA2/EBNA3 的分型（它们是 poly-A 转录本），但可能影响潜伏类型判定的完整性。

3. **重组株检测局限**：Pipeline 可检测 EBNA2/EBNA3 型别不一致的重组株，但复杂的重组模式（如多个断点）可能需要全基因组测序才能完全解析。

4. **混合感染**：若细胞系同时携带 I 型和 II 型 EBV，Pipeline 会标记混合信号，但无法从 RNA-seq 单独确定各型别的比例。

5. **DNA-based 确认**：对于发表级结果，建议用 EBNA2 靶向 PCR + Sanger 测序确认 RNA-seq 的分型结果，这仍是 EBV 分型的金标准 [15]。

---

## 8. 测试与验收标准

1. **I 型阳性对照**：使用已知 EBV-1 细胞系（如 B95-8 来源的 LCL）的 RNA-seq → Pipeline 判定为 Type 1，高置信度
2. **II 型阳性对照**：使用已知 EBV-2 细胞系（如 AG876）的 RNA-seq → Pipeline 判定为 Type 2，高置信度
3. **低输入测试**：下采样至 10% reads → Pipeline 要么正确判定（中置信度），要么报告数据不足（无误判）
4. **重组株测试**：混合 50% EBV-1 EBNA2 reads + 50% EBV-2 EBNA3 reads → Pipeline 标记重组
5. **输出完整性**：所有图、JSON、Markdown 报告均无错误生成

---

## 9. 参考文献

[1] Tzellos S, Farrell PJ. Epstein-Barr Virus Sequence Variation—Biology and Disease. Pathogens. 2012;1(2):156-175. doi:10.3390/pathogens1020156

[2] Adldinger HK et al. A putative transforming gene of Jijoye virus differs from that of Epstein-Barr virus prototypes. Virology. 1985;141:221-234.

[3] Sample J et al. Epstein-Barr virus types 1 and 2 differ in their EBNA-3A, EBNA-3B, and EBNA-3C genes. J Virol. 1990;64(9):4084-4092.

[4] Blazquez AC et al. A Comparative Genomic Analysis of Epstein-Barr Virus Strains with a Focus on EBV2 Variability. Int J Mol Sci. 2025;26(6):2708. doi:10.3390/ijms26062708

[5] Correia S et al. Sequence Variation of Epstein-Barr Virus: Viral Types, Geography, Codon Usage, and Diseases. J Virol. 2018;92(20):e01132-18. doi:10.1128/JVI.01132-18

[6] Romero-Masters JC et al. B cells infected with Type 2 Epstein-Barr virus (EBV) have increased NFATc1/NFATc2 activity and enhanced lytic gene expression in comparison to Type 1 EBV infection. PLoS Pathog. 2020;16(2):e1008365.

[7] Palser AL et al. Genome Diversity of Epstein-Barr Virus from Multiple Tumor Types and Normal Infection. J Virol. 2015;89(10):5222-5237. doi:10.1128/JVI.03614-14

[8] Chen S et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34(17):i884-i890.

[9] Kim D et al. HISAT2: Graph-based alignment of next-generation sequencing reads to a population of genomes. Nat Biotechnol. 2019;37:907-915.

[10] Lin Z et al. Quantitative and Qualitative RNA-Seq-Based Evaluation of Epstein-Barr Virus Transcription in Type I Latency Burkitt's Lymphoma Cells. J Virol. 2010;84(24):13093-13101.

[11] O'Grady T et al. Analysis of EBV Transcription Using High-Throughput RNA Sequencing. Methods Mol Biol. 2017;1532:131-153.

[12] Midgley RS et al. Novel intertypic recombinants of Epstein-Barr virus in the Chinese population. J Virol. 2000;74(3):1544-1548.

[13] Yao QY et al. Isolation of intertypic recombinants of Epstein-Barr virus from T-cell-immunocompromised individuals. J Virol. 1996;70(8):4895-4903.

[14] Liao Y et al. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014;30(7):923-930.

[15] Lin JC et al. Precision of genotyping of Epstein-Barr virus by polymerase chain reaction using three gene loci (EBNA-2, EBNA-3C, and EBER). Blood. 1993;81(12):3372.
