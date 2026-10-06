# 09 验收记录（advisor，gate1b）

方案：`00_design.md`。执行：两个 Opus subagent 并行（任务 1 Python；任务 2 R 主流程 + Python 并列表）。两个任务一轮通过，没有重派。验收方法同前：advisor 在分派前后用独立实现（Python 对任务 1；一套独立写的 R 流程对任务 2，同参数同 seed）算参考值，逐格比对。

## 输入与环境事实
- `设计v1_final_2026-10-06.md` 未上传到 `analysis/`：缺，按任务书的 Q2b/Q2c 描述执行。
- 三份新输入来自分支 `claude/determined-goldberg-atmb8n` commit 0e66840，已合并（该提交同时把 `data/pilot/gate1_proteins.tsv` 更新为 344 行，与 `analysis/gate1_proteins.tsv` md5 相同）。
- 独立 R 环境：micromamba，`<scratchpad>/envs/rdapar`，R 4.4.3、bioconductor-dapar 1.38.0、limma 3.62.1、impute 1.80.0、metap 1.1。`library(DAPAR)` 因依赖 DAPARdata 而失败，DAPARdata 无法获取（bioconductor.org、galaxy depot、GitHub archive 均被代理 403）；改用 `lazyLoad` 直接取 DAPAR 1.38 的函数源码，底层调用 limma / impute / stats。全局 Python 版本前后不变（2.2.3 / 1.26.4 / 3.8.4）。

## 任务 1 缺失核查（Q2b）— 合格
文件：`01_missing_classes.tsv`（480 行）、`01_missing_classes.md`、`01_missing_check.py`。
验收：每蛋白 A+B+C = 120；每蛋白 × fraction 的 A/B/C 与 advisor 参考逐格一致（0 处不符）；A 类 run 集合与 pg_matrix 四基因行的非空 run 集合完全一致；A 类中 n_precursors = 0 的行 0。

| gene | A（pg_matrix 有值） | B（无值但有前体） | C（无前体） |
|---|---|---|---|
| GRB2 | 90 | 3 | 27 |
| SHC1 | 115 | 1 | 4 |
| CBL | 59 | 12 | 49 |
| EGFR | 80 | 7 | 33 |

B 类 23 行的 n_precursors 全部 = 1；min_PG_qvalue ≤ 0.01 只有 1 例（CBL 8min FR1 Rep1，0.0093），其余 22 例 > 0.01；min_precursor_qvalue 全部 < 0.01。按 fraction：EGFR 的 C 集中在 FR1（19/20）、FR2（12/20）；GRB2 集中在 FR5（18/20）；CBL 在 FR5（16）、FR6（10）、FR2（9）、FR4（8）。

## 任务 2 作者流程分支（Q2c）— 合格
文件：`02_author_branch.R`、`02_author_compare.py`、`02_author_branch.tsv`（16 行）、`02_author_cells.tsv`（480 行）、`02_author_compare.tsv`（16 行）、`02_author_limma_FR<j>_<tp>vsCTRL.tsv`（24 个）、`02_author_merged_log2.tsv`（7831 行 × 120 列）、`02_author_detquant_values.tsv`（120 行）、`02_author_branch.md`。
验收：
- 处理顺序表每步带作者脚本行号（D: 29 处、T: 24 处）和 DAPAR 1.38 源码行号（17 处）；参数全部为作者值（th=3；LOESS fast span 0.7；KNN k=15 rowmax=colmax=0.99 maxp=1500，seed 固定 1；detQuant 2.5% 分位数 × 1；limma 非配对 OnevsOne + BH；sumlog）。
- 填补格数：每 fraction observed + KNN + detQuant = 保留行数 × 20（6 个 fraction 全部成立），与 advisor 的独立 R 流程逐数一致：

| FR | S3 保留行 | observed | KNN | MEC 块 | detQuant |
|---|---|---|---|---|---|
| FR1 | 3872 | 70902 | 5138 | 350 | 1400 |
| FR2 | 4485 | 80913 | 6747 | 510 | 2040 |
| FR3 | 4041 | 69153 | 8939 | 682 | 2728 |
| FR4 | 4612 | 84160 | 6404 | 419 | 1676 |
| FR5 | 4064 | 73668 | 5684 | 482 | 1928 |
| FR6 | 4207 | 74228 | 7840 | 518 | 2072 |

- 四蛋白 480 格的 value 与 source 与 advisor 参考逐格一致（0 处不符）：observed 335、KNN 49、detQuant 16、filtered_out 80。GRB2 在 FR5、CBL 在 FR5、EGFR 在 FR1 和 FR2 没过 S3（该条件内 ≥3 个非缺失的要求），这些 fraction 的格记 filtered_out（GRB2 20、CBL 20、EGFR 40）。
- 16 行结果：SHC1 四个时间点可算，MS / pval_combi / FDR 与 advisor 一致到全部有效位；GRB2、CBL、EGFR 四个时间点全部留空（见下）。并列表 16 行齐全，gate1 列与 `08_gate1/02_events.tsv` 一致。

| gene | tp | MS（作者流程） | FDR | max1→max2 | gate1 raw | gate1 centered |
|---|---|---|---|---|---|---|
| SHC1 | 2min | 0.0863 | 0.144 | FR4, FR1 | 0.0981 | 0.1013 |
| SHC1 | 8min | 0.0848 | 0.102 | FR1, FR4 | 0.0883 | 0.0951 |
| SHC1 | 20min | 0.0490 | 0.404 | FR1, FR4 | 0.0543 | 0.0679 |
| SHC1 | 90min | 0.0154 | 0.322 | FR1, FR6 | 0.0288 | 0.0371 |
| GRB2 / CBL / EGFR | 全部 | 空（某 fraction 没过 S3） | 空 | 空 | GRB2 8min 0.1298 / 0.1288，其余空 | |

作者标注条件（MS > 0.1 且 FDR < 0.05）在全部 1382 个六 fraction 齐全的蛋白里满足的数量：2min 0、8min 0、20min 3、90min 1；四个锚点蛋白均不满足。

## 改过什么
- 环境：原计划 `library(DAPAR)`，改为 lazyLoad 取函数（§0 已记录原因）；`mvFilter`、`impute.detQuant` 按旧版逻辑实现。
- subagent 自定的一处（advisor 认可，与 advisor 的独立实现相同）：某 fraction 没过 S3 的蛋白在合并表里为 NA，作者 T:47 的 `maxN` 在这种行上报错，因此 T1–T10 只在 120 列全非 NA 的 1382 行上做，其余蛋白 MS / p / FDR 留空；BH 的 n = 1382。
- 无重派；输出未被 advisor 手改。

## 没解决 / 需要知道的问题
1. **作者流程对三个锚点蛋白无输出**：GRB2（FR5）、CBL（FR5）、EGFR（FR1、FR2）在作者的 `mvFilter(atLeastOneCond, th=3)` 下被整个 fraction 剔除，Movement Score 在作者定义下不可算。若要给它们一个数，需要新的定义（如只在可用 fraction 上算，或缺 fraction 份额记 0）；本轮按约束未改参数、未加定义。
2. **缺失性质**：四蛋白无值的 136 个 (蛋白, run) 里 113 个 parquet 也没有前体（C 类），23 个有且只有 1 个前体但 PG q 值 > 0.01（22 例）或 ≤ 0.01（1 例）。
3. **与作者的实现差别**（`02_author_branch.md` ⑦ 共 11 条）：输入是 120 个单独搜库的合并矩阵；KNN seed 固定（行数 > maxp 时 impute.knn 走随机二分聚类，seed 影响结果）；LOESS 以 DAPAR 1.38 + limma 3.62.1 为准；`mvFilter`/`impute.detQuant` 旧版逻辑；作者 EGF 版脚本不在仓库，`Gene.names` 列用基因名代替；limma 表写成 tsv。
4. impute.knn 中因缺失率 > rowmax 退回列均值的行与 MEC 块一致，随后被置回 NA 由 detQuant 填；KNN 后仍为 NA 而退回列均值的格（FR1 9 … FR3 399）来源记为 KNN。
5. 合并表只写了作者 D:105 过滤后的 7831 行；6 个 fraction 都没过 S3 的 1776 行未写入。
