# 09_gate1b 任务 2（Q2c）：作者 DAPAR 流程分支

由 `02_author_branch.R`（主流程）生成；⑥ 节由 `02_author_compare.py` 填写。只报数。

## ① 输入与环境

- 输入：`data/pilot/proteome_pg_matrix_120.tsv`（10000 行 × 122 列；120 个 run 全部可由 07 `parse_design` 三个正则解析为 5 tp × 6 FR × 4 rep）。
- 并列表对照：`analysis/08_gate1/02_events.tsv`（表 E-C，由 `02_author_compare.py` 读取）。
- 作者脚本：`code/Protein contour/Phospho/SpatialProteoDynamics.github.io-v1.0/SpatialProteoDynamics-SpatialProteoDynamics.github.io-a6b8aac/DataProcessing/DAPAR_script_OSM.R`、`translocation_plots.R`（同目录）。
- R：R version 4.4.3 (2025-02-28)；Rscript = `<scratchpad>/envs/rdapar/bin/Rscript`（micromamba 环境，未装新包）。
- 包版本：limma 3.62.1；impute 1.80.0；metap 1.1；Biobase 2.66.0；DAPAR 1.38.0（未 `library`，见下）；dplyr 1.2.1、tidyr 1.3.2、stringr 1.6.0（仅被 DAPAR `limmaCompleteTest` 交叉核对调用）；MSnbase 2.32.0（环境内有，本脚本未用）。
- DAPAR 取法：`e <- new.env(); lazyLoad(file.path(.libPaths()[1], "DAPAR", "R", "DAPAR"), envir = e)`，再把 `e` 中闭包的环境设为 `e`（否则找不到 `pkgs.require` 等内部函数）。`library(DAPAR)` 因依赖 DAPARdata（代理 403）不可用。主流程底层直接调 limma / impute / stats；DAPAR 1.38 的 `LOESS`、`findMECBlock`、`reIntroduceMEC`、`getQuantile4Imp`、`limmaCompleteTest` 用于逐 fraction 交叉核对，`make.design.1`、`make.contrast` 用于构造 limma 设计矩阵与对比。
- DAPAR 源码行号引用 `<scratchpad>/dapar_1.38_src.txt`（下称 src:）；作者行号 D: = `DAPAR_script_OSM.R`，T: = `translocation_plots.R`。
- S0：原始 10000 行 → 去 `Genes` 空 4 行（D:10）→ 取第一个基因名（D:11-12；含 `;` 的 522 行）→ 去重复基因名 389 行（D:13，保留首次出现）→ **9607 行**进入流程。四蛋白对应蛋白组：GRB2 P62993、SHC1 P29353、CBL P22681、EGFR P00533。
- 非正值（按缺失）0 个；非空字符串解析失败 0 个。

## ② 处理顺序（实际使用的函数与参数）

| 步 | 作者行号 | DAPAR 1.38 源码行号 | 本轮实际函数与参数 |
|---|---|---|---|
| S0 | D:9-13 | — | `Genes != ""` → `gsub("^;","")` → `gsub(";.*","")` → `!duplicated()`；行 ID = 基因名 |
| S1 | D:21-28, D:40 | — | 每 fraction 20 列，按 [2min, 8min, 20min, 90min, CTRL] × Rep1–4 排（作者布局：4 个处理条件块在前，CTRL 在 17:20） |
| S2 | D:41 | — | `log2(x)`；≤ 0 按缺失（本数据 0 个） |
| S3 | D:42 | —（1.38 无 mvFilter） | 旧版 mvFilter 逻辑：保留至少一个条件中非缺失数 ≥ 3 的行（按旧版逻辑实现） |
| S4 | D:43-58 | — | 每对比 k：`KEEP` ⇔ CTRL 非缺失 ≥ 3 或 条件 k 非缺失 ≥ 3 |
| S5 | D:67 | wrapper.normalizeD src:25-46（LOESS 分支 src:42-43）；LOESS src:4-22（overall src:10-13） | `limma::normalizeCyclicLoess(x, method = "fast", span = 0.7)`（iterations = 3、weights = NULL 为库默认；NA 由 `loessFit` 排除） |
| S6 | D:68 | wrapper.impute.KNN src:49-73（逐条件 src:59-67；==0→NA src:68） | 逐条件 4 列 `impute::impute.knn(x, k = 15, rowmax = 0.99, colmax = 0.99, maxp = 1500, rng.seed = 1)`；之后 `x[x == 0] <- NA`（src:68） |
| S7a | D:69-70 | findMECBlock src:76-93；reIntroduceMEC src:96-105 | MEC = 某条件 4 个 rep 在 LOESS 后矩阵中全缺的 (行, 条件) 块；KNN 后把这些块重新置 NA |
| S7b | D:71-73 | getQuantile4Imp src:108-114（1.38 的 wrapper.impute.detQuant src:117-143 为 metacell 版，未用） | 每列 `stats::quantile(x, 0.025, na.rm = TRUE)`（type 7）× 1；旧版 `impute.detQuant(qData, values)`：每列剩余 NA ← 该列值（按旧版逻辑实现） |
| S8 | D:76-82 | limmaCompleteTest src:146-194（OnevsOne src:159-173）；make.design src:197-212；make.design.1 src:215-239；make.contrast src:242-301（contrast=1 src:277-286）；formatLimmaResult src:304-352 | 每对比 k：输入 = 条件 k 4 列 + CTRL 4 列，只取 KEEP 行；设计矩阵 = `make.design.1`（条件指示矩阵，非配对）；对比 = `make.contrast(contrast = 1)` → `( ConditionA )/ 1 -( ConditionB )/ 1`（tp − CTRL）；`limma::lmFit → limma::contrasts.fit → limma::eBayes`（默认参数）；`p.adjust(P, "BH", n = length(P))` 在该 fraction × 对比的 KEEP 行内 |
| S9 | D:89-107 | — | 6 个 fraction 的填补后 log2 矩阵按基因合并为 120 列（没过 S3 的块为 NA）；`rowSums(is.na) < 120` 过滤（D:105） |
| T1 | T:16-20 | — | 取该 tp 与 CTRL 各 24 列，`2^`（含填补值） |
| T2 | T:23-25 | — | 每 fraction 4 个 rep 的算术均值（线性） |
| T3 | T:30-31 | — | 份额 = fraction 均值 / 6 个 fraction 均值之和（tp、CTRL 分别） |
| T4 | T:33-34 | — | MS = \|份额_tp − 份额_CTRL\| 的 6 个值的均值 |
| T5 | T:36, T:38-47 | — | `match(max(x), x)`；`maxN`（原样）+ `match` |
| T6 | T:58-148 | — | fraction max1 / max2 在该 tp 的 limma 表中的 P_Value（作者第 3 列）；蛋白不在表 → 1 |
| T7 | T:150-153 | — | `metap::sumlog(c(p1, p2))$p` |
| T8 | T:157-171 | — | limegreen / dodgerblue / coral / gray，规则原样 |
| T9 | T:173 | — | `p.adjust(pval_combi, "BH")`，在全部可算蛋白上（见 ⑦ (g)） |
| T10 | T:179-184 | — | 标注条件 `FDR < 0.05 & MS > 0.1`（严格不等号，原样） |

## ③ 每 fraction 的过滤与填补统计

每 fraction 输入 9607 行 × 20 列。observed = LOESS 后仍为观测值的格；KNN = impute.knn 填补且未被 MEC 重置的格（含 impute.knn 内部退回列均值的格，见 KNN_colmean）；MEC_blocks = 某条件 4 个 rep 全缺的 (行, 条件) 块数；detQuant = 2.5% 分位数填补的格（= MEC 格 + impute.knn 输出为 0 被 src:68 置 NA 的格）。impute.knn 退回列均值的行（逐条件 5 次调用之和）分两类：rowmax 类 = 缺失比例 > rowmax = 0.99（即 4/4 缺失）的行，impute 发警告并用列均值填，随后 S7a 作为 MEC 重新置 NA、由 detQuant 填；newmiss 类 = KNN 计算后仍为 NA 的行（无警告），用列均值填，这些格留在最终矩阵中、来源记 KNN（KNN_colmean 列为其格数）。

| fraction | S3_kept | observed | KNN | KNN_colmean | MEC_blocks | detQuant | colmean_rows_rowmax | colmean_rows_newmiss | check |
|---|---|---|---|---|---|---|---|---|---|
| FR1 | 3872 | 70902 | 5138 | 9 | 350 | 1400 | 350 | 9 | 70902+5138+1400 = 77440 = 3872×20 ✓ |
| FR2 | 4485 | 80913 | 6747 | 131 | 510 | 2040 | 510 | 131 | 80913+6747+2040 = 89700 = 4485×20 ✓ |
| FR3 | 4041 | 69153 | 8939 | 399 | 682 | 2728 | 682 | 399 | 69153+8939+2728 = 80820 = 4041×20 ✓ |
| FR4 | 4612 | 84160 | 6404 | 22 | 419 | 1676 | 419 | 22 | 84160+6404+1676 = 92240 = 4612×20 ✓ |
| FR5 | 4064 | 73668 | 5684 | 121 | 482 | 1928 | 482 | 121 | 73668+5684+1928 = 81280 = 4064×20 ✓ |
| FR6 | 4207 | 74228 | 7840 | 11 | 518 | 2072 | 518 | 11 | 74228+7840+2072 = 84140 = 4207×20 ✓ |
| 合计 | 25281 | 453024 | 40752 | 693 | 2961 | 11844 | 2961 | 693 |  |

- rowmax 类退回列均值的行数 = MEC 块数（每 fraction 相等：是），且等于 impute.knn 警告 "N rows with more than 99 % entries missing; mean imputation used for these rows" 中 N 之和。
- detQuant 格 = 4 × MEC 块 + impute.knn 输出为 0 的格：FR1 1400 = 4×350 + 0；FR2 2040 = 4×510 + 0；FR3 2728 = 4×682 + 0；FR4 1676 = 4×419 + 0；FR5 1928 = 4×482 + 0；FR6 2072 = 4×518 + 0。impute.knn 输出为 0 的格每 fraction 为 0 / 0 / 0 / 0 / 0 / 0（只可能来自 knnimp.split 切出的 ≤ k 行小簇：库代码对其调 `meanimp(x[index, ])` 时缺失格已置 0、识别不到缺失；复刻计数小簇内缺失格 0 / 0 / 0 / 0 / 0 / 0），故 src:68 的 `== 0 → NA` 在本数据上无实际作用。
- S8 limma 行数（KEEP 行）：FR1 2min 3695 / 8min 3680 / 20min 3716 / 90min 3696；FR2 2min 4300 / 8min 4220 / 20min 4281 / 90min 4219；FR3 2min 3771 / 8min 3782 / 20min 3862 / 90min 3791；FR4 2min 4391 / 8min 4395 / 20min 4381 / 90min 4461；FR5 2min 3840 / 8min 3836 / 20min 3811 / 90min 3952；FR6 2min 3995 / 8min 3881 / 20min 3929 / 90min 3942。
- detQuant 值（每列 2.5% 分位数 × 1，log2）范围：FR1 19.780–20.226；FR2 19.600–20.128；FR3 18.221–19.477；FR4 19.020–19.950；FR5 18.738–19.825；FR6 18.562–19.323。120 个值见 `02_author_detquant_values.tsv`。
- S9：合并表（D:105 过滤后）7831 行 × 120 列（1776 行在 6 个 fraction 都没过 S3 而被去掉）；其中 120 列全非 NA（6 个 fraction 都过 S3）的 1382 行进入 T1–T10。
- S3 保留次数分布（9607 行中，在几个 fraction 过 S3）：0 个 1776 行；1 个 1892 行；2 个 1332 行；3 个 1315 行；4 个 1062 行；5 个 848 行；6 个 1382 行。

逐条件明细：

| fraction | condition | missing_cells | MEC_blocks | knn_meanimp_rows_rowmax | knn_warning_rows | knn_meanimp_rows_newmiss | knn_meanimp_cells_newmiss | knn_smallcluster_cells | knn_cluster_lines |
|---|---|---|---|---|---|---|---|---|---|
| FR1 | 2min | 1363 | 78 | 78 | 78 | 0 | 0 | 0 | 3 |
| FR1 | 8min | 1360 | 74 | 74 | 74 | 4 | 4 | 0 | 3 |
| FR1 | 20min | 1223 | 68 | 68 | 68 | 0 | 0 | 0 | 3 |
| FR1 | 90min | 1443 | 64 | 64 | 64 | 5 | 5 | 0 | 3 |
| FR1 | CTRL | 1149 | 66 | 66 | 66 | 0 | 0 | 0 | 3 |
| FR2 | 2min | 1568 | 92 | 92 | 92 | 16 | 16 | 0 | 5 |
| FR2 | 8min | 1980 | 150 | 150 | 150 | 1 | 1 | 0 | 3 |
| FR2 | 20min | 1740 | 89 | 89 | 89 | 70 | 70 | 0 | 5 |
| FR2 | 90min | 2036 | 128 | 128 | 128 | 44 | 44 | 0 | 4 |
| FR2 | CTRL | 1463 | 51 | 51 | 51 | 0 | 0 | 0 | 4 |
| FR3 | 2min | 2671 | 178 | 178 | 178 | 9 | 9 | 0 | 3 |
| FR3 | 8min | 2910 | 308 | 308 | 308 | 86 | 86 | 0 | 4 |
| FR3 | 20min | 2146 | 111 | 111 | 111 | 4 | 4 | 0 | 4 |
| FR3 | 90min | 1944 | 58 | 58 | 58 | 86 | 86 | 0 | 4 |
| FR3 | CTRL | 1996 | 27 | 27 | 27 | 214 | 214 | 0 | 3 |
| FR4 | 2min | 1651 | 105 | 105 | 105 | 2 | 2 | 0 | 4 |
| FR4 | 8min | 1744 | 126 | 126 | 126 | 0 | 0 | 0 | 4 |
| FR4 | 20min | 1856 | 109 | 109 | 109 | 0 | 0 | 0 | 5 |
| FR4 | 90min | 1162 | 43 | 43 | 43 | 0 | 0 | 0 | 4 |
| FR4 | CTRL | 1667 | 36 | 36 | 36 | 20 | 20 | 0 | 4 |
| FR5 | 2min | 1725 | 131 | 131 | 131 | 14 | 14 | 0 | 4 |
| FR5 | 8min | 1668 | 127 | 127 | 127 | 40 | 40 | 0 | 5 |
| FR5 | 20min | 1884 | 152 | 152 | 152 | 0 | 0 | 0 | 3 |
| FR5 | 90min | 1044 | 44 | 44 | 44 | 0 | 0 | 0 | 3 |
| FR5 | CTRL | 1291 | 28 | 28 | 28 | 67 | 67 | 0 | 4 |
| FR6 | 2min | 1636 | 53 | 53 | 53 | 8 | 8 | 0 | 4 |
| FR6 | 8min | 2456 | 177 | 177 | 177 | 2 | 2 | 0 | 3 |
| FR6 | 20min | 2064 | 129 | 129 | 129 | 0 | 0 | 0 | 3 |
| FR6 | 90min | 1964 | 108 | 108 | 108 | 1 | 1 | 0 | 3 |
| FR6 | CTRL | 1792 | 51 | 51 | 51 | 0 | 0 | 0 | 3 |

## ④ 四蛋白每 fraction 的状态

来源代码：O = observed，K = KNN（M = KNN 步骤中 impute.knn 退回列均值的格，tsv 中来源仍记 KNN），D = detQuant，- = filtered_out（该 fraction 没过 S3）；每格 4 个字符依次为 Rep1–Rep4。`n_obs` 为该 fraction 各条件原始非缺失数（CTRL/2/8/20/90）；过 S3 需至少一个条件 ≥ 3。KEEP 列为该蛋白是否进入该 fraction 对应对比（2/8/20/90min vs CTRL）的 limma 表。

**GRB2**（P62993）：原始非空格 90；observed 88、KNN 12（其中列均值退回 0）、detQuant 0、filtered_out 20（其中原始非空 2）。没过 S3 的 fraction：FR5。

| fraction | S3 | n_obs | CTRL | 2min | 8min | 20min | 90min | KEEP 2/8/20/90 |
|---|---|---|---|---|---|---|---|---|
| FR1 | 过 | 4/4/4/4/4 | OOOO | OOOO | OOOO | OOOO | OOOO | Y/Y/Y/Y |
| FR2 | 过 | 4/4/4/4/4 | OOOO | OOOO | OOOO | OOOO | OOOO | Y/Y/Y/Y |
| FR3 | 过 | 4/4/4/4/4 | OOOO | OOOO | OOOO | OOOO | OOOO | Y/Y/Y/Y |
| FR4 | 过 | 3/4/4/4/4 | OOKO | OOOO | OOOO | OOOO | OOOO | Y/Y/Y/Y |
| FR5 | 没过 | 1/0/1/0/0 | ---- | ---- | ---- | ---- | ---- | — |
| FR6 | 过 | 1/2/3/2/1 | KKKO | KKOO | OKOO | KOOK | KKOK | N/Y/N/N |

**SHC1**（P29353）：原始非空格 115；observed 115、KNN 5（其中列均值退回 0）、detQuant 0、filtered_out 0（其中原始非空 0）。没过 S3 的 fraction：无。

| fraction | S3 | n_obs | CTRL | 2min | 8min | 20min | 90min | KEEP 2/8/20/90 |
|---|---|---|---|---|---|---|---|---|
| FR1 | 过 | 4/4/4/4/4 | OOOO | OOOO | OOOO | OOOO | OOOO | Y/Y/Y/Y |
| FR2 | 过 | 4/4/4/4/4 | OOOO | OOOO | OOOO | OOOO | OOOO | Y/Y/Y/Y |
| FR3 | 过 | 4/4/4/4/4 | OOOO | OOOO | OOOO | OOOO | OOOO | Y/Y/Y/Y |
| FR4 | 过 | 4/4/4/4/4 | OOOO | OOOO | OOOO | OOOO | OOOO | Y/Y/Y/Y |
| FR5 | 过 | 1/4/4/4/3 | KKKO | OOOO | OOOO | OOOO | KOOO | Y/Y/Y/Y |
| FR6 | 过 | 4/4/4/4/3 | OOOO | OOOO | OOOO | OOOO | OOOK | Y/Y/Y/Y |

**CBL**（P22681）：原始非空格 59；observed 56、KNN 28（其中列均值退回 0）、detQuant 16、filtered_out 20（其中原始非空 3）。没过 S3 的 fraction：FR5。

| fraction | S3 | n_obs | CTRL | 2min | 8min | 20min | 90min | KEEP 2/8/20/90 |
|---|---|---|---|---|---|---|---|---|
| FR1 | 过 | 4/4/2/4/4 | OOOO | OOOO | KKOO | OOOO | OOOO | Y/Y/Y/Y |
| FR2 | 过 | 4/2/2/0/1 | OOOO | KOOK | KOKO | DDDD | KOKK | Y/Y/Y/Y |
| FR3 | 过 | 4/1/2/2/2 | OOOO | KOKK | KOKO | KOOK | OOKK | Y/Y/Y/Y |
| FR4 | 过 | 1/4/4/3/0 | KKKO | OOOO | OOOO | OOKO | DDDD | Y/Y/Y/N |
| FR5 | 没过 | 0/2/1/0/0 | ---- | ---- | ---- | ---- | ---- | — |
| FR6 | 过 | 0/1/3/2/0 | DDDD | OKKK | OKOO | KOOK | DDDD | N/Y/N/N |

**EGFR**（P00533）：原始非空格 80；observed 76、KNN 4（其中列均值退回 1）、detQuant 0、filtered_out 40（其中原始非空 4）。没过 S3 的 fraction：FR1、FR2。

| fraction | S3 | n_obs | CTRL | 2min | 8min | 20min | 90min | KEEP 2/8/20/90 |
|---|---|---|---|---|---|---|---|---|
| FR1 | 没过 | 0/0/0/0/0 | ---- | ---- | ---- | ---- | ---- | — |
| FR2 | 没过 | 2/0/0/2/0 | ---- | ---- | ---- | ---- | ---- | — |
| FR3 | 过 | 3/4/4/4/4 | OOOM | OOOO | OOOO | OOOO | OOOO | Y/Y/Y/Y |
| FR4 | 过 | 4/4/4/4/4 | OOOO | OOOO | OOOO | OOOO | OOOO | Y/Y/Y/Y |
| FR5 | 过 | 4/3/4/3/3 | OOOO | KOOO | OOOO | OOKO | KOOO | Y/Y/Y/Y |
| FR6 | 过 | 4/4/4/4/4 | OOOO | OOOO | OOOO | OOOO | OOOO | Y/Y/Y/Y |

## ⑤ 结果表（`02_author_branch.tsv`）

T1–T10 的蛋白集合 = 合并表中 6 个 fraction 都过 S3 的 1382 行（BH 的 n）。没过 S3 的蛋白（某 fraction 为 NA）在作者代码中无法计算（见 ⑦ (g)），本表留空。MS、份额 4 位小数；p 3 位有效数字；全精度值见 tsv。
四蛋白是否在这 1382 行内：GRB2 否（FR5 没过 S3）；SHC1 是；CBL 否（FR5 没过 S3）；EGFR 否（FR1、FR2 没过 S3）。

| gene | timepoint | MS | max1 | max2 | share_max1_tp | share_max1_ctrl | dir | p_max1 | p_max2 | pval_combi | FDR | color | MS>0.1 | FDR<0.05 | n_FDR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GRB2 | 2min |  |  |  |  |  |  |  |  |  |  |  |  |  | 1382 |
| GRB2 | 8min |  |  |  |  |  |  |  |  |  |  |  |  |  | 1382 |
| GRB2 | 20min |  |  |  |  |  |  |  |  |  |  |  |  |  | 1382 |
| GRB2 | 90min |  |  |  |  |  |  |  |  |  |  |  |  |  | 1382 |
| SHC1 | 2min | 0.0863 | FR4 | FR1 | 0.4139 | 0.1691 | up | 0.00872 | 0.0105 | 0.000939 | 0.144 | limegreen | no | no | 1382 |
| SHC1 | 8min | 0.0848 | FR1 | FR4 | 0.3603 | 0.5909 | down | 0.00537 | 0.0173 | 0.000955 | 0.102 | limegreen | no | no | 1382 |
| SHC1 | 20min | 0.0490 | FR1 | FR4 | 0.4759 | 0.5909 | down | 0.108 | 0.134 | 0.0759 | 0.404 | limegreen | no | no | 1382 |
| SHC1 | 90min | 0.0154 | FR1 | FR6 | 0.5571 | 0.5909 | down | 0.0138 | 0.0642 | 0.00712 | 0.322 | dodgerblue | no | no | 1382 |
| CBL | 2min |  |  |  |  |  |  |  |  |  |  |  |  |  | 1382 |
| CBL | 8min |  |  |  |  |  |  |  |  |  |  |  |  |  | 1382 |
| CBL | 20min |  |  |  |  |  |  |  |  |  |  |  |  |  | 1382 |
| CBL | 90min |  |  |  |  |  |  |  |  |  |  |  |  |  | 1382 |
| EGFR | 2min |  |  |  |  |  |  |  |  |  |  |  |  |  | 1382 |
| EGFR | 8min |  |  |  |  |  |  |  |  |  |  |  |  |  | 1382 |
| EGFR | 20min |  |  |  |  |  |  |  |  |  |  |  |  |  | 1382 |
| EGFR | 90min |  |  |  |  |  |  |  |  |  |  |  |  |  | 1382 |

全蛋白计数（每个 tp）：

| timepoint | n_proteins_in_FDR | n_MS_gt_0.1 | n_FDR_lt_0.05 | n_label_T181 | n_p1_eq_1 | n_p2_eq_1 |
|---|---|---|---|---|---|---|
| 2min | 1382 | 35 | 7 | 0 | 10 | 7 |
| 8min | 1382 | 33 | 3 | 0 | 9 | 15 |
| 20min | 1382 | 33 | 7 | 3 | 10 | 9 |
| 90min | 1382 | 29 | 5 | 1 | 4 | 7 |

（n_label_T181 = 满足作者 T:181 标注条件 `FDR < 0.05 & MS > 0.1` 的蛋白数；n_p1_eq_1 / n_p2_eq_1 = max1 / max2 fraction 的 limma 表中没有该蛋白而取 p = 1 的蛋白数。）

## ⑥ 并列表（`02_author_compare.tsv`）

<!-- SECTION6_BEGIN -->
MS_author_pipeline 来自本轮 `02_author_branch.tsv`；MS_gate1_raw / MS_gate1_centered / MS_perrep_mean_gate1_raw 原样抄 `analysis/08_gate1/02_events.tsv` 表 E-C 的 `MS_author`（raw / centered）与 `MS_perrep_mean`（raw）；空 = 原表为空串（gate1 未算；MS_author 为空的原因是该 tp 或 CTRL 有整个 fraction 4 个 rep 全缺，MS_perrep_mean 为空的原因是没有 rep 在两边 6 个 fraction 都有值，见 08_gate1/02_events.md ②）；`无` = gate1 E-C 无该基因（EGFR）。n_obs / n_imputed 为该蛋白在该 tp（或 CTRL）24 格（6 FR × 4 rep）中来源为 observed / KNN + detQuant 的格数；filtered_out 格两者都不计（最后一列为 filtered_out 格数，只在 md 中列出，tsv 无此列）。MS 4 位小数，全精度见 tsv。

| gene | timepoint | MS_author_pipeline | MS_gate1_raw | MS_gate1_centered | MS_perrep_mean_gate1_raw | n_obs_cells_tp | n_obs_cells_ctrl | n_imputed_cells_tp | n_imputed_cells_ctrl | filtered_out_tp/ctrl |
|---|---|---|---|---|---|---|---|---|---|---|
| GRB2 | 2min |  |  |  |  | 18 | 16 | 2 | 4 | 4/4 |
| GRB2 | 8min |  | 0.1298 | 0.1288 |  | 19 | 16 | 1 | 4 | 4/4 |
| GRB2 | 20min |  |  |  |  | 18 | 16 | 2 | 4 | 4/4 |
| GRB2 | 90min |  |  |  |  | 17 | 16 | 3 | 4 | 4/4 |
| SHC1 | 2min | 0.0863 | 0.0981 | 0.1013 | 0.0910 | 24 | 21 | 0 | 3 | 0/0 |
| SHC1 | 8min | 0.0848 | 0.0883 | 0.0951 | 0.0524 | 24 | 21 | 0 | 3 | 0/0 |
| SHC1 | 20min | 0.0490 | 0.0543 | 0.0679 | 0.0469 | 24 | 21 | 0 | 3 | 0/0 |
| SHC1 | 90min | 0.0154 | 0.0288 | 0.0371 |  | 22 | 21 | 2 | 3 | 0/0 |
| CBL | 2min |  |  |  |  | 12 | 13 | 8 | 7 | 4/4 |
| CBL | 8min |  |  |  |  | 13 | 13 | 7 | 7 | 4/4 |
| CBL | 20min |  |  |  |  | 11 | 13 | 9 | 7 | 4/4 |
| CBL | 90min |  |  |  |  | 7 | 13 | 13 | 7 | 4/4 |
| EGFR | 2min |  | 无 | 无 | 无 | 15 | 15 | 1 | 1 | 8/8 |
| EGFR | 8min |  | 无 | 无 | 无 | 16 | 15 | 0 | 1 | 8/8 |
| EGFR | 20min |  | 无 | 无 | 无 | 15 | 15 | 1 | 1 | 8/8 |
| EGFR | 90min |  | 无 | 无 | 无 | 15 | 15 | 1 | 1 | 8/8 |

### 自检（Python 部分）

| check | result | value |
|---|---|---|
| 02_author_compare.tsv：16 行，列 = 方案列，(gene, tp) 齐全 | PASS | 16 行 |
| gate1 三列与 08_gate1/02_events.tsv（E-C）逐格一致（EGFR 不在 E-C → 无） | PASS | E-C 基因 CBL,GRB2,SHC1；不一致 0 |
| MS_author_pipeline 与 02_author_branch.tsv 一致 | PASS | 不一致 0 |
| n_obs + n_imputed + filtered_out = 24（每 tp / CTRL） | PASS |  |
| 02_author_cells 值 = 02_author_merged_log2 对应格（Python 独立回读；filtered_out ⇔ 空） | PASS | 480 格，不一致 0 |
| 四蛋白 observed 格 = pg_matrix 原始非空格 − filtered_out 格（逐格集合相等） | PASS | GRB2 原始非空 90，observed 88，filtered_out 中原始非空 2；SHC1 原始非空 115，observed 115，filtered_out 中原始非空 0；CBL 原始非空 59，observed 56，filtered_out 中原始非空 3；EGFR 原始非空 80，observed 76，filtered_out 中原始非空 4 |
| T2–T4 用合并矩阵独立复算 MS：16 行一致（不在完整行内的蛋白 tsv 为空） | PASS | max\|ΔMS\| = 1.5e-15，max\|Δ份额\| = 4.1e-15，max1/max2 不一致 0 |
| T6–T9 独立复算 p_max1 / p_max2 / pval_combi / pval_combi_FDR 一致 | PASS | max\|Δp\|（pval_combi 为相对差）= 4.4e-15，max\|ΔFDR\| = 1.4e-15，BH 的 n = [1382] |
| 全局 Python 版本 = pandas 2.2.3 / numpy 1.26.4 / matplotlib 3.8.4（本脚本运行时） | PASS | pandas 2.2.3 / numpy 1.26.4 / matplotlib 3.8.4 |

- 全局 Python（本脚本运行时）：pandas 2.2.3 / numpy 1.26.4 / matplotlib 3.8.4。

<!-- SECTION6_END -->

## ⑦ 与作者实现的差别清单

- (a) 输入：本轮是 120 个单独 DIA-NN 搜库的 pg_matrix 合并表（`proteome_pg_matrix_120.tsv`，线性 PG 定量），不是作者一次 Spectronaut 合并搜库的蛋白报告（D:8-9）。S0 的去空、取第一个基因名、去重按 D:10-13 照做（本数据无以 `;` 开头的基因名）。
- (b) 条件：EGF 五个时间点（CTRL, 2min, 8min, 20min, 90min）；对比 4 个（2/8/20/90min vs CTRL）。作者 OSM 脚本为 5 个渗透压条件；translocation 作者只做 2min vs CTRL，本轮对四个 EGF 时间点各做一遍（把 T:16 的 `"2min"` 换成各 tp，其余代码不变）。
- (c) KNN 随机种子：作者 `rng.seed = sample(seq_len(1000), 1)`（src:65，随机，未设种子）；本轮每次调用固定 `rng.seed = 1`。每个条件的行数 > maxp = 1500，impute 走二分聚类（用随机起点），故种子影响 KNN 结果；作者的结果本身不可复现。
- (d) `mvFilter`、`impute.detQuant`：DAPAR 1.38 已无（改为 metacell 机制），按旧版逻辑实现：mvFilter(atLeastOneCond, th = 3) = 至少一个条件非缺失数 ≥ 3 保留；impute.detQuant(qData, values) = 每列剩余 NA 用该列的 shiftedImpVal 替换。
- (e) LOESS：以 DAPAR 1.38 源码为准（src:10-13：`limma::normalizeCyclicLoess(x = qData, method = "fast", span = span)`，span 默认 0.7），配 limma 3.62.1。作者 2020 年所用 DAPAR / limma 版本的 LOESS 实现是否完全相同无法核对。
- (f) T1 的 `2^` 输入是填补后的 log2 值（含 KNN 与 detQuant 填补值），与作者相同（作者 T:7 读的 merged_table 即 D:93-107 的填补后矩阵）。
- (g) 【方案未定义处】没过 S3 的 (蛋白, fraction) 块在合并表中为 NA（D:93 `merge(..., all = T)`）。作者 T:24-34 对含 NA 的行得到全 NA 的份额与 MS，T:47 的 `maxN` 随即报错（`sort.int: index 5 outside bounds`，本环境实测），即作者代码在含此类行的表上无法运行完。本轮把 T1–T10 限于 120 列全非 NA 的 1382 行（BH 的 n = 1382），其余蛋白的 MS / p / FDR 留空。四蛋白中受影响的见 ④、⑤。
- (h) T:56 作者用 `full_table$Gene.names` 作行名；OSM 版 D:107 写出的 merged_table 不含 `Gene.names` 列（基因名在行名中），作者 EGF 版 DAPAR 脚本不在仓库。本轮按方案 T6 用基因名（合并表行 ID）匹配 limma 表的行名（作者 csv 的 `X` 列）。
- (i) DAPAR 1.38 的 `wrapper.impute.KNN` 另有按 metacell 把 "Missing MEC" 置 NA 并更新 metacell（src:69-71）；本轮不建 metacell，用作者显式的 findMECBlock / reIntroduceMEC（D:69-70）实现同一效果。src:68 的 `== 0 → NA` 照 1.38 保留（本数据上 0 格，无实际作用，见 ③）；作者所用旧版 wrapper.impute.KNN 是否含此行无法离线核对。
- (j) limma：主流程直接调 `limma::lmFit → contrasts.fit → eBayes`，设计矩阵与对比字符串用 DAPAR 1.38 的 `make.design.1` / `make.contrast` 生成；与 DAPAR 1.38 `limmaCompleteTest` 整体调用的结果逐 fraction × 对比核对一致（见自检）。BH 照 D:81 在该 fraction × 对比的 KEEP 行内做。
- (k) 输出形式：limma 表按方案写成 tsv（列 gene, logFC, P_Value, BH），对应作者 D:82 的 csv（行名 = 基因名、logFC、P_Value、BH_pvalue）；T6 取的 "第 3 列" = P_Value（未校正）。

## ⑧ 跳过项

- D:83-85 的 .rnk 文件、D:89-91 / D:97-99 的每 fraction 条件均值 / 2^ / SD 表、D:101-102 的 `cbind(mv_keep, prot_imp_2)` 表：不影响 T1–T10，未写出。
- T:179-184 的 ggplot 火山图：未画；标注条件 `FDR < 0.05 & MS > 0.1` 以计数与 ⑤ 的 `MS_gt_0.1`、`FDR_lt_0.05` 列给出。
- `设计v1_final_2026-10-06.md`：缺，按 00_design 任务书执行。
- 无其它跳过；无新装 R 包。

## 自检（R 部分）

| check | result | value |
|---|---|---|
| 引用行号与作者脚本 / DAPAR 源码逐行一致 | PASS | DAPAR_script_OSM.R 29 处、translocation_plots.R 24 处、dapar_1.38_src.txt 17 处；不一致 0 |
| 每 fraction：observed + KNN + detQuant = S3 保留行数 × 20 | PASS | FR1 70902+5138+1400=3872×20；FR2 80913+6747+2040=4485×20；FR3 69153+8939+2728=4041×20；FR4 84160+6404+1676=4612×20；FR5 73668+5684+1928=4064×20；FR6 74228+7840+2072=4207×20 |
| impute.knn 退回列均值的行（缺失 > rowmax）= MEC 块数（每 fraction） | PASS | 与 impute.knn 警告中的行数一致 |
| detQuant 格 = MEC 格 + KNN 输出为 0 的格（src:68） | PASS |  |
| KNN_colmean 格（impute.knn 内 KNN 后仍为 NA → 列均值）全部来源记 KNN、格数 = 复刻计数 | PASS | FR1 9 / FR2 131 / FR3 399 / FR4 22 / FR5 121 / FR6 11 |
| impute.knn 逐行复刻（计数用）与 impute::impute.knn 输出完全相同 | PASS | 30 次调用 |
| S5 = DAPAR 1.38 LOESS(type='overall')（identical） | PASS | 6 fraction |
| S7a = DAPAR 1.38 findMECBlock / reIntroduceMEC（identical） | PASS | 6 fraction |
| S7b 分位数 = DAPAR 1.38 getQuantile4Imp（identical） | PASS | 120 列 |
| S8 直接 limma = DAPAR 1.38 limmaCompleteTest（logFC、P_Value） | PASS | max\|差\| = 0.0e+00 |
| observed 格的最终值 = LOESS 后观测值 | PASS | max\|差\| = 0.0e+00 |
| 四蛋白原始非空格数 = GRB2 90 / SHC1 115 / CBL 59 / EGFR 80 | PASS | GRB2 90 / SHC1 115 / CBL 59 / EGFR 80 |
| 四蛋白 observed 格数 = 原始非空格数 − 没过 S3 的 fraction 中的非空格数 | PASS | GRB2 88 = 90 − 2；SHC1 115 = 115 − 0；CBL 56 = 59 − 3；EGFR 76 = 80 − 4 |
| 02_author_cells 的值与 02_author_merged_log2 一致（filtered_out ⇔ 空） | PASS | 480 格，不一致 0 |
| branch 16 行；n_proteins_in_FDR = 合并表中 120 列全非 NA 的行数 | PASS | 1382 |
| pval_combi_FDR 的 BH 在全部可算蛋白上做（非四蛋白） | PASS | 2min 1382 / 8min 1382 / 20min 1382 / 90min 1382 |
| 作者 maxN 对含 NA 的行报错（确认 T1–T10 限于完整行的原因） | PASS | sort.int: index 5 outside bounds |
| 复跑一致：R 产物 md5 与上一次运行相同（rng.seed = 1） | PASS | 与上一次运行逐文件 md5 比较：28 / 28 相同 |

- 复跑：与上一次运行逐文件 md5 比较：28 / 28 相同。

交叉核对明细（每 fraction）：

| fraction | loess_identical | knn_trace_identical | findMEC_identical | reIntroduceMEC_identical | getQuantile4Imp_identical | limma_vs_limmaCompleteTest_maxdiff | observed_vs_loess_maxdiff |
|---|---|---|---|---|---|---|---|
| FR1 | TRUE | TRUE | TRUE | TRUE | TRUE | 0.0e+00 | 0.0e+00 |
| FR2 | TRUE | TRUE | TRUE | TRUE | TRUE | 0.0e+00 | 0.0e+00 |
| FR3 | TRUE | TRUE | TRUE | TRUE | TRUE | 0.0e+00 | 0.0e+00 |
| FR4 | TRUE | TRUE | TRUE | TRUE | TRUE | 0.0e+00 | 0.0e+00 |
| FR5 | TRUE | TRUE | TRUE | TRUE | TRUE | 0.0e+00 | 0.0e+00 |
| FR6 | TRUE | TRUE | TRUE | TRUE | TRUE | 0.0e+00 | 0.0e+00 |

