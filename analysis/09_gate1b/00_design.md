# 09_gate1b 设计方案（advisor 写，subagent 执行）

沿用前几轮通用约定：读表 `sep='\t', dtype=str, keep_default_na=False`，空串 = 缺失；run 名解析用 07 `parse_design` 三个正则；只读 `data/`、`code/`；输出全部写 `analysis/09_gate1b/`；每个任务保存可复跑脚本；只报数，不解释生物学；缺输入写"缺 X，跳过"，不造数据；subagent 不做 git 提交；**任何 R / Python 包只装在独立环境里，不改全局 Python（pandas 2.2.3 / numpy 1.26.4 / matplotlib 3.8.4）**。

## 0. 输入现状（advisor 已核，2026-10-06）
| 文件 | 状态 |
|---|---|
| `data/pilot/gate1_missing_check.tsv` | 480 行 = 4 蛋白 × 120 run；列 `run, gene, n_precursors, min_PG_qvalue, min_precursor_qvalue`；`n_precursors = 0` 的行 113（CBL 49 / EGFR 33 / GRB2 27 / SHC1 4），这些行 q 值为空串。 |
| `analysis/gate1_proteins.tsv` | 344 行（上一轮已用；`data/pilot/gate1_proteins.tsv` 仍是 0 行）。四蛋白有值 run 数 GRB2 90 / SHC1 115 / CBL 59 / EGFR 80，与 `proteome_pg_matrix_120.tsv` 中同基因行的非空 run 完全一致（advisor 已核）。 |
| `data/pilot/proteome_pg_matrix_120.tsv` | 10000 行 × 122 列（`Protein.Group, Genes` + 120 run）；run 名全部可解析（5 tp × 6 FR × 4 rep）；值为线性 PG 定量（26434.5–5.8e10，中位 1.24e7），整体缺失率 60.0%；`Genes` 空 4 行、重复 24 个基因名；四蛋白各一行：GRB2 P62993、SHC1 P29353、CBL P22681、EGFR P00533。 |
| `设计v1_final_2026-10-06.md` | **未上传到 `analysis/`**：缺，按本任务书的 Q2b/Q2c 描述执行。 |
| 作者脚本 | `code/Protein contour/Phospho/SpatialProteoDynamics.github.io-v1.0/SpatialProteoDynamics-SpatialProteoDynamics.github.io-a6b8aac/DataProcessing/DAPAR_script_OSM.R`（107 行）、`translocation_plots.R`（185 行）。 |

**环境（advisor 已建）**：`/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad/envs/rdapar`（micromamba，R 4.4.3，bioconductor-dapar 1.38.0、limma 3.62.1、impute 1.80.0、MSnbase 2.32.0、r-metap 1.1、r-imp4p 1.3）。**`library(DAPAR)` 不能用**：它依赖 DAPARdata，而 bioconductor.org / galaxy depot / GitHub archive 都被代理 403 拦截。可用方式：`e <- new.env(); lazyLoad(file.path(.libPaths()[1], "DAPAR", "R", "DAPAR"), envir = e)` 直接取 DAPAR 1.38 的函数源码（advisor 已 deparse 到 `<scratchpad>/dapar_1.38_src.txt`，430 行），其底层只调用 limma / impute / stats。作者 2020 年用的 DAPAR（约 1.18–1.20）中的 `mvFilter`、`impute.detQuant` 在 1.38 已不存在（改为 metacell 机制），按下面 §1 给出的旧版等价逻辑在 R 里实现，并在 md 注明"按旧版逻辑实现"。R 脚本用 `<env>/bin/Rscript` 跑。Python 等价实现只作备用（sklearn / statsmodels / inmoose 的 wheel 已下载在 `<scratchpad>/pyfallback/`）。

---

## 1. 作者处理流程（逐步抄录，行号为仓库文件行号）

### 1.1 `DAPAR_script_OSM.R`（按 fraction 分别处理）
| 步 | 作者代码 | 内容 | 本轮的等价实现（R） |
|---|---|---|---|
| S0 读表与去重 | `:9-13` | 读 Spectronaut 蛋白报告；`PG.Genes != ""`（:10）；去掉开头 `;`（:11）；只留第一个基因名 `gsub(";.*","")`（:12）；`!duplicated(PG.Genes)`（:13）保留首次出现 | 对 `proteome_pg_matrix_120.tsv`：去掉 `Genes` 空的 4 行；`Genes` 按 `;` 取第一个；重复基因名保留首次出现（报去掉的行数）。行标识用 `Genes`（作者用基因名作 ID，:31-35）。 |
| S1 分 fraction | `:21-28`、`:40` | `fractions <- FR1..FR6`；`for j in 1:6` 取列名含 `FRj` 的样本（:24）；分组表（:26-28）按 Condition | 每个 fraction 20 列（5 条件 × 4 rep），Condition = timepoint（CTRL, 2min, 8min, 20min, 90min），rep 从 run 名取。 |
| S2 log2 | `:41` `createMSnset(..., logData=T)` | 线性值取 log2 | `log2(x)`；0 或非正值按缺失。 |
| S3 缺失过滤 | `:42` `mvFilter(prot, type="atLeastOneCond", th=3)` | 旧版 DAPAR `mvFilter`：保留"至少一个条件里非缺失值个数 ≥ 3"的蛋白行 | 旧版等价：`keep <- apply(sapply(conds, function(c) rowSums(!is.na(X[, cond==c])) >= 3), 1, any)`。报每个 fraction 保留行数。 |
| S4 对比保留标记 | `:44-58` | `mv`：每条件非缺失个数；`mv_keep[,k] = "KEEP"` 若 CTRL ≥ 3 **或** 条件 k ≥ 3（:48, :51, :54, :57） | 对四个对比（2min/8min/20min/90min vs CTRL）各算 KEEP 标记，用于 S8 limma 的行筛选。 |
| S5 LOESS 归一化 | `:67` `wrapper.normalizeD(prot_mv, method="LOESS", type="overall")` | DAPAR 1.38 `LOESS(qData, conds, type="overall", span=0.7)` = `limma::normalizeCyclicLoess(x = qData, method = "fast", span = 0.7)`（见 `dapar_1.38_src.txt` :3-23） | 直接调 `limma::normalizeCyclicLoess(X, method="fast", span=0.7)`（含 NA；limma 对 NA 的处理照库默认）。记录 limma 版本。 |
| S6 KNN 填补 | `:68` `wrapper.impute.KNN(prot_norm, 15)` | DAPAR 1.38 `wrapper.impute.KNN`：**按条件分别**对该条件的 4 列调 `impute::impute.knn(X[, ind], k = K, rowmax = 0.99, colmax = 0.99, maxp = 1500, rng.seed = sample(1:1000, 1))`（src :48-74）；K = 15 为邻居数（邻居 = 其他蛋白行） | 同参数逐条件调 `impute::impute.knn`；`rng.seed` 作者为随机，本轮固定 `rng.seed = 1` 并在 md 写明；impute.knn 对缺失率 > rowmax 的行会改用列均值填补，记录这类行数。 |
| S7 整组缺失（MEC）恢复 + 低值填补 | `:69-70` `findMECBlock(prot_norm)` / `reIntroduceMEC(prot_imp_KNN, MEC)`；`:71-73` `getQuantile4Imp(qData)$shiftedImpVal` → `impute.detQuant(qData, values)` | MEC = 某条件 4 个 rep 全缺的 (蛋白, 条件) 块（src :75-94），KNN 后把这些块重新置 NA（src :95-106）；`getQuantile4Imp(qdata, qval = 0.025, factor = 1)` = 每列 2.5% 分位数 × 1（src :107-115）；旧版 `impute.detQuant(qData, values)` = 每列剩余 NA 用该列的 shiftedImpVal 替换 | 同；报每列的 2.5% 分位数（24 个 = 6 FR × ... 不，每 fraction 20 列，共 120 个值）。 |
| S8 limma | `:76-82` | 每个对比 k：输入 = 条件 k 的 4 列 + CTRL 4 列（:77），只取 KEEP 行（:78）；`limmaCompleteTest(input, conditions, comp.type="OnevsOne")`（:80）= `lmFit` + `makeContrasts` + `contrasts.fit` + `eBayes`（src :145-195；design = 条件指示矩阵 `make.design.1`，**非配对**，src :214-240）；BH 校正（:81） | 同：每 fraction × 对比做 `lmFit(X, design) → contrasts.fit → eBayes`，取 logFC 与 P.Value，BH 在该 fraction × 对比的 KEEP 行内做。写出与作者同格式的 `EGF_prot_limma_COND<tp>vsC_FR<j>.csv` 风格表（logFC, P_Value, BH）。 |
| S9 合并 | `:89-107` | 条件均值（:89）、2^（:90）、SD（:91）；`merge_table` 合并各 fraction 填补后的 log2 值（:93-95）；`rowSums(is.na) < 120` 过滤（:105）；写 `merged_table`（:107） | 合并 6 个 fraction 填补后的 log2 矩阵为 120 列；没过 S3 的 (蛋白, fraction) 块为 NA。 |

### 1.2 `translocation_plots.R`（作者只做 2 min vs CTRL；本轮对四个 EGF 时间点各做一遍）
| 步 | 作者代码 | 内容 |
|---|---|---|
| T0 | `:7-13` | 读 merged_table；读 6 个 fraction 的 limma csv（2min vs C） |
| T1 | `:16-20` | 取 2min 列与 CTRL 列；`2^`（回线性） |
| T2 | `:23-25` | `FR` 标签（每 fraction 4 列）；**每 fraction 4 个 rep 取均值** |
| T3 | `:30-31` | 份额 = fraction 均值 / 6 个 fraction 均值之和（分别对 tp 和 CTRL） |
| T4 | `:33-34` | `MS = |份额_tp − 份额_CTRL|`（6 个）；**Movement Score = 6 个的均值** |
| T5 | `:36`, `:38-47` | `MS_max1` = 差最大的 fraction 编号（`match(max)`）；`MS_max2` = 第二大的（`maxN` + `match`，最大值并列时等于 max1） |
| T6 | `:58-148` | p1 = fraction max1 的 limma p（第 3 列 = P_Value；蛋白不在该表里 → p = 1，:64-65）；p2 = fraction max2 的同样取法 |
| T7 | `:150-153` | `metap::sumlog(c(p1, p2))$p`（Fisher 合并） |
| T8 | `:157-171` | 方向配色：max1/max2 落在 {1,2}↔{3,4} → limegreen（Cyt↔Mem）；{1,2}↔{5,6} → dodgerblue（Cyt↔Nuc）；{3,4}↔{5,6} → coral（Mem↔Nuc）；其他 gray |
| T9 | `:173` | `pval_combi_FDR = p.adjust(pval_combi, "BH")`（在全部蛋白上） |
| T10 | `:179-184` | 标注条件 `pval_combi_FDR < 0.05 & MS_2min_mean > 0.1` |

本轮与作者的已知差别（md 必须逐条写）：(a) 输入是 120 个单独 DIA-NN 搜库的 pg_matrix 合并表，不是一次 Spectronaut 合并搜库的报告；(b) 条件是 EGF 五个时间点，对比 4 个；(c) KNN 随机种子固定为 1；(d) `mvFilter`、`impute.detQuant` 按旧版逻辑实现；(e) DAPAR 1.38 的 `LOESS` 与作者版本是否相同只能以 1.38 源码为准，md 写明；(f) 2^ 的输入是填补后 log2 值（含填补值），与作者相同。

---

## 任务 1：缺失核查（Q2b）→ `01_missing_classes.md` + `01_missing_classes.tsv` + `01_missing_check.py`
### 定义
对每个 (gene ∈ {GRB2, SHC1, CBL, EGFR}, run ∈ 120)：
- **A** = `analysis/gate1_proteins.tsv` 中该 (gene, run) 有 quantity（= pg_matrix 有值）。
- **B** = 非 A 且 `gate1_missing_check.tsv` 的 `n_precursors > 0`（parquet 有前体但 pg_matrix 无值）；记 `min_PG_qvalue`、`min_precursor_qvalue`、`n_precursors`。
- **C** = 非 A 且 `n_precursors == 0`。
- 另核：A 类行的 `n_precursors` 是否都 > 0（advisor：0 例外）；B 类再按 `min_PG_qvalue ≤ 0.01` / `> 0.01` 拆分。
### 输出
- `01_missing_classes.tsv`（480 行）：`gene, run, timepoint, fraction, rep, class, n_precursors, min_PG_qvalue, min_precursor_qvalue, pg_quantity`。
- `01_missing_classes.md`：① 输入 ② 定义 ③ 每蛋白 A/B/C 计数（合计 120）④ 每蛋白 × fraction 的 A/B/C 表（6 行 × 3 列 + 行合计 20）⑤ B 类清单（gene, run, n_precursors, min_PG_qvalue, min_precursor_qvalue）及 q ≤ 0.01 / > 0.01 的计数 ⑥ 跳过项。
### 验收（advisor 参考）
A/B/C：GRB2 90/3/27；SHC1 115/1/4；CBL 59/12/49；EGFR 80/7/33；每蛋白合计 120；B 类中 q ≤ 0.01：CBL 1，其余 0。每 fraction 行合计 = 20。

---

## 任务 2：作者流程分支（Q2c）→ `02_author_branch.md` + `02_author_branch.tsv` + `02_author_cells.tsv` + `02_author_compare.tsv` + `02_author_branch.R`（主流程，R）+ `02_author_compare.py`（并列表，Python）
### 执行
按 §1.1 S0–S9 在完整矩阵上跑（全部约 10000 行蛋白都进流程，limma 和 BH 在全体蛋白上做，不只四个蛋白），再按 §1.2 T1–T10 对四个 EGF 时间点各算一遍，取出 GRB2 / SHC1 / CBL / EGFR。参数全部固定为作者值：`th = 3`、`span = 0.7`、`K = 15, rowmax = 0.99, colmax = 0.99, maxp = 1500`、`qval = 0.025, factor = 1`、BH。**不得为了结果调整任何参数。**
### 输出
- `02_author_branch.tsv`（4 蛋白 × 4 tp = 16 行）：`gene, timepoint, MS_author_pipeline, MS_max1, MS_max2, share_max1_tp, share_max1_ctrl, direction_max1(up/down), p_max1, p_max2, pval_combi, pval_combi_FDR, color_category(limegreen/dodgerblue/coral/gray), MS_gt_0.1, FDR_lt_0.05, n_proteins_in_FDR`。
- `02_author_cells.tsv`（4 蛋白 × 120 格 = 480 行）：`gene, timepoint, rep, fraction, value_log2_final, source(observed / KNN / detQuant / filtered_out)`。`filtered_out` = 该蛋白在该 fraction 没过 S3（整 fraction 无值进入后续）。
- `02_author_compare.tsv`：四蛋白 × 4 tp 的并列表：`gene, timepoint, MS_author_pipeline, MS_gate1_raw, MS_gate1_centered, MS_perrep_mean_gate1_raw, n_obs_cells_tp, n_obs_cells_ctrl, n_imputed_cells_tp, n_imputed_cells_ctrl`，其中 gate1 的值直接抄 `analysis/08_gate1/02_events.tsv`（表 E-C，EGFR 在 gate1 没有 E-C，填"无"）。
- 中间产物（便于核对，均写入 `analysis/09_gate1b/`）：`02_author_limma_FR<j>_<tp>vsCTRL.tsv`（24 个；列 gene, logFC, P_Value, BH），`02_author_merged_log2.tsv`（合并后的 120 列矩阵，全部蛋白），`02_author_detquant_values.tsv`（120 列的 2.5% 分位数）。
- `02_author_branch.md`：① 输入与环境（R、各包版本、lazyLoad 方式）② §1 的处理顺序表（可直接引用 00_design §1，但每步要写实际用的函数与参数）③ 每 fraction 的过滤与填补统计：S3 保留行数、KNN 填补格数、MEC 块数、detQuant 填补格数、impute.knn 退回列均值的行数 ④ 四蛋白每 fraction 的状态（过/没过 S3；每 (tp, rep) 的 source）⑤ 结果表（02_author_branch.tsv 渲染）⑥ 并列表 ⑦ 与作者实现的差别清单（§1 末尾 (a)–(f) 逐条确认）⑧ 跳过项。
### 验收
1. 处理顺序表每步带作者脚本行号与 DAPAR 源码行号。
2. 填补格数报出：每 fraction 的 KNN 格数 + detQuant 格数 + 观测格数 = 保留行数 × 20。
3. 四蛋白 480 格的 source 与 `02_author_merged_log2.tsv` 一致（observed 格的 log2 值 = LOESS 后的观测值，非 observed 格有填补值）；四蛋白在 pg_matrix 的原始非空格数（GRB2 90 / SHC1 115 / CBL 59 / EGFR 80）= observed 格数（若该 fraction 没过 S3，则这些观测格记 filtered_out，md 说明）。
4. 并列表 16 行齐全；gate1 的值与 `08_gate1/02_events.tsv` 一致。
5. advisor 抽查：对 1 个 fraction 用独立 R 代码复算 S2–S7 后的矩阵（同 seed），四蛋白的 20 格一致（1e-6）；T2–T4 用合并矩阵复算 MS，16 行一致。
6. 全局 Python 版本前后不变。
