# 08_gate1 设计方案（advisor 写，subagent 执行）

沿用 `analysis/06_pilot/00_design.md` 的通用约定（读表 `sep='\t', dtype=str, keep_default_na=False`，空串 = 缺失；run 名解析用 07 `parse_design` 的三个正则 `_(2min|8min|20min|90min|CTRL)_`、`_(FR\d)_`、`_(Rep\d)`；只读 `data/`、`code/`；输出全部写 `analysis/08_gate1/`；每个任务保存可复跑脚本；只报数和方向计数，不解释生物学；缺输入写"缺 X，跳过"，不造数据；subagent 不做 git 提交；不装包，不改全局环境）。

## 0. 输入现状（advisor 已核，2026-10-06）
| 文件 | 状态 |
|---|---|
| `data/pilot/gate1_proteins.tsv` | **只有表头，0 行数据**（33 字节：`run\tgenes\tprotein_group\tquantity`）。来自分支 `claude/determined-goldberg-atmb8n` commit bee801a "Add gate1 protein quantities"，`analysis/gate1_proteins.tsv` 是同一文件的副本（md5 相同）。已合并进当前分支。 |
| `analysis/07_pilot2/01_run_factors.tsv` | 240 行，proteome 120 行可用。列 `layer, run, timepoint, fraction, rep, median_log2, a0, delta, a, n_equations, max_residual`。`a` = 07 实际使用的 log2 加性因子（`num = raw × 2^a`；`a = anchor − m_07`，其中 `m_07` 是 07 实际减掉的 run 中位数，`anchor` 是 layer 常数）。 |
| 作者代码 | `code/Protein contour/Phospho/SpatialProteoDynamics.github.io-v1.0/SpatialProteoDynamics-SpatialProteoDynamics.github.io-a6b8aac/DataProcessing/translocation_plots.R`（185 行）、`DAPAR_script_OSM.R`（107 行）。 |

**结论：三个任务的全部数值结果本轮都无法产生。** 执行要求：脚本按本方案写好并对着真实文件跑；真实文件 0 行时，每个输出文件照常生成，表只含表头，md 的每个结果节写"缺数据：`gate1_proteins.tsv` 0 行，跳过"，并报"四个蛋白各在 0/120 个 run 测到"。不得用任何合成数据写进 `analysis/08_gate1/`；允许在 scratchpad 用明显的合成夹具做代码冒烟测试，只在 md 的"跳过项"里写一句"脚本在 scratchpad 合成夹具上跑通"，不写任何数字。数据到位后直接复跑脚本即可。

## 1. 作者 Movement Score 的精确定义（原文照抄，行号为仓库文件行号）

`translocation_plots.R`：
```r
19  prot_2min<-2^prot_2min
20  prot_ctrl<-2^prot_ctrl
23  FR<-c(rep("FR1",4), rep("FR2", 4),  rep("FR3", 4),  rep("FR4", 4),  rep("FR5", 4),  rep("FR6", 4))
24  prot_2min_mean<-t(apply(prot_2min, 1, function(x) tapply(x,FR,function(x) {mean(x)})))
25  prot_CTRL_mean<-t(apply(prot_ctrl, 1, function(x) tapply(x,FR,function(x) {mean(x)})))
30  prot_2min_scaled<-t(apply(prot_2min_mean_raw, 1, function(x) {x/sum(x)}))
31  prot_Ctrl_scaled<-t(apply(prot_CTRL_mean_raw, 1, function(x) {x/sum(x)}))
33  MS_2min<-abs(prot_2min_scaled-prot_Ctrl_scaled)
34  MS_2min_mean<-apply(MS_2min, 1, function(x) mean(x))
36  MS_max1<-(apply(MS_2min, 1, function(x) {match(max(x),x)}))
47  MS_max2<-(apply(MS_2min, 1, function(x) {match(maxN(x),x)}))
181   geom_text_repel(aes(label=Gene), data=temp_table[(temp_table$pval_combi_FDR<0.05) & (temp_table$MS_2min_mean>0.1),])
```
即：① log2 值还原为线性（:19-20）；② 每个 fraction 内 4 个重复取**均值**（:23-25），得 6 个 fraction 均值；③ 份额 = 该 fraction 均值 / 6 个均值之和（:30-31）；④ 每个 fraction 的 |份额_tp − 份额_CTRL|（:33）；⑤ **Movement Score = 这 6 个绝对差的均值**（:34）；⑥ 差最大和次大的 fraction 记为 MS_max1 / MS_max2（:36, :47），用于方向配色（:157-171）和 limma p 值合并（:58-153）；⑦ 图上标注阈值 `MS_2min_mean > 0.1` 且 `pval_combi_FDR < 0.05`（:181）。作者只对 2 min vs CTRL 算。
作者输入的 `merged_table_HeLa_EGF_prot.txt` 来自 `DAPAR_script_OSM.R` 的流程（按 fraction 分别处理）：log2（:41 `logData=T`）→ 缺失过滤 `atLeastOneCond, th=3`（:42）→ **LOESS overall 归一化**（:67）→ KNN(15) 填补 + MEC + detQuant（:68-73）→ limma（:76-82）→ 合并（:93-107）。本轮不复现 DAPAR 的归一化和填补，只用 raw 与 centered 两版。

## 2. 共用定义
- 蛋白集：GRB2、SHC1、CBL、EGFR。匹配 `genes` 列（分号分隔的多基因按 `;` 拆开后任一等于目标基因即匹配，记录原始 `genes` 字符串）。
- 一个基因对应多个 `protein_group`：全部保留；`primary` = 出现 run 数最多的 group（并列取 `protein_group` 字符串排序第一，并记录并列）；其余为 `secondary`，写附表。
- run 解析：`timepoint ∈ {CTRL,2min,8min,20min,90min}`，`fraction ∈ FR1..FR6`，`rep ∈ Rep1..Rep4`；期望 120 个 proteome run；解析失败的 run 计数并列出。
- `quantity` 转 float；空串/非数视为缺失。
- 两版归一化（2026-10-06 确认）：`raw` = quantity；`centered` = quantity × 2^a（a 取 `01_run_factors.tsv` 中 layer=proteome、同 run 的 `a`；等价于 quantity / 2^(m_07 − anchor)，anchor 为常数）。run 名匹配用 (timepoint, fraction, rep) 三元组。
- 份额：对每个 (gene, protein_group, timepoint, rep)：`n_fractions_present` = 6 个 fraction 中有值的个数；`share_f = q_f / Σ_{present} q_f`；`Cyt = share_FR1 + share_FR2`、`Mem = share_FR3 + share_FR4`、`Nuc = share_FR5 + share_FR6`；某区室两个 fraction 都缺 → 该区室格留空；`n_fractions_present < 6` 的行在 md 单列"缺 fraction 清单"（gene, tp, rep, 缺哪些）；不填补。
- 区室均值/标准差：按 timepoint 对 4 个重复的份额取 `mean`、`sd(ddof=1)`，n 不足 4 时如实写 n。
- Δ份额 = `mean_tp − mean_CTRL`。
- 方向计数（2min、8min 各一次）：按 rep 编号配对，`Mem_tp,r − Mem_CTRL,r > 0` 且 `Cyt_tp,r − Cyt_CTRL,r < 0` 的 rep 数，写成 `n/4`（分母 = 两边都有份额的 rep 数，若不足 4 写实际分母并注明）；另报"仅 Mem 上升"、"仅 Cyt 下降"的 n/4。
- Movement Score（作者定义，§1）：对每个 timepoint：fraction 均值 `m_f,tp = mean over reps of q_f`（线性值；缺失的 rep 不计入均值，与作者在填补后无缺失的前提不同，需在 md 注明）；`s_f,tp = m_f,tp / Σ_f m_f,tp`；`MS_tp = mean_f |s_f,tp − s_f,CTRL|`；同时报 `MS_max1`、`MS_max2`（fraction 编号）。阈值按作者 :181 的严格 `> 0.1`（2026-10-06 确认；原稿 ≥ 0.1 作废）。raw、centered 各算。另加一个对照列 `MS_perrep_mean`：按 rep 配对算 `mean_f |s_f,tp,r − s_f,CTRL,r|` 再对 rep 取均值（非作者定义，仅对照）。
- t 检验：Mem 份额 tp 组（4 个 rep）vs CTRL 组（4 个 rep），`scipy.stats.ttest_ind(equal_var=False)`（Welch），双侧；表中注明 n = 4 / 4，只作参考。

---

## 任务 1：份额表 → `01_shares.tsv` + `01_shares_secondary.tsv` + `01_groups.tsv` + `01_shares.md` + `01_shares.py`
- `01_groups.tsv`：`gene, protein_group, genes_raw, n_runs_detected, is_primary`（每个 gene×group 一行；gene 没测到则一行 `n_runs_detected = 0`）。
- `01_shares.tsv`（primary group）与 `01_shares_secondary.tsv`（其余 group），列固定：`gene, protein_group, timepoint, rep, normalization, n_fractions_present, Cyt, Mem, Nuc`。行 = 4 gene × 5 tp × 4 rep × 2 normalization = 160（数据齐全时）；某 (tp, rep) 一个 fraction 都没有则不出行，并在 md 记录。
- `01_shares.md`：① 输入（文件行数、解析成功 run 数/120、四个蛋白各在多少 run 测到 /120、group 表）② 定义 ③ 缺 fraction 清单 ④ 份额表（primary，两版）⑤ 附表说明 ⑥ 跳过项。
- 验收：Cyt+Mem+Nuc = 1 ± 1e-9（三格都非空的行）；raw 与 centered 行数相同；n_fractions_present ∈ 1..6；0 行输入时所有表只有表头、md 写明 0/120。

## 任务 2：三个易位事件 → `02_events.md` + `02_events.tsv` + `02_events.png` + `02_events.py`
（依赖 `01_shares.tsv`；EGFR 总量一节由任务 3 产出后并入本 md 末尾。）
- 对 GRB2、SHC1、CBL（primary group）、两版归一化：
  - 表 E-A：`gene, normalization, timepoint, compartment, n, mean, sd, delta_vs_ctrl`（3×2×5×3 = 90 行）。
  - 表 E-B：`gene, normalization, timepoint(2min|8min), n_mem_up_cyt_down, n_mem_up, n_cyt_down, denominator`（3×2×2 = 12 行；格式 `n/4`）。
  - 表 E-C：`gene, normalization, timepoint, MS_author, MS_max1, MS_max2, MS_perrep_mean, gt_0.1(yes/no)`（3×2×4 = 24 行，timepoint ∈ EGF 四个）。
  - 表 E-D：`gene, normalization, timepoint, mem_mean_tp, mem_mean_ctrl, t, p_welch, n_tp, n_ctrl`（24 行）。
- `02_events.tsv`：长格式 `table, gene, normalization, timepoint, key, value`，涵盖 E-A..E-D 全部数值。
- `02_events.png`：2 行（raw / centered）× 9 列（GRB2-Cyt, GRB2-Mem, GRB2-Nuc, SHC1-…, CBL-…），横轴 CTRL,2min,8min,20min,90min 等距，纵轴份额（0–1），每子图 4 条重复线；画图前调用 `dataviz` skill；dpi ≥ 150。0 行输入时仍出图框，子图内写"no data"。
- `02_events.md`：① 输入 ② 定义（含 §1 作者定义引用）③ E-A ④ E-B ⑤ E-C（含 ≥0.1 的时间点清单）⑥ E-D ⑦ EGFR 总量（抄 `03_egfr_total.md` 的结果节）⑧ 跳过项。
- 验收：E-A 的 mean 可由 01_shares.tsv 复算；E-B 的分母与 01_shares 中该 gene 两边都有值的 rep 数一致；E-C 的 MS 用 fraction 级均值算，不是用区室份额；p 值表注明 n=4。

## 任务 3：EGFR 总量 → `03_egfr_total.tsv` + `03_egfr_total.md` + `03_egfr_total.py`
（不依赖任务 1，可并行。）
- 对 EGFR（primary group；secondary 也算，放附表）：每个 (timepoint, rep)：`total = Σ_{present} q_f`（raw、centered 两版），`n_fractions_present`。
- 比值：`ratio_90_vs_ctrl_r = total_90min,r / total_CTRL,r`，按 rep 配对，4 个各报；两边任一缺则留空。
- `03_egfr_total.tsv`：`gene, protein_group, is_primary, timepoint, rep, normalization, n_fractions_present, total_quantity`；比值另存 `03_egfr_ratio.tsv`：`gene, protein_group, normalization, rep, total_90min, total_CTRL, ratio`。
- `03_egfr_total.md`：① 输入（EGFR 测到 run 数 /120）② 定义 ③ 总量表 ④ 90 min / CTRL 比值表（4 个 rep 各一行，两版）⑤ 跳过项。
- 验收：total 等于 01_shares 口径下同 (tp, rep) 的 6 个 quantity 之和（advisor 复算）；缺 fraction 的行标出。

## 验收总则
四个蛋白各在多少个 run 测到（/120）必须报；份额行之和 = 1（浮点误差内）；raw、centered 两版都齐；方向计数 `n/4`；Movement Score 定义原文与行号见 §1。0 行输入时上述项以"缺数据，跳过"形式出现，脚本须能在数据到位后直接复跑。
