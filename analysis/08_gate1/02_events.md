# 08_gate1 任务 2：GRB2 / SHC1 / CBL 三个易位事件

生成脚本：`analysis/08_gate1/02_events.py`（按 `analysis/08_gate1/00_design.md` §0、§1、§2、任务 2 执行；`python3`，任意目录可复跑，输出确定）。只报数和方向计数。

输出：`02_events.tsv`（长格式，E-A..E-D 全部数值）、`02_events.png`、本文件。

## ① 输入

- `data/pilot/gate1_proteins.tsv`：数据行 0（不含表头）；distinct run 0；解析成功 0/120；解析失败 0；多 run 同三元组 0；值冲突格 0。
- `analysis/07_pilot2/01_run_factors.tsv`：layer=proteome 120 行，可用 (timepoint, fraction, rep) → a 120 个。
- `analysis/08_gate1/01_groups.tsv`：0 行；primary（is_primary=yes）：无。
- `analysis/08_gate1/01_shares.tsv`：0 行；其中 GRB2 0、SHC1 0、CBL 0 行。
- `analysis/08_gate1/03_egfr_total.md`：⑦ 节来源。
- **四个蛋白各在多少 run 测到（/120，任一 protein_group 有有效 quantity 的三元组数，口径同任务 1）**：GRB2 0/120，SHC1 0/120，CBL 0/120，EGFR 0/120。

## ② 定义

### 作者 Movement Score（`00_design.md` §1，原文与行号）

`code/Protein contour/Phospho/SpatialProteoDynamics.github.io-v1.0/SpatialProteoDynamics-SpatialProteoDynamics.github.io-a6b8aac/DataProcessing/translocation_plots.R`（运行时逐行核对，见 ⑧ 自检）：

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
181   geom_text_repel(aes(label=Gene), data=temp_table[(temp_table$pval_combi_FDR<0.05) & (temp_table$MS_2min_mean>0.1),])+
```
（:181 行尾的 `+` 是 ggplot 链式续行符，`00_design.md` §1 的引用省略了它；其余 11 行与 §1 引用逐字相同。）

即：① log2 值还原为线性（:19-20）；② 每个 fraction 内 4 个重复取**均值**（:23-25），得 6 个 fraction 均值；③ 份额 = 该 fraction 均值 / 6 个均值之和（:30-31）；④ 每个 fraction 的 |份额_tp − 份额_CTRL|（:33）；⑤ **Movement Score = 这 6 个绝对差的均值**（:34）；⑥ 差最大和次大的 fraction 记为 MS_max1 / MS_max2（:36, :47）；⑦ 图上标注阈值 `MS_2min_mean > 0.1` 且 `pval_combi_FDR < 0.05`（:181）。作者只对 2 min vs CTRL 算。

### 共用定义（`00_design.md` §2，照抄任务 2 用到的条目）

- 份额：对每个 (gene, protein_group, timepoint, rep)：`n_fractions_present` = 6 个 fraction 中有值的个数；`share_f = q_f / Σ_{present} q_f`；`Cyt = share_FR1 + share_FR2`、`Mem = share_FR3 + share_FR4`、`Nuc = share_FR5 + share_FR6`；某区室两个 fraction 都缺 → 该区室格留空；不填补。（取自 `01_shares.tsv`，primary group。）
- 区室均值/标准差：按 timepoint 对 4 个重复的份额取 `mean`、`sd(ddof=1)`，n 不足 4 时如实写 n。
- Δ份额 = `mean_tp − mean_CTRL`。
- 方向计数（2min、8min 各一次）：按 rep 编号配对，`Mem_tp,r − Mem_CTRL,r > 0` 且 `Cyt_tp,r − Cyt_CTRL,r < 0` 的 rep 数，写成 `n/4`（分母 = 两边都有份额的 rep 数，若不足 4 写实际分母并注明）；另报"仅 Mem 上升"、"仅 Cyt 下降"的 n/4。
- Movement Score（作者定义，§1）：对每个 timepoint：fraction 均值 `m_f,tp = mean over reps of q_f`（线性值；缺失的 rep 不计入均值，与作者在填补后无缺失的前提不同）；`s_f,tp = m_f,tp / Σ_f m_f,tp`；`MS_tp = mean_f |s_f,tp − s_f,CTRL|`；同时报 `MS_max1`、`MS_max2`（fraction 编号）。阈值 0.1（:181）。raw、centered 各算。另加对照列 `MS_perrep_mean`：按 rep 配对算 `mean_f |s_f,tp,r − s_f,CTRL,r|` 再对 rep 取均值（非作者定义，仅对照）。
- t 检验：Mem 份额 tp 组（4 个 rep）vs CTRL 组（4 个 rep），`scipy.stats.ttest_ind(equal_var=False)`（Welch），双侧；n = 4 / 4，只作参考。
- 两版归一化：`raw` = quantity；`centered` = quantity × 2^a（a 取 `01_run_factors.tsv` 中 layer=proteome、同 (timepoint, fraction, rep) 的 `a`）。fraction 级数值的 run 解析、基因匹配、值冲突规则与任务 1 相同；primary group 取 `01_groups.tsv` 中 `is_primary = yes`。

### 脚本补充规则（方案未写明或有歧义处；未替换方案定义）

- **E-B 分母（歧义）**：方案写"两边都有份额的 rep 数"。本脚本取 tp 与 CTRL 两边 `Mem`、`Cyt` 四格都非空的 rep；三个计数都在这同一组 rep 上数，格式 `n/分母`；分母 = 0 时三格留空。差值 = 0 不计入"上升"或"下降"。
- **"仅 Mem 上升" / "仅 Cyt 下降"（歧义）**：按表 E-B 列名 `n_mem_up` / `n_cyt_down` 理解为单看一项（`Mem_tp,r − Mem_CTRL,r > 0`，不看 Cyt；`Cyt_tp,r − Cyt_CTRL,r < 0`，不看 Mem）。若"仅"指排他（Mem 升且 Cyt 不降），其值 = `n_mem_up − n_mem_up_cyt_down`（同一组 rep），`n_cyt_down` 同理。
- **E-C 缺整个 fraction（方案未定义）**：tp 或 CTRL 中任一 fraction 在 4 个 rep 全缺，或 6 个均值之和 = 0 → `MS_author`、`MS_max1`、`MS_max2`、`ge_0.1` 留空，在 ⑤ 列出；不在少于 6 个 fraction 上算。
- **MS_max2**：照作者 `maxN`（:38-45）取第二大的值，再 `match` 第一个等于它的位置；若最大值并列，`MS_max2` 与 `MS_max1` 相同（作者代码行为）。写成 `FRk`，对应作者的整数 k。
- **MS_perrep_mean**：某 rep 只有在 tp 与 CTRL 两边 6 个 fraction 都有值时才计入；计入 rep 数不足 4 的在 ⑤ 列出。
- **阈值（歧义）**：方案列名 `ge_0.1` 与"≥0.1 的时间点清单"用 ≥ 0.1；作者 :181 为严格 `> 0.1`。本脚本按方案用 ≥；两者只在 MS 恰为 0.1 时不同。作者的 `pval_combi_FDR < 0.05` 条件（limma + sumlog + BH）本轮不复现。
- quantity 按线性值使用（§2"线性值"），不做 2^x（作者 :19-20 的 2^ 针对其 log2 输入）。
- E-A：CTRL 行的 `delta_vs_ctrl` = 0；n < 2 时 sd 留空；n = 0 时 mean、delta 留空。
- E-D：tp 组在前、CTRL 组在后（t > 0 表示 tp 的 Mem 均值更高）；任一组 n < 2 或结果为 nan 时 t、p 留空；只用非空 Mem 份额，n_tp、n_ctrl 为实际个数。
- 02_events.tsv：`key` 在 E-A 中为 `<区室>_<统计量>`（如 `Mem_mean`、`Mem_n`），其余表为列名；数值 12 位有效数字，空值写空串。md 表中份额、均值、sd、MS 保留 4 位小数，t 保留 3 位，p 保留 3 位有效数字。
- 图：份额取自 `01_shares.tsv`；某 rep 某 tp 份额缺 → 该点断开；配色为 dataviz 参考调色板分类色 1–4（固定顺序），另以 marker 形状区分 rep。

## ③ 表 E-A：区室份额均值 / 标准差 / Δ份额（primary group，两版）

缺数据：`gate1_proteins.tsv` 0 行，跳过。

## ④ 表 E-B：方向计数（Mem 上升且 Cyt 下降，按 rep 配对 vs CTRL）

缺数据：`gate1_proteins.tsv` 0 行，跳过。

## ⑤ 表 E-C：Movement Score（作者定义，fraction 级重复均值；对照列 MS_perrep_mean）

缺数据：`gate1_proteins.tsv` 0 行，跳过。

## ⑥ 表 E-D：Mem 份额 tp vs CTRL，Welch t 检验（n = 4 / 4，只作参考）

缺数据：`gate1_proteins.tsv` 0 行，跳过。

## ⑦ EGFR 总量（抄自 `03_egfr_total.md` 的 ③④ 两节）

以下原样抄自 `analysis/08_gate1/03_egfr_total.md`（只把标题降一级，正文未改）。

### ③ EGFR 总量表

EGFR 在 0/120 个 run 测到。

缺数据：`gate1_proteins.tsv` 0 行，跳过

### ④ EGFR 90 min / CTRL 比值表

EGFR 在 0/120 个 run 测到。

缺数据：`gate1_proteins.tsv` 0 行，跳过

## ⑧ 跳过项

- 缺数据：`gate1_proteins.tsv` 0 行，跳过：③ E-A、④ E-B、⑤ E-C、⑥ E-D 均无数值；`02_events.tsv` 只有表头；`02_events.png` 只有 2 × 9 图框，子图内写 "no data"。数据到位、任务 1 复跑后直接复跑本脚本即可。
- 四个蛋白各在 run 测到：GRB2 0/120，SHC1 0/120，CBL 0/120，EGFR 0/120。
- 不复现 DAPAR 的 LOESS 归一化与 KNN/MEC/detQuant 填补（`00_design.md` §1），只用 raw 与 centered 两版；不填补缺失。
- 不复现作者 :58-153 的 limma p 值合并与 :173 的 BH 校正（无 `pval_combi_FDR`）；E-C 只按 MS 阈值列出。
- 脚本在 scratchpad 合成夹具上跑通。

### 自检（方案验收项）

| 检查 | 结果 | 值 |
|---|---|---|
| 输入 01_shares.tsv / 01_groups.tsv 列 = 任务 1 方案列 | PASS |  |
| 01_groups.tsv：每个基因至多一个 is_primary=yes | PASS |  |
| 01_shares.tsv 的 protein_group = 01_groups.tsv 的 primary | 不适用 | 0 行 |
| 01_shares.tsv 与 fraction 级重算份额一致（GRB2/SHC1/CBL） | 不适用 | 0 行 |
| ⑦ 节 = 03_egfr_total.md 的 ③④ 原文（仅标题降一级） | PASS | analysis/08_gate1/03_egfr_total.md |
| §1 引用的 R 原文与 translocation_plots.R 行号逐行一致 | PASS | 12 行 |
| 02_events.tsv 列 = table, gene, normalization, timepoint, key, value | PASS |  |
| 0 行输入：E-A..E-D 全为 0 行，02_events.tsv 只有表头 | PASS | E-A 0 / E-B 0 / E-C 0 / E-D 0；tsv 0 行 |
| 02_events.tsv 覆盖 E-A..E-D 全部数值（行数 = 4·EA + 4·EB + 5·EC + 6·ED，逐格回读一致） | PASS | 0 / 0，不一致 0 |
| E-A：n / mean / sd(ddof=1) 由 01_shares.tsv 复算一致 | 不适用 | 0 行 |
| E-B：分母 = 01_shares 中两边 Mem、Cyt 都有值的 rep 数 | 不适用 | 0 行 |
| E-C：MS_author 由 fraction 级（6 个 fraction）重复均值复算一致 | 不适用 | 0 行 |
| E-D：t / p_welch 与 Welch 公式手算一致（双侧） | 不适用 | 0 行 |
| 02_events.png：子图数 = 18（2 行 × 9 列），dpi ≥ 150 | PASS | 子图 18，dpi 200，3700×1120 px |
| 02_events.png：无数据时 18 个子图均写 'no data' | PASS | 'no data' 子图 18 |

