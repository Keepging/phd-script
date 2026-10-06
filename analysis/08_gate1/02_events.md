# 08_gate1 任务 2：GRB2 / SHC1 / CBL 三个易位事件

生成脚本：`analysis/08_gate1/02_events.py`（按 `analysis/08_gate1/00_design.md` §0、§1、§2、任务 2 执行；`python3`，任意目录可复跑，输出确定）。只报数和方向计数。

输出：`02_events.tsv`（长格式，E-A..E-D 全部数值）、`02_events.png`、本文件。

## ① 输入

- `analysis/gate1_proteins.tsv`：数据行 344（不含表头）；distinct run 119；解析成功 119/120；解析失败 0；多 run 同三元组 0；值冲突格 0。
- `analysis/07_pilot2/01_run_factors.tsv`：layer=proteome 120 行，可用 (timepoint, fraction, rep) → a 120 个。
- `analysis/08_gate1/01_groups.tsv`：4 行；primary（is_primary=yes）：GRB2 `P62993`、SHC1 `P29353`、CBL `P22681`、EGFR `P00533`。
- `analysis/08_gate1/01_shares.tsv`：160 行；其中 GRB2 40、SHC1 40、CBL 40 行。
- `analysis/08_gate1/03_egfr_total.md`：⑦ 节来源。
- **四个蛋白各在多少 run 测到（/120，任一 protein_group 有有效 quantity 的三元组数，口径同任务 1）**：GRB2 90/120，SHC1 115/120，CBL 59/120，EGFR 80/120。

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
- **E-C 缺整个 fraction（方案未定义）**：tp 或 CTRL 中任一 fraction 在 4 个 rep 全缺，或 6 个均值之和 = 0 → `MS_author`、`MS_max1`、`MS_max2`、`gt_0.1` 留空，在 ⑤ 列出；不在少于 6 个 fraction 上算。
- **MS_max2**：照作者 `maxN`（:38-45）取第二大的值，再 `match` 第一个等于它的位置；若最大值并列，`MS_max2` 与 `MS_max1` 相同（作者代码行为）。写成 `FRk`，对应作者的整数 k。
- **MS_perrep_mean**：某 rep 只有在 tp 与 CTRL 两边 6 个 fraction 都有值时才计入；计入 rep 数不足 4 的在 ⑤ 列出。
- **阈值**：按作者 :181 的严格 `> 0.1`（2026-10-06 确认，原方案的 ≥ 0.1 已改）。作者的 `pval_combi_FDR < 0.05` 条件（limma + sumlog + BH）本轮不复现。
- quantity 按线性值使用（§2"线性值"），不做 2^x（作者 :19-20 的 2^ 针对其 log2 输入）。
- E-A：CTRL 行的 `delta_vs_ctrl` = 0；n < 2 时 sd 留空；n = 0 时 mean、delta 留空。
- E-D：tp 组在前、CTRL 组在后（t > 0 表示 tp 的 Mem 均值更高）；任一组 n < 2 或结果为 nan 时 t、p 留空；只用非空 Mem 份额，n_tp、n_ctrl 为实际个数。
- 02_events.tsv：`key` 在 E-A 中为 `<区室>_<统计量>`（如 `Mem_mean`、`Mem_n`），其余表为列名；数值 12 位有效数字，空值写空串。md 表中份额、均值、sd、MS 保留 4 位小数，t 保留 3 位，p 保留 3 位有效数字。
- 图：份额取自 `01_shares.tsv`；某 rep 某 tp 份额缺 → 该点断开；配色为 dataviz 参考调色板分类色 1–4（固定顺序），另以 marker 形状区分 rep。

## ③ 表 E-A：区室份额均值 / 标准差 / Δ份额（primary group，两版）

| gene | normalization | timepoint | compartment | n | mean | sd | delta_vs_ctrl |
|---|---|---|---|---|---|---|---|
| GRB2 | raw | CTRL | Cyt | 4 | 0.7758 | 0.1874 | 0.0000 |
| GRB2 | raw | CTRL | Mem | 4 | 0.1872 | 0.1165 | 0.0000 |
| GRB2 | raw | CTRL | Nuc | 1 | 0.1480 |  | 0.0000 |
| GRB2 | raw | 2min | Cyt | 4 | 0.3454 | 0.0598 | -0.4304 |
| GRB2 | raw | 2min | Mem | 4 | 0.6368 | 0.0765 | 0.4496 |
| GRB2 | raw | 2min | Nuc | 2 | 0.0356 | 0.0093 | -0.1124 |
| GRB2 | raw | 8min | Cyt | 4 | 0.3236 | 0.0644 | -0.4521 |
| GRB2 | raw | 8min | Mem | 4 | 0.6208 | 0.0847 | 0.4336 |
| GRB2 | raw | 8min | Nuc | 3 | 0.0741 | 0.0308 | -0.0739 |
| GRB2 | raw | 20min | Cyt | 4 | 0.4416 | 0.0749 | -0.3341 |
| GRB2 | raw | 20min | Mem | 4 | 0.5085 | 0.0703 | 0.3213 |
| GRB2 | raw | 20min | Nuc | 2 | 0.0997 | 0.0510 | -0.0483 |
| GRB2 | raw | 90min | Cyt | 4 | 0.6716 | 0.0386 | -0.1042 |
| GRB2 | raw | 90min | Mem | 4 | 0.3181 | 0.0350 | 0.1309 |
| GRB2 | raw | 90min | Nuc | 1 | 0.0411 |  | -0.1069 |
| GRB2 | centered | CTRL | Cyt | 4 | 0.7208 | 0.1615 | 0.0000 |
| GRB2 | centered | CTRL | Mem | 4 | 0.2416 | 0.0924 | 0.0000 |
| GRB2 | centered | CTRL | Nuc | 1 | 0.1504 |  | 0.0000 |
| GRB2 | centered | 2min | Cyt | 4 | 0.2852 | 0.0461 | -0.4356 |
| GRB2 | centered | 2min | Mem | 4 | 0.6852 | 0.0790 | 0.4436 |
| GRB2 | centered | 2min | Nuc | 2 | 0.0593 | 0.0088 | -0.0911 |
| GRB2 | centered | 8min | Cyt | 4 | 0.2575 | 0.0336 | -0.4633 |
| GRB2 | centered | 8min | Mem | 4 | 0.6572 | 0.0691 | 0.4156 |
| GRB2 | centered | 8min | Nuc | 3 | 0.1137 | 0.0383 | -0.0367 |
| GRB2 | centered | 20min | Cyt | 4 | 0.3437 | 0.0702 | -0.3771 |
| GRB2 | centered | 20min | Mem | 4 | 0.5884 | 0.0583 | 0.3468 |
| GRB2 | centered | 20min | Nuc | 2 | 0.1357 | 0.0615 | -0.0147 |
| GRB2 | centered | 90min | Cyt | 4 | 0.6014 | 0.0430 | -0.1194 |
| GRB2 | centered | 90min | Mem | 4 | 0.3841 | 0.0394 | 0.1425 |
| GRB2 | centered | 90min | Nuc | 1 | 0.0580 |  | -0.0924 |
| SHC1 | raw | CTRL | Cyt | 4 | 0.7330 | 0.1308 | 0.0000 |
| SHC1 | raw | CTRL | Mem | 4 | 0.2177 | 0.0814 | 0.0000 |
| SHC1 | raw | CTRL | Nuc | 4 | 0.0493 | 0.0538 | 0.0000 |
| SHC1 | raw | 2min | Cyt | 4 | 0.4693 | 0.0411 | -0.2636 |
| SHC1 | raw | 2min | Mem | 4 | 0.4655 | 0.0398 | 0.2478 |
| SHC1 | raw | 2min | Nuc | 4 | 0.0651 | 0.0130 | 0.0158 |
| SHC1 | raw | 8min | Cyt | 4 | 0.4792 | 0.0331 | -0.2537 |
| SHC1 | raw | 8min | Mem | 4 | 0.4316 | 0.0540 | 0.2139 |
| SHC1 | raw | 8min | Nuc | 4 | 0.0891 | 0.0230 | 0.0398 |
| SHC1 | raw | 20min | Cyt | 4 | 0.6003 | 0.0557 | -0.1326 |
| SHC1 | raw | 20min | Mem | 4 | 0.2989 | 0.0606 | 0.0812 |
| SHC1 | raw | 20min | Nuc | 4 | 0.1008 | 0.0614 | 0.0515 |
| SHC1 | raw | 90min | Cyt | 4 | 0.6681 | 0.0253 | -0.0648 |
| SHC1 | raw | 90min | Mem | 4 | 0.2745 | 0.0289 | 0.0568 |
| SHC1 | raw | 90min | Nuc | 4 | 0.0573 | 0.0234 | 0.0080 |
| SHC1 | centered | CTRL | Cyt | 4 | 0.6718 | 0.1190 | 0.0000 |
| SHC1 | centered | CTRL | Mem | 4 | 0.2634 | 0.0836 | 0.0000 |
| SHC1 | centered | CTRL | Nuc | 4 | 0.0648 | 0.0500 | 0.0000 |
| SHC1 | centered | 2min | Cyt | 4 | 0.3974 | 0.0252 | -0.2744 |
| SHC1 | centered | 2min | Mem | 4 | 0.5023 | 0.0343 | 0.2389 |
| SHC1 | centered | 2min | Nuc | 4 | 0.1003 | 0.0187 | 0.0355 |
| SHC1 | centered | 8min | Cyt | 4 | 0.3949 | 0.0233 | -0.2769 |
| SHC1 | centered | 8min | Mem | 4 | 0.4747 | 0.0416 | 0.2114 |
| SHC1 | centered | 8min | Nuc | 4 | 0.1303 | 0.0389 | 0.0655 |
| SHC1 | centered | 20min | Cyt | 4 | 0.4943 | 0.0481 | -0.1775 |
| SHC1 | centered | 20min | Mem | 4 | 0.3471 | 0.0639 | 0.0837 |
| SHC1 | centered | 20min | Nuc | 4 | 0.1586 | 0.0900 | 0.0938 |
| SHC1 | centered | 90min | Cyt | 4 | 0.5918 | 0.0416 | -0.0800 |
| SHC1 | centered | 90min | Mem | 4 | 0.3240 | 0.0342 | 0.0606 |
| SHC1 | centered | 90min | Nuc | 4 | 0.0841 | 0.0309 | 0.0193 |
| CBL | raw | CTRL | Cyt | 4 | 0.6200 | 0.4126 | 0.0000 |
| CBL | raw | CTRL | Mem | 4 | 0.3800 | 0.4126 | 0.0000 |
| CBL | raw | CTRL | Nuc | 0 |  |  |  |
| CBL | raw | 2min | Cyt | 4 | 0.1223 | 0.0231 | -0.4977 |
| CBL | raw | 2min | Mem | 4 | 0.6983 | 0.1189 | 0.3183 |
| CBL | raw | 2min | Nuc | 3 | 0.2392 | 0.0135 |  |
| CBL | raw | 8min | Cyt | 3 | 0.1541 | 0.0820 | -0.4659 |
| CBL | raw | 8min | Mem | 4 | 0.7961 | 0.1105 | 0.4161 |
| CBL | raw | 8min | Nuc | 3 | 0.1178 | 0.0577 |  |
| CBL | raw | 20min | Cyt | 4 | 0.1297 | 0.1036 | -0.4903 |
| CBL | raw | 20min | Mem | 4 | 0.7993 | 0.0858 | 0.4193 |
| CBL | raw | 20min | Nuc | 2 | 0.1420 | 0.0977 |  |
| CBL | raw | 90min | Cyt | 4 | 0.7483 | 0.4614 | 0.1283 |
| CBL | raw | 90min | Mem | 2 | 0.5033 | 0.6208 | 0.1233 |
| CBL | raw | 90min | Nuc | 0 |  |  |  |
| CBL | centered | CTRL | Cyt | 4 | 0.5441 | 0.3666 | 0.0000 |
| CBL | centered | CTRL | Mem | 4 | 0.4559 | 0.3666 | 0.0000 |
| CBL | centered | CTRL | Nuc | 0 |  |  |  |
| CBL | centered | 2min | Cyt | 4 | 0.0903 | 0.0087 | -0.4538 |
| CBL | centered | 2min | Mem | 4 | 0.6677 | 0.1721 | 0.2118 |
| CBL | centered | 2min | Nuc | 3 | 0.3227 | 0.0789 |  |
| CBL | centered | 8min | Cyt | 3 | 0.1049 | 0.0716 | -0.4392 |
| CBL | centered | 8min | Mem | 4 | 0.7931 | 0.1339 | 0.3372 |
| CBL | centered | 8min | Nuc | 3 | 0.1710 | 0.0868 |  |
| CBL | centered | 20min | Cyt | 4 | 0.0923 | 0.0855 | -0.4518 |
| CBL | centered | 20min | Mem | 4 | 0.8379 | 0.0697 | 0.3820 |
| CBL | centered | 20min | Nuc | 2 | 0.1396 | 0.0938 |  |
| CBL | centered | 90min | Cyt | 4 | 0.7297 | 0.4642 | 0.1855 |
| CBL | centered | 90min | Mem | 2 | 0.5407 | 0.5950 | 0.0848 |
| CBL | centered | 90min | Nuc | 0 |  |  |  |

n < 4 的格：24：GRB2/raw/CTRL/Nuc n=1；GRB2/raw/2min/Nuc n=2；GRB2/raw/8min/Nuc n=3；GRB2/raw/20min/Nuc n=2；GRB2/raw/90min/Nuc n=1；GRB2/centered/CTRL/Nuc n=1；GRB2/centered/2min/Nuc n=2；GRB2/centered/8min/Nuc n=3；GRB2/centered/20min/Nuc n=2；GRB2/centered/90min/Nuc n=1；CBL/raw/CTRL/Nuc n=0；CBL/raw/2min/Nuc n=3；CBL/raw/8min/Cyt n=3；CBL/raw/8min/Nuc n=3；CBL/raw/20min/Nuc n=2；CBL/raw/90min/Mem n=2；CBL/raw/90min/Nuc n=0；CBL/centered/CTRL/Nuc n=0；CBL/centered/2min/Nuc n=3；CBL/centered/8min/Cyt n=3；CBL/centered/8min/Nuc n=3；CBL/centered/20min/Nuc n=2；CBL/centered/90min/Mem n=2；CBL/centered/90min/Nuc n=0。

## ④ 表 E-B：方向计数（Mem 上升且 Cyt 下降，按 rep 配对 vs CTRL）

| gene | normalization | timepoint | n_mem_up_cyt_down | n_mem_up | n_cyt_down | denominator |
|---|---|---|---|---|---|---|
| GRB2 | raw | 2min | 4/4 | 4/4 | 4/4 | 4 |
| GRB2 | raw | 8min | 4/4 | 4/4 | 4/4 | 4 |
| GRB2 | centered | 2min | 4/4 | 4/4 | 4/4 | 4 |
| GRB2 | centered | 8min | 4/4 | 4/4 | 4/4 | 4 |
| SHC1 | raw | 2min | 4/4 | 4/4 | 4/4 | 4 |
| SHC1 | raw | 8min | 4/4 | 4/4 | 4/4 | 4 |
| SHC1 | centered | 2min | 4/4 | 4/4 | 4/4 | 4 |
| SHC1 | centered | 8min | 4/4 | 4/4 | 4/4 | 4 |
| CBL | raw | 2min | 3/4 | 3/4 | 3/4 | 4 |
| CBL | raw | 8min | 2/3 | 2/3 | 2/3 | 3 |
| CBL | centered | 2min | 3/4 | 3/4 | 3/4 | 4 |
| CBL | centered | 8min | 2/3 | 2/3 | 2/3 | 3 |

分母不足 4 的行：2：CBL/raw/8min 分母 3；CBL/centered/8min 分母 3。

## ⑤ 表 E-C：Movement Score（作者定义，fraction 级重复均值；对照列 MS_perrep_mean）

| gene | normalization | timepoint | MS_author | MS_max1 | MS_max2 | MS_perrep_mean | gt_0.1 |
|---|---|---|---|---|---|---|---|
| GRB2 | raw | 2min |  |  |  |  |  |
| GRB2 | raw | 8min | 0.1298 | FR4 | FR1 |  | yes |
| GRB2 | raw | 20min |  |  |  |  |  |
| GRB2 | raw | 90min |  |  |  |  |  |
| GRB2 | centered | 2min |  |  |  |  |  |
| GRB2 | centered | 8min | 0.1288 | FR4 | FR1 |  | yes |
| GRB2 | centered | 20min |  |  |  |  |  |
| GRB2 | centered | 90min |  |  |  |  |  |
| SHC1 | raw | 2min | 0.0981 | FR4 | FR1 | 0.0910 | no |
| SHC1 | raw | 8min | 0.0883 | FR4 | FR1 | 0.0524 | no |
| SHC1 | raw | 20min | 0.0543 | FR4 | FR5 | 0.0469 | no |
| SHC1 | raw | 90min | 0.0288 | FR5 | FR4 |  | no |
| SHC1 | centered | 2min | 0.1013 | FR4 | FR1 | 0.0750 | yes |
| SHC1 | centered | 8min | 0.0951 | FR4 | FR1 | 0.0659 | no |
| SHC1 | centered | 20min | 0.0679 | FR4 | FR6 | 0.0389 | no |
| SHC1 | centered | 90min | 0.0371 | FR5 | FR4 |  | no |
| CBL | raw | 2min |  |  |  |  |  |
| CBL | raw | 8min |  |  |  |  |  |
| CBL | raw | 20min |  |  |  |  |  |
| CBL | raw | 90min |  |  |  |  |  |
| CBL | centered | 2min |  |  |  |  |  |
| CBL | centered | 8min |  |  |  |  |  |
| CBL | centered | 20min |  |  |  |  |  |
| CBL | centered | 90min |  |  |  |  |  |

**MS_author > 0.1 的 (gene, normalization, timepoint)**：3 个：(GRB2, raw, 8min)；(GRB2, centered, 8min)；(SHC1, centered, 2min)。

- MS_author 留空（tp 或 CTRL 缺整个 fraction）：14：GRB2/raw/2min（缺 2min/FR5）；GRB2/raw/20min（缺 20min/FR5）；GRB2/raw/90min（缺 90min/FR5）；GRB2/centered/2min（缺 2min/FR5）；GRB2/centered/20min（缺 20min/FR5）；GRB2/centered/90min（缺 90min/FR5）；CBL/raw/2min（缺 CTRL/FR5, CTRL/FR6）；CBL/raw/8min（缺 CTRL/FR5, CTRL/FR6）；CBL/raw/20min（缺 CTRL/FR5, CTRL/FR6, 20min/FR2, 20min/FR5）；CBL/raw/90min（缺 CTRL/FR5, CTRL/FR6, 90min/FR4, 90min/FR5, 90min/FR6）；CBL/centered/2min（缺 CTRL/FR5, CTRL/FR6）；CBL/centered/8min（缺 CTRL/FR5, CTRL/FR6）；CBL/centered/20min（缺 CTRL/FR5, CTRL/FR6, 20min/FR2, 20min/FR5）；CBL/centered/90min（缺 CTRL/FR5, CTRL/FR6, 90min/FR4, 90min/FR5, 90min/FR6）。
- fraction 均值所用 rep 数不足 4（缺失 rep 不计入均值）：52：GRB2/raw/CTRL/FR4 n=3；GRB2/raw/CTRL/FR5 n=1；GRB2/raw/CTRL/FR6 n=1；GRB2/raw/2min/FR6 n=2；GRB2/raw/8min/FR5 n=1；GRB2/raw/8min/FR6 n=3；GRB2/raw/20min/FR6 n=2；GRB2/raw/90min/FR6 n=1；GRB2/centered/CTRL/FR4 n=3；GRB2/centered/CTRL/FR5 n=1；GRB2/centered/CTRL/FR6 n=1；GRB2/centered/2min/FR6 n=2；GRB2/centered/8min/FR5 n=1；GRB2/centered/8min/FR6 n=3；GRB2/centered/20min/FR6 n=2；GRB2/centered/90min/FR6 n=1；SHC1/raw/CTRL/FR5 n=1；SHC1/raw/90min/FR5 n=3；SHC1/raw/90min/FR6 n=3；SHC1/centered/CTRL/FR5 n=1；SHC1/centered/90min/FR5 n=3；SHC1/centered/90min/FR6 n=3；CBL/raw/CTRL/FR4 n=1；CBL/raw/2min/FR2 n=2；CBL/raw/2min/FR3 n=1；CBL/raw/2min/FR5 n=2；CBL/raw/2min/FR6 n=1；CBL/raw/8min/FR1 n=2；CBL/raw/8min/FR2 n=2；CBL/raw/8min/FR3 n=2；CBL/raw/8min/FR5 n=1；CBL/raw/8min/FR6 n=3；CBL/raw/20min/FR3 n=2；CBL/raw/20min/FR4 n=3；CBL/raw/20min/FR6 n=2；CBL/raw/90min/FR2 n=1；CBL/raw/90min/FR3 n=2；CBL/centered/CTRL/FR4 n=1；CBL/centered/2min/FR2 n=2；CBL/centered/2min/FR3 n=1；CBL/centered/2min/FR5 n=2；CBL/centered/2min/FR6 n=1；CBL/centered/8min/FR1 n=2；CBL/centered/8min/FR2 n=2；CBL/centered/8min/FR3 n=2；CBL/centered/8min/FR5 n=1；CBL/centered/8min/FR6 n=3；CBL/centered/20min/FR3 n=2；CBL/centered/20min/FR4 n=3；CBL/centered/20min/FR6 n=2；CBL/centered/90min/FR2 n=1；CBL/centered/90min/FR3 n=2。
- MS_perrep_mean 计入 rep 数不足 4：24：GRB2/raw/2min n=0；GRB2/raw/8min n=0；GRB2/raw/20min n=0；GRB2/raw/90min n=0；GRB2/centered/2min n=0；GRB2/centered/8min n=0；GRB2/centered/20min n=0；GRB2/centered/90min n=0；SHC1/raw/2min n=1；SHC1/raw/8min n=1；SHC1/raw/20min n=1；SHC1/raw/90min n=0；SHC1/centered/2min n=1；SHC1/centered/8min n=1；SHC1/centered/20min n=1；SHC1/centered/90min n=0；CBL/raw/2min n=0；CBL/raw/8min n=0；CBL/raw/20min n=0；CBL/raw/90min n=0；CBL/centered/2min n=0；CBL/centered/8min n=0；CBL/centered/20min n=0；CBL/centered/90min n=0。

## ⑥ 表 E-D：Mem 份额 tp vs CTRL，Welch t 检验（n = 4 / 4，只作参考）

| gene | normalization | timepoint | mem_mean_tp | mem_mean_ctrl | t | p_welch | n_tp | n_ctrl |
|---|---|---|---|---|---|---|---|---|
| GRB2 | raw | 2min | 0.6368 | 0.1872 | 6.453 | 0.00116 | 4 | 4 |
| GRB2 | raw | 8min | 0.6208 | 0.1872 | 6.022 | 0.00131 | 4 | 4 |
| GRB2 | raw | 20min | 0.5085 | 0.1872 | 4.723 | 0.00542 | 4 | 4 |
| GRB2 | raw | 90min | 0.3181 | 0.1872 | 2.153 | 0.107 | 4 | 4 |
| GRB2 | centered | 2min | 0.6852 | 0.2416 | 7.299 | 0.000377 | 4 | 4 |
| GRB2 | centered | 8min | 0.6572 | 0.2416 | 7.205 | 0.000511 | 4 | 4 |
| GRB2 | centered | 20min | 0.5884 | 0.2416 | 6.349 | 0.00137 | 4 | 4 |
| GRB2 | centered | 90min | 0.3841 | 0.2416 | 2.837 | 0.0463 | 4 | 4 |
| SHC1 | raw | 2min | 0.4655 | 0.2177 | 5.469 | 0.00424 | 4 | 4 |
| SHC1 | raw | 8min | 0.4316 | 0.2177 | 4.378 | 0.00651 | 4 | 4 |
| SHC1 | raw | 20min | 0.2989 | 0.2177 | 1.599 | 0.165 | 4 | 4 |
| SHC1 | raw | 90min | 0.2745 | 0.2177 | 1.315 | 0.263 | 4 | 4 |
| SHC1 | centered | 2min | 0.5023 | 0.2634 | 5.288 | 0.00622 | 4 | 4 |
| SHC1 | centered | 8min | 0.4747 | 0.2634 | 4.526 | 0.00848 | 4 | 4 |
| SHC1 | centered | 20min | 0.3471 | 0.2634 | 1.590 | 0.166 | 4 | 4 |
| SHC1 | centered | 90min | 0.3240 | 0.2634 | 1.343 | 0.251 | 4 | 4 |
| CBL | raw | 2min | 0.6983 | 0.3800 | 1.482 | 0.222 | 4 | 4 |
| CBL | raw | 8min | 0.7961 | 0.3800 | 1.948 | 0.135 | 4 | 4 |
| CBL | raw | 20min | 0.7993 | 0.3800 | 1.990 | 0.133 | 4 | 4 |
| CBL | raw | 90min | 0.5033 | 0.3800 | 0.254 | 0.83 | 2 | 4 |
| CBL | centered | 2min | 0.6677 | 0.4559 | 1.046 | 0.351 | 4 | 4 |
| CBL | centered | 8min | 0.7931 | 0.4559 | 1.728 | 0.163 | 4 | 4 |
| CBL | centered | 20min | 0.8379 | 0.4559 | 2.047 | 0.127 | 4 | 4 |
| CBL | centered | 90min | 0.5407 | 0.4559 | 0.185 | 0.877 | 2 | 4 |

`scipy.stats.ttest_ind(equal_var=False)`，双侧，未做多重检验校正；每组 n ≤ 4，p 值只作参考。

## ⑦ EGFR 总量（抄自 `03_egfr_total.md` 的 ③④ 两节）

以下原样抄自 `analysis/08_gate1/03_egfr_total.md`（只把标题降一级，正文未改）。

### ③ EGFR 总量表

EGFR 在 80/120 个 run 测到。

gene = EGFR；primary protein_group = `P00533`（genes = `EGFR`，测到 80/120 run）。total_quantity = 该 (timepoint, rep) 下有值 fraction 的 quantity 之和，不填补；raw = quantity，centered = quantity × 2^a（a 来自 07 `01_run_factors.tsv`, proteome）。`有缺 fraction` = 是 表示 n_fractions_present < 6。数值为 6 位有效数字，完整精度见 `03_egfr_total.tsv`。

**raw**

| timepoint | rep | n_fractions_present | total_quantity | 缺的 fraction | 有缺 fraction |
|---|---|---|---|---|---|
| CTRL | Rep1 | 5 | 1.60857e+08 | FR1 | 是 |
| CTRL | Rep2 | 5 | 2.69796e+08 | FR1 | 是 |
| CTRL | Rep3 | 4 | 1.70902e+08 | FR1,FR2 | 是 |
| CTRL | Rep4 | 3 | 5.84781e+07 | FR1,FR2,FR3 | 是 |
| 2min | Rep1 | 3 | 1.60789e+08 | FR1,FR2,FR5 | 是 |
| 2min | Rep2 | 4 | 3.02834e+08 | FR1,FR2 | 是 |
| 2min | Rep3 | 4 | 2.17546e+08 | FR1,FR2 | 是 |
| 2min | Rep4 | 4 | 2.3859e+08 | FR1,FR2 | 是 |
| 8min | Rep1 | 4 | 1.87922e+08 | FR1,FR2 | 是 |
| 8min | Rep2 | 4 | 1.56687e+08 | FR1,FR2 | 是 |
| 8min | Rep3 | 4 | 1.59447e+08 | FR1,FR2 | 是 |
| 8min | Rep4 | 4 | 2.05479e+08 | FR1,FR2 | 是 |
| 20min | Rep1 | 4 | 1.25251e+08 | FR1,FR2 | 是 |
| 20min | Rep2 | 5 | 2.00439e+08 | FR1 | 是 |
| 20min | Rep3 | 4 | 1.77175e+08 | FR1,FR5 | 是 |
| 20min | Rep4 | 4 | 1.17991e+08 | FR1,FR2 | 是 |
| 90min | Rep1 | 3 | 7.57694e+07 | FR1,FR2,FR5 | 是 |
| 90min | Rep2 | 4 | 6.67118e+07 | FR1,FR2 | 是 |
| 90min | Rep3 | 4 | 6.79139e+07 | FR1,FR2 | 是 |
| 90min | Rep4 | 4 | 8.58314e+07 | FR1,FR2 | 是 |

**centered**

| timepoint | rep | n_fractions_present | total_quantity | 缺的 fraction | 有缺 fraction |
|---|---|---|---|---|---|
| CTRL | Rep1 | 5 | 1.82886e+08 | FR1 | 是 |
| CTRL | Rep2 | 5 | 2.68884e+08 | FR1 | 是 |
| CTRL | Rep3 | 4 | 2.1167e+08 | FR1,FR2 | 是 |
| CTRL | Rep4 | 3 | 5.06689e+07 | FR1,FR2,FR3 | 是 |
| 2min | Rep1 | 3 | 1.88419e+08 | FR1,FR2,FR5 | 是 |
| 2min | Rep2 | 4 | 3.08706e+08 | FR1,FR2 | 是 |
| 2min | Rep3 | 4 | 2.76379e+08 | FR1,FR2 | 是 |
| 2min | Rep4 | 4 | 2.74894e+08 | FR1,FR2 | 是 |
| 8min | Rep1 | 4 | 2.16183e+08 | FR1,FR2 | 是 |
| 8min | Rep2 | 4 | 1.85197e+08 | FR1,FR2 | 是 |
| 8min | Rep3 | 4 | 1.74191e+08 | FR1,FR2 | 是 |
| 8min | Rep4 | 4 | 2.75765e+08 | FR1,FR2 | 是 |
| 20min | Rep1 | 4 | 1.33487e+08 | FR1,FR2 | 是 |
| 20min | Rep2 | 5 | 2.4581e+08 | FR1 | 是 |
| 20min | Rep3 | 4 | 2.47208e+08 | FR1,FR5 | 是 |
| 20min | Rep4 | 4 | 1.35271e+08 | FR1,FR2 | 是 |
| 90min | Rep1 | 3 | 7.47593e+07 | FR1,FR2,FR5 | 是 |
| 90min | Rep2 | 4 | 7.30967e+07 | FR1,FR2 | 是 |
| 90min | Rep3 | 4 | 7.61634e+07 | FR1,FR2 | 是 |
| 90min | Rep4 | 4 | 8.58346e+07 | FR1,FR2 | 是 |

缺 fraction 清单（primary，n_fractions_present < 6）：40 行。

| normalization | timepoint | rep | 缺的 fraction |
|---|---|---|---|
| raw | CTRL | Rep1 | FR1 |
| raw | CTRL | Rep2 | FR1 |
| raw | CTRL | Rep3 | FR1,FR2 |
| raw | CTRL | Rep4 | FR1,FR2,FR3 |
| raw | 2min | Rep1 | FR1,FR2,FR5 |
| raw | 2min | Rep2 | FR1,FR2 |
| raw | 2min | Rep3 | FR1,FR2 |
| raw | 2min | Rep4 | FR1,FR2 |
| raw | 8min | Rep1 | FR1,FR2 |
| raw | 8min | Rep2 | FR1,FR2 |
| raw | 8min | Rep3 | FR1,FR2 |
| raw | 8min | Rep4 | FR1,FR2 |
| raw | 20min | Rep1 | FR1,FR2 |
| raw | 20min | Rep2 | FR1 |
| raw | 20min | Rep3 | FR1,FR5 |
| raw | 20min | Rep4 | FR1,FR2 |
| raw | 90min | Rep1 | FR1,FR2,FR5 |
| raw | 90min | Rep2 | FR1,FR2 |
| raw | 90min | Rep3 | FR1,FR2 |
| raw | 90min | Rep4 | FR1,FR2 |
| centered | CTRL | Rep1 | FR1 |
| centered | CTRL | Rep2 | FR1 |
| centered | CTRL | Rep3 | FR1,FR2 |
| centered | CTRL | Rep4 | FR1,FR2,FR3 |
| centered | 2min | Rep1 | FR1,FR2,FR5 |
| centered | 2min | Rep2 | FR1,FR2 |
| centered | 2min | Rep3 | FR1,FR2 |
| centered | 2min | Rep4 | FR1,FR2 |
| centered | 8min | Rep1 | FR1,FR2 |
| centered | 8min | Rep2 | FR1,FR2 |
| centered | 8min | Rep3 | FR1,FR2 |
| centered | 8min | Rep4 | FR1,FR2 |
| centered | 20min | Rep1 | FR1,FR2 |
| centered | 20min | Rep2 | FR1 |
| centered | 20min | Rep3 | FR1,FR5 |
| centered | 20min | Rep4 | FR1,FR2 |
| centered | 90min | Rep1 | FR1,FR2,FR5 |
| centered | 90min | Rep2 | FR1,FR2 |
| centered | 90min | Rep3 | FR1,FR2 |
| centered | 90min | Rep4 | FR1,FR2 |

附表（secondary protein_group）：0 个（无）。

### ④ EGFR 90 min / CTRL 比值表

EGFR 在 80/120 个 run 测到。

gene = EGFR；primary protein_group = `P00533`。ratio_r = total_90min,r / total_CTRL,r，按 rep 编号配对；任一侧 total 缺则留空。n_fr_* = 对应 total 的 n_fractions_present（< 6 即该侧有缺 fraction）。total 为 6 位有效数字，ratio 为 4 位小数，完整精度见 `03_egfr_ratio.tsv`。

**raw**

| rep | total_90min | n_fr_90min | total_CTRL | n_fr_CTRL | ratio (90min/CTRL) | 备注 |
|---|---|---|---|---|---|---|
| Rep1 | 7.57694e+07 | 3 | 1.60857e+08 | 5 | 0.4710 |  |
| Rep2 | 6.67118e+07 | 4 | 2.69796e+08 | 5 | 0.2473 |  |
| Rep3 | 6.79139e+07 | 4 | 1.70902e+08 | 4 | 0.3974 |  |
| Rep4 | 8.58314e+07 | 4 | 5.84781e+07 | 3 | 1.4678 |  |

**centered**

| rep | total_90min | n_fr_90min | total_CTRL | n_fr_CTRL | ratio (90min/CTRL) | 备注 |
|---|---|---|---|---|---|---|
| Rep1 | 7.47593e+07 | 3 | 1.82886e+08 | 5 | 0.4088 |  |
| Rep2 | 7.30967e+07 | 4 | 2.68884e+08 | 5 | 0.2719 |  |
| Rep3 | 7.61634e+07 | 4 | 2.1167e+08 | 4 | 0.3598 |  |
| Rep4 | 8.58346e+07 | 4 | 5.06689e+07 | 3 | 1.6940 |  |

附表（secondary protein_group）：0 个（无）。

## ⑧ 跳过项

- 不复现 DAPAR 的 LOESS 归一化与 KNN/MEC/detQuant 填补（`00_design.md` §1），只用 raw 与 centered 两版；不填补缺失。
- 不复现作者 :58-153 的 limma p 值合并与 :173 的 BH 校正（无 `pval_combi_FDR`）；E-C 只按 MS 阈值列出。
- 脚本在 scratchpad 合成夹具上跑通。

### 自检（方案验收项）

| 检查 | 结果 | 值 |
|---|---|---|
| 输入 01_shares.tsv / 01_groups.tsv 列 = 任务 1 方案列 | PASS |  |
| 01_groups.tsv：每个基因至多一个 is_primary=yes | PASS |  |
| 01_shares.tsv 的 protein_group = 01_groups.tsv 的 primary | PASS | 不一致 0 行 |
| 01_shares.tsv 与 fraction 级重算份额一致（GRB2/SHC1/CBL） | PASS | 比较 308 格，max|差| = 5.0e-13，一侧缺 0 |
| ⑦ 节 = 03_egfr_total.md 的 ③④ 原文（仅标题降一级） | PASS | analysis/08_gate1/03_egfr_total.md |
| §1 引用的 R 原文与 translocation_plots.R 行号逐行一致 | PASS | 12 行 |
| 02_events.tsv 列 = table, gene, normalization, timepoint, key, value | PASS |  |
| 表行数 E-A 90 / E-B 12 / E-C 24 / E-D 24 | PASS | E-A 90 / E-B 12 / E-C 24 / E-D 24 |
| 02_events.tsv 覆盖 E-A..E-D 全部数值（行数 = 4·EA + 4·EB + 5·EC + 6·ED，逐格回读一致） | PASS | 672 / 672，不一致 0 |
| E-A：n / mean / sd(ddof=1) 由 01_shares.tsv 复算一致 | PASS | n 不一致 0，max|Δmean| = 5.0e-13，max|Δsd| = 4.3e-13 |
| E-A：delta_vs_ctrl = mean_tp − mean_CTRL | PASS | max|差| = 1.0e-12 |
| E-B：分母 = 01_shares 中两边 Mem、Cyt 都有值的 rep 数；三个计数复算一致 | PASS | 12/12 行一致 |
| E-C：MS_author 由 fraction 级（6 个 fraction）重复均值复算一致 | PASS | 10 行有值，不一致 0 |
| E-C：MS_max1 / MS_max2 ∈ FR1..FR6，且为 |Δs| 最大 / 次大的 fraction | PASS | 不一致 0 |
| E-C：MS_author 不等于用 3 个区室份额算的同式值 | PASS | 相等 0 |
| E-C：gt_0.1 = (MS_author > 0.1) | PASS | 不一致 0 |
| E-D：t / p_welch 与 Welch 公式手算一致（双侧） | PASS | 不一致 0 |
| 02_events.png：子图数 = 18（2 行 × 9 列），dpi ≥ 150 | PASS | 子图 18，dpi 200，3700×1120 px |
| 02_events.png：有数据的子图每图 4 条重复线 | PASS | 有数据子图 18，'no data' 子图 0 |

