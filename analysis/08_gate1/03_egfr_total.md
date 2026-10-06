# 08_gate1 任务 3：EGFR 总量

生成脚本：`analysis/08_gate1/03_egfr_total.py`（按 `analysis/08_gate1/00_design.md` §0、§2、任务 3 执行）。

## ① 输入

- `data/pilot/gate1_proteins.tsv`：数据行 0（不含表头）；列 `run, genes, protein_group, quantity`。
- `analysis/07_pilot2/01_run_factors.tsv`：layer=proteome 120 行，可用 (timepoint, fraction, rep) → a 120 个；异常（a 非数或三元组重复）0 个。
- distinct run 字符串 0；解析成功 0，落入 (timepoint, fraction, rep) 设计格 0/120；解析失败 0。
- `genes` 含 EGFR 的行 0；其中 run 解析失败 0、quantity 缺失 0、quantity ≤ 0 0；同一 (protein_group, tp, fraction, rep) 多条同值记录 0 格（计一次），不同值冲突 0 格（不计入，见下）。
- **EGFR 在 0/120 个 run 测到**（任一 EGFR protein_group 的 quantity 有值）。

EGFR protein_group 表：

缺数据：`gate1_proteins.tsv` 0 行，跳过（无 EGFR protein_group）。

## ② 定义

- 匹配：`genes` 按 `;` 拆开（去首尾空白）后任一项等于 `EGFR` 即匹配，记录原始 `genes` 字符串。
- 多个 `protein_group` 全部保留；primary = 测到 run 数最多的 group（测到 = run 解析成功、落入 120 设计格且 quantity 有值），并列取 `protein_group` 字符串排序第一并记录并列；其余为 secondary，列附表。
- run 解析：07 `parse_design`（`code/Protein contour/Zhihan/Test/07_main_analysis.py:32-38`）三个正则 `_(2min|8min|20min|90min|CTRL)_`（忽略大小写，结果统一写成 CTRL/2min/8min/20min/90min）、`_(FR\d)_`、`_(Rep\d)`；fraction ∉ FR1..FR6 或 rep ∉ Rep1..Rep4 计为解析失败。
- `quantity` 转 float；空串 / 非数 / 非有限值视为缺失。
- 两版归一化：`raw` = quantity；`centered` = quantity × 2^a，a 取 `analysis/07_pilot2/01_run_factors.tsv` 中 layer=proteome、同 (timepoint, fraction, rep) 的 `a`（等价于 quantity / 2^(m_07 − anchor)，anchor 为常数）。
- 总量：对每个 (protein_group, timepoint, rep, normalization)，`n_fractions_present` = FR1..FR6 中有值的个数；`total_quantity = Σ_present q_f`；不填补；`n_fractions_present = 0` 时 `total_quantity` 留空。
- 比值：`ratio_r = total_90min,r / total_CTRL,r`，按 rep 编号配对，4 个 rep 各一行；任一侧 total 缺则 ratio 留空（CTRL total = 0 时也留空并注明）。

## ③ EGFR 总量表

EGFR 在 0/120 个 run 测到。

缺数据：`gate1_proteins.tsv` 0 行，跳过

## ④ EGFR 90 min / CTRL 比值表

EGFR 在 0/120 个 run 测到。

缺数据：`gate1_proteins.tsv` 0 行，跳过

## ⑤ 跳过项

- 缺数据：`gate1_proteins.tsv` 0 行，跳过。③ 总量表、④ 比值表均无数值；`03_egfr_total.tsv`、`03_egfr_ratio.tsv` 只有表头。数据到位后直接复跑本脚本即可。
- 不复现 DAPAR 的 LOESS 归一化与 KNN/MEC/detQuant 填补（00_design.md §1），只给 raw 与 centered 两版。
- 脚本在 scratchpad 合成夹具上跑通。
