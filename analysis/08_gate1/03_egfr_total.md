# 08_gate1 任务 3：EGFR 总量

生成脚本：`analysis/08_gate1/03_egfr_total.py`（按 `analysis/08_gate1/00_design.md` §0、§2、任务 3 执行）。

## ① 输入

- `analysis/gate1_proteins.tsv`：数据行 344（不含表头）；列 `run, genes, protein_group, quantity`。
- `analysis/07_pilot2/01_run_factors.tsv`：layer=proteome 120 行，可用 (timepoint, fraction, rep) → a 120 个；异常（a 非数或三元组重复）0 个。
- distinct run 字符串 119；解析成功 119，落入 (timepoint, fraction, rep) 设计格 119/120；解析失败 0。
- `genes` 含 EGFR 的行 80；其中 run 解析失败 0、quantity 缺失 0、quantity ≤ 0 0；同一 (protein_group, tp, fraction, rep) 多条同值记录 0 格（计一次），不同值冲突 0 格（不计入，见下）。
- **EGFR 在 80/120 个 run 测到**（任一 EGFR protein_group 的 quantity 有值）。

EGFR protein_group 表：

| gene | protein_group | genes_raw | n_runs_detected | is_primary |
|---|---|---|---|---|
| EGFR | P00533 | EGFR | 80/120 | yes |

## ② 定义

- 匹配：`genes` 按 `;` 拆开（去首尾空白）后任一项等于 `EGFR` 即匹配，记录原始 `genes` 字符串。
- 多个 `protein_group` 全部保留；primary = 测到 run 数最多的 group（测到 = run 解析成功、落入 120 设计格且 quantity 有值），并列取 `protein_group` 字符串排序第一并记录并列；其余为 secondary，列附表。
- run 解析：07 `parse_design`（`code/Protein contour/Zhihan/Test/07_main_analysis.py:32-38`）三个正则 `_(2min|8min|20min|90min|CTRL)_`（忽略大小写，结果统一写成 CTRL/2min/8min/20min/90min）、`_(FR\d)_`、`_(Rep\d)`；fraction ∉ FR1..FR6 或 rep ∉ Rep1..Rep4 计为解析失败。
- `quantity` 转 float；空串 / 非数 / 非有限值视为缺失。
- 两版归一化：`raw` = quantity；`centered` = quantity × 2^a，a 取 `analysis/07_pilot2/01_run_factors.tsv` 中 layer=proteome、同 (timepoint, fraction, rep) 的 `a`（等价于 quantity / 2^(m_07 − anchor)，anchor 为常数）。
- 总量：对每个 (protein_group, timepoint, rep, normalization)，`n_fractions_present` = FR1..FR6 中有值的个数；`total_quantity = Σ_present q_f`；不填补；`n_fractions_present = 0` 时 `total_quantity` 留空。
- 比值：`ratio_r = total_90min,r / total_CTRL,r`，按 rep 编号配对，4 个 rep 各一行；任一侧 total 缺则 ratio 留空（CTRL total = 0 时也留空并注明）。

## ③ EGFR 总量表

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

## ④ EGFR 90 min / CTRL 比值表

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

## ⑤ 跳过项

- 无整节跳过。
- 不复现 DAPAR 的 LOESS 归一化与 KNN/MEC/detQuant 填补（00_design.md §1），只给 raw 与 centered 两版。
- 脚本在 scratchpad 合成夹具上跑通。
