# 02 归一化趋势：run 级 median_log2

## ① 输入

| 文件 | 行数(不含表头) | 列 |
|---|---|---|
| `data/pilot/run_medians.tsv` | 240 | `run`, `n_precursors`, `median_log2` |

- 解析成功 240 行；解析失败 0 行。
- `phospho`：120 行。
- `proteome`：120 行。
- 输出：`02_run_medians.tsv`（240 行）、`02_run_medians_summary.tsv`（60 行）、`02_run_medians.png`、`02_run_medians_phospho.png`、`02_run_medians_proteome.png`、`02_run_medians.py`。

## ② 定义

（照抄 `analysis/06_pilot/00_design.md` 任务 2「解析」「偏移定义」两节）

**解析（照 07 `parse_design`，`:32-38`）**
- `layer` = `run` 在第一个 `/` 之前的部分（phospho | proteome）。
- `timepoint` = 正则 `_(2min|8min|20min|90min|CTRL)_`（忽略大小写）；`fraction` = `_(FR\d)_`；`rep` = `_(Rep\d)`。解析失败的行计数并列出（预期 0）。
- 预期 2 layer × 5 tp × 6 FR × 4 rep = 240，每个 cell 恰 4 行。

**偏移定义（两种都算，都只是计数）**
- `anchor_{layer}` = 该 layer 120 个 `median_log2` 的中位数（07 `factors()` `:239-242` 的 anchor）。`offset_anchor` = `median_log2 − anchor_layer`。
- `offset_ctrl` = `median_log2 − median_log2(同 layer、同 fraction、同 rep 的 CTRL run)`；CTRL 行本身留空。
- "方向一致" := 同一 layer×fraction×timepoint 的 4 个重复的 offset 符号全同（全 >0 或全 <0；出现 0 记为不一致，并单独计数）。

**summary 列**：sd 用样本标准差（ddof=1）。CTRL 行的 `*_ctrl` 列留空。

执行细节（不改定义，仅说明实现）：
- 中位数用 07 `median()`（`:100-101`）的同一写法；偶数个取中间两值的平均。
- 正负号在未取整的 float 差值上判断；`n_zero_*` 计 offset 恰等于 0 的行。
- `consistent_*` 取值 `yes` / `no`；CTRL 行 `consistent_ctrl` 为空。
- 输出数值保留 6 位小数（`mean_n_precursors` 2 位）；`median_log2`、`n_precursors` 原样照抄输入字符串。

## ③ anchor 值

| layer | n_runs | anchor_layer |
|---|---|---|
| phospho | 120 | 18.845850 |
| proteome | 120 | 22.671000 |

## ④ summary 表（全文，60 行；与 `02_run_medians_summary.tsv` 相同）

| layer | fraction | timepoint | n_reps | mean_median_log2 | sd_median_log2 | min_median_log2 | max_median_log2 | mean_n_precursors | n_pos_anchor | n_neg_anchor | n_zero_anchor | consistent_anchor | n_pos_ctrl | n_neg_ctrl | n_zero_ctrl | consistent_ctrl |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| phospho | FR1 | CTRL | 4 | 19.818225 | 0.935866 | 18.886700 | 20.936700 | 5400.75 | 4 | 0 | 0 | yes |  |  |  |  |
| phospho | FR1 | 2min | 4 | 20.874600 | 0.286563 | 20.647100 | 21.287300 | 5811.75 | 4 | 0 | 0 | yes | 3 | 1 | 0 | no |
| phospho | FR1 | 8min | 4 | 20.632125 | 0.223107 | 20.376200 | 20.912300 | 6061.75 | 4 | 0 | 0 | yes | 3 | 1 | 0 | no |
| phospho | FR1 | 20min | 4 | 20.320600 | 0.702456 | 19.575000 | 21.184900 | 5064.00 | 4 | 0 | 0 | yes | 3 | 1 | 0 | no |
| phospho | FR1 | 90min | 4 | 20.036925 | 0.330312 | 19.578900 | 20.360900 | 4886.75 | 4 | 0 | 0 | yes | 2 | 2 | 0 | no |
| phospho | FR2 | CTRL | 4 | 19.354500 | 0.675947 | 18.421800 | 19.885200 | 6495.25 | 3 | 1 | 0 | no |  |  |  |  |
| phospho | FR2 | 2min | 4 | 19.532275 | 0.203011 | 19.292100 | 19.788800 | 6783.75 | 4 | 0 | 0 | yes | 2 | 2 | 0 | no |
| phospho | FR2 | 8min | 4 | 19.214125 | 0.599885 | 18.319700 | 19.566400 | 6127.00 | 3 | 1 | 0 | no | 1 | 3 | 0 | no |
| phospho | FR2 | 20min | 4 | 19.516000 | 0.124746 | 19.365300 | 19.647500 | 5782.00 | 4 | 0 | 0 | yes | 2 | 2 | 0 | no |
| phospho | FR2 | 90min | 4 | 19.482850 | 0.553253 | 18.692200 | 19.982100 | 5747.75 | 3 | 1 | 0 | no | 2 | 2 | 0 | no |
| phospho | FR3 | CTRL | 4 | 17.395550 | 0.551194 | 16.854500 | 18.097700 | 3674.25 | 0 | 4 | 0 | yes |  |  |  |  |
| phospho | FR3 | 2min | 4 | 16.997850 | 0.279040 | 16.734000 | 17.348100 | 3032.25 | 0 | 4 | 0 | yes | 1 | 3 | 0 | no |
| phospho | FR3 | 8min | 4 | 17.203425 | 0.298463 | 16.918100 | 17.498000 | 3440.25 | 0 | 4 | 0 | yes | 1 | 3 | 0 | no |
| phospho | FR3 | 20min | 4 | 17.120825 | 0.248013 | 16.791400 | 17.390700 | 3335.25 | 0 | 4 | 0 | yes | 2 | 2 | 0 | no |
| phospho | FR3 | 90min | 4 | 17.925350 | 0.453043 | 17.521800 | 18.542500 | 3051.00 | 0 | 4 | 0 | yes | 3 | 1 | 0 | no |
| phospho | FR4 | CTRL | 4 | 18.645375 | 0.715608 | 17.609000 | 19.236000 | 5106.00 | 2 | 2 | 0 | no |  |  |  |  |
| phospho | FR4 | 2min | 4 | 18.646750 | 0.208946 | 18.428000 | 18.861000 | 5272.50 | 1 | 3 | 0 | no | 2 | 2 | 0 | no |
| phospho | FR4 | 8min | 4 | 18.521200 | 0.331526 | 18.035700 | 18.782000 | 4592.75 | 0 | 4 | 0 | yes | 1 | 3 | 0 | no |
| phospho | FR4 | 20min | 4 | 18.557500 | 0.352738 | 18.261700 | 19.065200 | 4438.75 | 1 | 3 | 0 | no | 2 | 2 | 0 | no |
| phospho | FR4 | 90min | 4 | 18.753500 | 0.270100 | 18.479600 | 19.095300 | 4556.25 | 1 | 3 | 0 | no | 3 | 1 | 0 | no |
| phospho | FR5 | CTRL | 4 | 18.991650 | 0.693195 | 18.178700 | 19.815200 | 5856.25 | 2 | 2 | 0 | no |  |  |  |  |
| phospho | FR5 | 2min | 4 | 19.137950 | 0.188015 | 18.951000 | 19.311900 | 6162.00 | 4 | 0 | 0 | yes | 3 | 1 | 0 | no |
| phospho | FR5 | 8min | 4 | 19.487650 | 0.432666 | 19.114500 | 20.078100 | 5852.50 | 4 | 0 | 0 | yes | 3 | 1 | 0 | no |
| phospho | FR5 | 20min | 4 | 19.321975 | 0.399594 | 18.746300 | 19.599700 | 5550.75 | 3 | 1 | 0 | no | 3 | 1 | 0 | no |
| phospho | FR5 | 90min | 4 | 19.489750 | 0.283126 | 19.228300 | 19.836400 | 5733.75 | 4 | 0 | 0 | yes | 3 | 1 | 0 | no |
| phospho | FR6 | CTRL | 4 | 18.385000 | 0.997594 | 17.705000 | 19.867400 | 4387.75 | 1 | 3 | 0 | no |  |  |  |  |
| phospho | FR6 | 2min | 4 | 18.020900 | 0.131434 | 17.871000 | 18.167700 | 4419.50 | 0 | 4 | 0 | yes | 2 | 2 | 0 | no |
| phospho | FR6 | 8min | 4 | 17.857550 | 0.328178 | 17.518800 | 18.284700 | 3897.25 | 0 | 4 | 0 | yes | 2 | 2 | 0 | no |
| phospho | FR6 | 20min | 4 | 18.037500 | 0.335128 | 17.555900 | 18.333800 | 4131.00 | 0 | 4 | 0 | yes | 3 | 1 | 0 | no |
| phospho | FR6 | 90min | 4 | 18.364375 | 0.435469 | 17.791000 | 18.805900 | 4064.50 | 0 | 4 | 0 | yes | 3 | 1 | 0 | no |
| proteome | FR1 | CTRL | 4 | 23.060475 | 0.055954 | 22.976600 | 23.091300 | 24577.00 | 4 | 0 | 0 | yes |  |  |  |  |
| proteome | FR1 | 2min | 4 | 23.089175 | 0.073224 | 23.018900 | 23.190800 | 22882.75 | 4 | 0 | 0 | yes | 2 | 2 | 0 | no |
| proteome | FR1 | 8min | 4 | 23.055800 | 0.081728 | 22.962400 | 23.161600 | 23059.25 | 4 | 0 | 0 | yes | 1 | 3 | 0 | no |
| proteome | FR1 | 20min | 4 | 23.086450 | 0.086949 | 22.988500 | 23.166400 | 23561.50 | 4 | 0 | 0 | yes | 2 | 2 | 0 | no |
| proteome | FR1 | 90min | 4 | 23.045375 | 0.042778 | 22.997200 | 23.089600 | 23595.25 | 4 | 0 | 0 | yes | 2 | 2 | 0 | no |
| proteome | FR2 | CTRL | 4 | 22.942100 | 0.279674 | 22.553100 | 23.164400 | 34900.25 | 3 | 1 | 0 | no |  |  |  |  |
| proteome | FR2 | 2min | 4 | 22.987500 | 0.149338 | 22.801700 | 23.159600 | 34157.50 | 4 | 0 | 0 | yes | 3 | 1 | 0 | no |
| proteome | FR2 | 8min | 4 | 22.804950 | 0.108752 | 22.691800 | 22.932500 | 33104.50 | 4 | 0 | 0 | yes | 1 | 3 | 0 | no |
| proteome | FR2 | 20min | 4 | 22.829150 | 0.125871 | 22.666400 | 22.950500 | 33251.50 | 3 | 1 | 0 | no | 1 | 3 | 0 | no |
| proteome | FR2 | 90min | 4 | 22.844275 | 0.127086 | 22.730000 | 22.975700 | 32354.50 | 4 | 0 | 0 | yes | 1 | 3 | 0 | no |
| proteome | FR3 | CTRL | 4 | 22.293050 | 0.268216 | 22.053200 | 22.676700 | 26587.75 | 1 | 3 | 0 | no |  |  |  |  |
| proteome | FR3 | 2min | 4 | 21.952300 | 0.138788 | 21.806600 | 22.141000 | 22130.00 | 0 | 4 | 0 | yes | 0 | 4 | 0 | yes |
| proteome | FR3 | 8min | 4 | 21.919800 | 0.047642 | 21.865500 | 21.973500 | 21511.50 | 0 | 4 | 0 | yes | 0 | 4 | 0 | yes |
| proteome | FR3 | 20min | 4 | 21.947475 | 0.061991 | 21.875100 | 22.020000 | 23336.25 | 0 | 4 | 0 | yes | 0 | 4 | 0 | yes |
| proteome | FR3 | 90min | 4 | 22.142525 | 0.154416 | 22.046500 | 22.372900 | 24567.75 | 0 | 4 | 0 | yes | 1 | 3 | 0 | no |
| proteome | FR4 | CTRL | 4 | 22.837950 | 0.142249 | 22.696100 | 23.012400 | 32624.50 | 4 | 0 | 0 | yes |  |  |  |  |
| proteome | FR4 | 2min | 4 | 22.813675 | 0.117206 | 22.715800 | 22.975300 | 30088.00 | 4 | 0 | 0 | yes | 1 | 3 | 0 | no |
| proteome | FR4 | 8min | 4 | 22.704375 | 0.025318 | 22.675000 | 22.735700 | 29239.25 | 4 | 0 | 0 | yes | 1 | 3 | 0 | no |
| proteome | FR4 | 20min | 4 | 22.731900 | 0.065344 | 22.641800 | 22.795500 | 29102.75 | 3 | 1 | 0 | no | 1 | 3 | 0 | no |
| proteome | FR4 | 90min | 4 | 22.773425 | 0.140575 | 22.642400 | 22.928200 | 31636.25 | 2 | 2 | 0 | no | 1 | 3 | 0 | no |
| proteome | FR5 | CTRL | 4 | 22.627400 | 0.211595 | 22.363700 | 22.849500 | 37097.00 | 2 | 2 | 0 | no |  |  |  |  |
| proteome | FR5 | 2min | 4 | 22.412950 | 0.094063 | 22.281800 | 22.493100 | 34459.75 | 0 | 4 | 0 | yes | 0 | 4 | 0 | yes |
| proteome | FR5 | 8min | 4 | 22.471900 | 0.086023 | 22.370400 | 22.554200 | 34237.25 | 0 | 4 | 0 | yes | 1 | 3 | 0 | no |
| proteome | FR5 | 20min | 4 | 22.338400 | 0.072572 | 22.250200 | 22.424100 | 33256.00 | 0 | 4 | 0 | yes | 1 | 3 | 0 | no |
| proteome | FR5 | 90min | 4 | 22.566000 | 0.097788 | 22.503300 | 22.710700 | 36386.75 | 1 | 3 | 0 | no | 1 | 3 | 0 | no |
| proteome | FR6 | CTRL | 4 | 22.213250 | 0.383745 | 21.913400 | 22.741200 | 30474.75 | 1 | 3 | 0 | no |  |  |  |  |
| proteome | FR6 | 2min | 4 | 21.957125 | 0.182277 | 21.770200 | 22.137200 | 29837.75 | 0 | 4 | 0 | yes | 0 | 4 | 0 | yes |
| proteome | FR6 | 8min | 4 | 21.865800 | 0.075745 | 21.756700 | 21.928100 | 27598.25 | 0 | 4 | 0 | yes | 0 | 4 | 0 | yes |
| proteome | FR6 | 20min | 4 | 21.961175 | 0.080001 | 21.860200 | 22.049500 | 28291.75 | 0 | 4 | 0 | yes | 2 | 2 | 0 | no |
| proteome | FR6 | 90min | 4 | 22.073150 | 0.124382 | 21.966600 | 22.250600 | 28644.75 | 0 | 4 | 0 | yes | 2 | 2 | 0 | no |

## ⑤ 一致性计数

每个 layer×fraction：`consistent_anchor=yes` 的时间点数（/5，含 CTRL）、`consistent_ctrl=yes` 的时间点数（/4，仅 EGF 时间点）；`zero_cells_*` = 该 layer×fraction 中出现 ≥1 个 offset=0 的时间点数（这些 cell 记为不一致）。

| layer | fraction | consistent_anchor=yes | timepoints(anchor) | consistent_ctrl=yes | timepoints(ctrl) | zero_cells_anchor | zero_cells_ctrl |
|---|---|---|---|---|---|---|---|
| phospho | FR1 | 5/5 | CTRL,2min,8min,20min,90min | 0/4 | — | 0 | 0 |
| phospho | FR2 | 2/5 | 2min,20min | 0/4 | — | 0 | 0 |
| phospho | FR3 | 5/5 | CTRL,2min,8min,20min,90min | 0/4 | — | 0 | 0 |
| phospho | FR4 | 1/5 | 8min | 0/4 | — | 0 | 0 |
| phospho | FR5 | 3/5 | 2min,8min,90min | 0/4 | — | 0 | 0 |
| phospho | FR6 | 4/5 | 2min,8min,20min,90min | 0/4 | — | 0 | 0 |
| proteome | FR1 | 5/5 | CTRL,2min,8min,20min,90min | 0/4 | — | 0 | 0 |
| proteome | FR2 | 3/5 | 2min,8min,90min | 0/4 | — | 0 | 0 |
| proteome | FR3 | 4/5 | 2min,8min,20min,90min | 3/4 | 2min,8min,20min | 0 | 0 |
| proteome | FR4 | 3/5 | CTRL,2min,8min | 0/4 | — | 0 | 0 |
| proteome | FR5 | 3/5 | 2min,8min,20min | 1/4 | 2min | 0 | 0 |
| proteome | FR6 | 4/5 | 2min,8min,20min,90min | 2/4 | 2min,8min | 0 | 0 |

| 范围 | consistent_anchor=yes | consistent_ctrl=yes | zero_cells_anchor | zero_cells_ctrl |
|---|---|---|---|---|
| phospho | 20/30 | 0/24 | 0 | 0 |
| proteome | 22/30 | 6/24 | 0 | 0 |
| 总计 | 42/60 | 6/48 | 0 | 0 |

- offset_anchor 恰为 0 的行数：0（/240）；offset_ctrl 恰为 0 的行数：0（/192）。

## ⑥ 解析失败列表 / 跳过 / 缺失项

- 解析失败：0 行（列表为空）。
- 缺 CTRL 对照导致 offset_ctrl 无法计算的 EGF 行：0。
- 重复的 CTRL 键（layer, fraction, rep）：0。
- n_reps ≠ 4 的 cell：0。
