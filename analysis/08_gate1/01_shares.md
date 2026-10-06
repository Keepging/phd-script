# 08 任务 1：GRB2 / SHC1 / CBL / EGFR 区室份额表

脚本：`analysis/08_gate1/01_shares.py`（`python3`，任意目录可复跑，输出确定）。只报计数和数值。

## ① 输入

- `analysis/gate1_proteins.tsv`：344 行（不含表头），列 `run, genes, protein_group, quantity`；读法 `pd.read_csv(sep='\t', dtype=str, keep_default_na=False)`，空串 = 缺失。
- `analysis/07_pilot2/01_run_factors.tsv`：240 行，layer=proteome 120 行，有效 (timepoint, fraction, rep) → a：120 个。
- distinct run：119；解析成功 119/120；解析失败 0；覆盖的 (timepoint, fraction, rep) 三元组 119/120；多个 run 映射到同一三元组：0；在 run_factors 中找不到 a 的三元组：0。
- 匹配目标基因的行：344（按基因：GRB2 90、SHC1 115、CBL 59、EGFR 80）；其中 run 解析失败 0、quantity 缺失（空串/非数）0、quantity ≤ 0 0。
- **四个蛋白各在多少 run 测到（/120，任一 protein_group 有有效 quantity 的三元组数）**：GRB2 90/120，SHC1 115/120，CBL 59/120，EGFR 80/120。

group 表（`01_groups.tsv`）：

| gene | protein_group | genes_raw | n_runs_detected | is_primary |
|---|---|---|---|---|
| GRB2 | P62993 | GRB2 | 90 | yes |
| SHC1 | P29353 | SHC1 | 115 | yes |
| CBL | P22681 | CBL | 59 | yes |
| EGFR | P00533 | EGFR | 80 | yes |

primary 并列：0 个基因。

## ② 定义

照抄 `analysis/08_gate1/00_design.md` §2 / 任务 1：

- 蛋白集：GRB2、SHC1、CBL、EGFR。匹配 `genes` 列（分号分隔的多基因按 `;` 拆开后任一等于目标基因即匹配，记录原始 `genes` 字符串）。
- 一个基因对应多个 `protein_group`：全部保留；`primary` = 出现 run 数最多的 group（并列取 `protein_group` 字符串排序第一，并记录并列）；其余为 `secondary`，写附表。
- run 解析：`timepoint ∈ {CTRL,2min,8min,20min,90min}`，`fraction ∈ FR1..FR6`，`rep ∈ Rep1..Rep4`；期望 120 个 proteome run；解析失败的 run 计数并列出。正则为 07 `parse_design`（`code/Protein contour/Zhihan/Test/07_main_analysis.py:32-38`）的 `_(2min|8min|20min|90min|CTRL)_`（re.I）、`_(FR\d)_`、`_(Rep\d)`。
- `quantity` 转 float；空串/非数视为缺失。
- 两版归一化：`raw` = quantity；`centered` = quantity × 2^a（a 取 `01_run_factors.tsv` 中 layer=proteome、同 run 的 `a`）。run 名匹配用 (timepoint, fraction, rep) 三元组。
- 份额：对每个 (gene, protein_group, timepoint, rep)：`n_fractions_present` = 6 个 fraction 中有值的个数；`share_f = q_f / Σ_{present} q_f`；`Cyt = share_FR1 + share_FR2`、`Mem = share_FR3 + share_FR4`、`Nuc = share_FR5 + share_FR6`；某区室两个 fraction 都缺 → 该区室格留空；`n_fractions_present < 6` 的行单列"缺 fraction 清单"；不填补。raw、centered 各算一次（`normalization` 列）。
- `01_shares.tsv`（primary）/ `01_shares_secondary.tsv`（其余 group）列：`gene, protein_group, timepoint, rep, normalization, n_fractions_present, Cyt, Mem, Nuc`；某 (tp, rep) 一个 fraction 都没有则不出行，并在 md 记录。
- `01_groups.tsv` 列：`gene, protein_group, genes_raw, n_runs_detected, is_primary`（每个 gene×group 一行；gene 没测到则一行 `n_runs_detected = 0`）。

脚本补充规则（方案未写明处，非替换定义）：

- timepoint 正则带 re.I，匹配到的写法统一成 `CTRL/2min/8min/20min/90min`；`FR\d`、`Rep\d` 匹配到但不在 FR1..FR6 / Rep1..Rep4 的 run 记为解析失败（原因"超出设计范围"）。
- 基因 token 去首尾空白后与目标基因做区分大小写的完全相等比较；同一行 genes 中重复的同一目标基因只计一次。
- `nan`/`inf` 字符串按"非数"处理为缺失；quantity ≤ 0 保留原值参与计算，计数见 ①，若有份额落在 [0,1] 外在 ⑥ 列出；某 (tp, rep) 有值 fraction 之和 = 0 时三格留空并在 ⑥ 列出。
- `n_runs_detected` = 该 gene×group 有 ≥1 条有效 quantity 的 (timepoint, fraction, rep) 三元组数（只计解析成功的 run）；基因级"测到 run 数"= 该基因任一 group 有有效 quantity 的三元组数。
- `genes_raw`：该 gene×group 出现过的全部原始 `genes` 字符串，排序后以 ` | ` 连接；`is_primary` 取 `yes`/`no`（gene 没测到的那一行留空）。
- 同一 (gene, protein_group, timepoint, fraction, rep) 有多条有效 quantity：值全相同取该值；值不同记为冲突，该格当缺失并在 ⑥ 列出。
- 方案两处对 0 行输入的要求不一致（`01_groups.tsv`"gene 没测到则一行 n_runs_detected = 0" vs 验收"0 行输入时所有表只有表头"）：本脚本在输入 0 行时按验收写只有表头的 `01_groups.tsv`，四个基因的 0/120 写在 ①；输入非 0 行时，没测到的基因各写一行 `n_runs_detected = 0`。
- tsv 中份额保留 12 位小数；md 表中保留 4 位。行序：gene（GRB2, SHC1, CBL, EGFR）→ protein_group → timepoint（CTRL, 2min, 8min, 20min, 90min）→ rep → normalization（raw, centered）。

## ③ 缺 fraction 清单

`n_fractions_present < 6` 的行：126（primary 126、secondary 0）。

| gene | protein_group | role | timepoint | rep | normalization | n_fractions_present | 缺的 fraction | 留空区室 |
|---|---|---|---|---|---|---|---|---|
| GRB2 | P62993 | primary | CTRL | Rep1 | raw | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | CTRL | Rep1 | centered | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | CTRL | Rep2 | raw | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | CTRL | Rep2 | centered | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | CTRL | Rep3 | raw | 3 | FR4,FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | CTRL | Rep3 | centered | 3 | FR4,FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 2min | Rep1 | raw | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 2min | Rep1 | centered | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 2min | Rep2 | raw | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 2min | Rep2 | centered | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 2min | Rep3 | raw | 5 | FR5 |  |
| GRB2 | P62993 | primary | 2min | Rep3 | centered | 5 | FR5 |  |
| GRB2 | P62993 | primary | 2min | Rep4 | raw | 5 | FR5 |  |
| GRB2 | P62993 | primary | 2min | Rep4 | centered | 5 | FR5 |  |
| GRB2 | P62993 | primary | 8min | Rep1 | raw | 5 | FR5 |  |
| GRB2 | P62993 | primary | 8min | Rep1 | centered | 5 | FR5 |  |
| GRB2 | P62993 | primary | 8min | Rep2 | raw | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 8min | Rep2 | centered | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 8min | Rep4 | raw | 5 | FR5 |  |
| GRB2 | P62993 | primary | 8min | Rep4 | centered | 5 | FR5 |  |
| GRB2 | P62993 | primary | 20min | Rep1 | raw | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 20min | Rep1 | centered | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 20min | Rep2 | raw | 5 | FR5 |  |
| GRB2 | P62993 | primary | 20min | Rep2 | centered | 5 | FR5 |  |
| GRB2 | P62993 | primary | 20min | Rep3 | raw | 5 | FR5 |  |
| GRB2 | P62993 | primary | 20min | Rep3 | centered | 5 | FR5 |  |
| GRB2 | P62993 | primary | 20min | Rep4 | raw | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 20min | Rep4 | centered | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 90min | Rep1 | raw | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 90min | Rep1 | centered | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 90min | Rep2 | raw | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 90min | Rep2 | centered | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 90min | Rep3 | raw | 5 | FR5 |  |
| GRB2 | P62993 | primary | 90min | Rep3 | centered | 5 | FR5 |  |
| GRB2 | P62993 | primary | 90min | Rep4 | raw | 4 | FR5,FR6 | Nuc |
| GRB2 | P62993 | primary | 90min | Rep4 | centered | 4 | FR5,FR6 | Nuc |
| SHC1 | P29353 | primary | CTRL | Rep1 | raw | 5 | FR5 |  |
| SHC1 | P29353 | primary | CTRL | Rep1 | centered | 5 | FR5 |  |
| SHC1 | P29353 | primary | CTRL | Rep2 | raw | 5 | FR5 |  |
| SHC1 | P29353 | primary | CTRL | Rep2 | centered | 5 | FR5 |  |
| SHC1 | P29353 | primary | CTRL | Rep3 | raw | 5 | FR5 |  |
| SHC1 | P29353 | primary | CTRL | Rep3 | centered | 5 | FR5 |  |
| SHC1 | P29353 | primary | 90min | Rep1 | raw | 5 | FR5 |  |
| SHC1 | P29353 | primary | 90min | Rep1 | centered | 5 | FR5 |  |
| SHC1 | P29353 | primary | 90min | Rep4 | raw | 5 | FR6 |  |
| SHC1 | P29353 | primary | 90min | Rep4 | centered | 5 | FR6 |  |
| CBL | P22681 | primary | CTRL | Rep1 | raw | 3 | FR4,FR5,FR6 | Nuc |
| CBL | P22681 | primary | CTRL | Rep1 | centered | 3 | FR4,FR5,FR6 | Nuc |
| CBL | P22681 | primary | CTRL | Rep2 | raw | 3 | FR4,FR5,FR6 | Nuc |
| CBL | P22681 | primary | CTRL | Rep2 | centered | 3 | FR4,FR5,FR6 | Nuc |
| CBL | P22681 | primary | CTRL | Rep3 | raw | 3 | FR4,FR5,FR6 | Nuc |
| CBL | P22681 | primary | CTRL | Rep3 | centered | 3 | FR4,FR5,FR6 | Nuc |
| CBL | P22681 | primary | CTRL | Rep4 | raw | 4 | FR5,FR6 | Nuc |
| CBL | P22681 | primary | CTRL | Rep4 | centered | 4 | FR5,FR6 | Nuc |
| CBL | P22681 | primary | 2min | Rep1 | raw | 3 | FR2,FR3,FR5 |  |
| CBL | P22681 | primary | 2min | Rep1 | centered | 3 | FR2,FR3,FR5 |  |
| CBL | P22681 | primary | 2min | Rep2 | raw | 4 | FR5,FR6 | Nuc |
| CBL | P22681 | primary | 2min | Rep2 | centered | 4 | FR5,FR6 | Nuc |
| CBL | P22681 | primary | 2min | Rep3 | raw | 4 | FR3,FR6 |  |
| CBL | P22681 | primary | 2min | Rep3 | centered | 4 | FR3,FR6 |  |
| CBL | P22681 | primary | 2min | Rep4 | raw | 3 | FR2,FR3,FR6 |  |
| CBL | P22681 | primary | 2min | Rep4 | centered | 3 | FR2,FR3,FR6 |  |
| CBL | P22681 | primary | 8min | Rep1 | raw | 3 | FR1,FR2,FR3 | Cyt |
| CBL | P22681 | primary | 8min | Rep1 | centered | 3 | FR1,FR2,FR3 | Cyt |
| CBL | P22681 | primary | 8min | Rep2 | raw | 3 | FR1,FR5,FR6 | Nuc |
| CBL | P22681 | primary | 8min | Rep2 | centered | 3 | FR1,FR5,FR6 | Nuc |
| CBL | P22681 | primary | 8min | Rep3 | raw | 3 | FR2,FR3,FR5 |  |
| CBL | P22681 | primary | 8min | Rep3 | centered | 3 | FR2,FR3,FR5 |  |
| CBL | P22681 | primary | 8min | Rep4 | raw | 5 | FR5 |  |
| CBL | P22681 | primary | 8min | Rep4 | centered | 5 | FR5 |  |
| CBL | P22681 | primary | 20min | Rep1 | raw | 2 | FR2,FR3,FR5,FR6 | Nuc |
| CBL | P22681 | primary | 20min | Rep1 | centered | 2 | FR2,FR3,FR5,FR6 | Nuc |
| CBL | P22681 | primary | 20min | Rep2 | raw | 4 | FR2,FR5 |  |
| CBL | P22681 | primary | 20min | Rep2 | centered | 4 | FR2,FR5 |  |
| CBL | P22681 | primary | 20min | Rep3 | raw | 3 | FR2,FR4,FR5 |  |
| CBL | P22681 | primary | 20min | Rep3 | centered | 3 | FR2,FR4,FR5 |  |
| CBL | P22681 | primary | 20min | Rep4 | raw | 2 | FR2,FR3,FR5,FR6 | Nuc |
| CBL | P22681 | primary | 20min | Rep4 | centered | 2 | FR2,FR3,FR5,FR6 | Nuc |
| CBL | P22681 | primary | 90min | Rep1 | raw | 2 | FR2,FR4,FR5,FR6 | Nuc |
| CBL | P22681 | primary | 90min | Rep1 | centered | 2 | FR2,FR4,FR5,FR6 | Nuc |
| CBL | P22681 | primary | 90min | Rep2 | raw | 3 | FR4,FR5,FR6 | Nuc |
| CBL | P22681 | primary | 90min | Rep2 | centered | 3 | FR4,FR5,FR6 | Nuc |
| CBL | P22681 | primary | 90min | Rep3 | raw | 1 | FR2,FR3,FR4,FR5,FR6 | Mem,Nuc |
| CBL | P22681 | primary | 90min | Rep3 | centered | 1 | FR2,FR3,FR4,FR5,FR6 | Mem,Nuc |
| CBL | P22681 | primary | 90min | Rep4 | raw | 1 | FR2,FR3,FR4,FR5,FR6 | Mem,Nuc |
| CBL | P22681 | primary | 90min | Rep4 | centered | 1 | FR2,FR3,FR4,FR5,FR6 | Mem,Nuc |
| EGFR | P00533 | primary | CTRL | Rep1 | raw | 5 | FR1 |  |
| EGFR | P00533 | primary | CTRL | Rep1 | centered | 5 | FR1 |  |
| EGFR | P00533 | primary | CTRL | Rep2 | raw | 5 | FR1 |  |
| EGFR | P00533 | primary | CTRL | Rep2 | centered | 5 | FR1 |  |
| EGFR | P00533 | primary | CTRL | Rep3 | raw | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | CTRL | Rep3 | centered | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | CTRL | Rep4 | raw | 3 | FR1,FR2,FR3 | Cyt |
| EGFR | P00533 | primary | CTRL | Rep4 | centered | 3 | FR1,FR2,FR3 | Cyt |
| EGFR | P00533 | primary | 2min | Rep1 | raw | 3 | FR1,FR2,FR5 | Cyt |
| EGFR | P00533 | primary | 2min | Rep1 | centered | 3 | FR1,FR2,FR5 | Cyt |
| EGFR | P00533 | primary | 2min | Rep2 | raw | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 2min | Rep2 | centered | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 2min | Rep3 | raw | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 2min | Rep3 | centered | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 2min | Rep4 | raw | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 2min | Rep4 | centered | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 8min | Rep1 | raw | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 8min | Rep1 | centered | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 8min | Rep2 | raw | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 8min | Rep2 | centered | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 8min | Rep3 | raw | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 8min | Rep3 | centered | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 8min | Rep4 | raw | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 8min | Rep4 | centered | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 20min | Rep1 | raw | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 20min | Rep1 | centered | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 20min | Rep2 | raw | 5 | FR1 |  |
| EGFR | P00533 | primary | 20min | Rep2 | centered | 5 | FR1 |  |
| EGFR | P00533 | primary | 20min | Rep3 | raw | 4 | FR1,FR5 |  |
| EGFR | P00533 | primary | 20min | Rep3 | centered | 4 | FR1,FR5 |  |
| EGFR | P00533 | primary | 20min | Rep4 | raw | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 20min | Rep4 | centered | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 90min | Rep1 | raw | 3 | FR1,FR2,FR5 | Cyt |
| EGFR | P00533 | primary | 90min | Rep1 | centered | 3 | FR1,FR2,FR5 | Cyt |
| EGFR | P00533 | primary | 90min | Rep2 | raw | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 90min | Rep2 | centered | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 90min | Rep3 | raw | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 90min | Rep3 | centered | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 90min | Rep4 | raw | 4 | FR1,FR2 | Cyt |
| EGFR | P00533 | primary | 90min | Rep4 | centered | 4 | FR1,FR2 | Cyt |


一个 fraction 都没有、因而不出行的 (gene, group, tp, rep, normalization)：0 条。

## ④ 份额表（primary，raw 与 centered 并列）

`01_shares.tsv`：160 行（raw 80、centered 80；数据齐全时应为 160）。按基因：GRB2 40、SHC1 40、CBL 40、EGFR 40。

| gene | protein_group | timepoint | rep | n_fr raw | Cyt raw | Mem raw | Nuc raw | n_fr centered | Cyt centered | Mem centered | Nuc centered |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GRB2 | P62993 | CTRL | Rep1 | 4 | 0.8372 | 0.1628 |  | 4 | 0.7670 | 0.2330 |  |
| GRB2 | P62993 | CTRL | Rep2 | 4 | 0.8376 | 0.1624 |  | 4 | 0.7613 | 0.2387 |  |
| GRB2 | P62993 | CTRL | Rep3 | 3 | 0.9264 | 0.0736 |  | 3 | 0.8656 | 0.1344 |  |
| GRB2 | P62993 | CTRL | Rep4 | 6 | 0.5018 | 0.3502 | 0.1480 | 6 | 0.4894 | 0.3602 | 0.1504 |
| GRB2 | P62993 | 2min | Rep1 | 4 | 0.3454 | 0.6546 |  | 4 | 0.2695 | 0.7305 |  |
| GRB2 | P62993 | 2min | Rep2 | 4 | 0.2627 | 0.7373 |  | 4 | 0.2269 | 0.7731 |  |
| GRB2 | P62993 | 2min | Rep3 | 5 | 0.4017 | 0.5693 | 0.0290 | 5 | 0.3231 | 0.6239 | 0.0530 |
| GRB2 | P62993 | 2min | Rep4 | 5 | 0.3719 | 0.5860 | 0.0422 | 5 | 0.3212 | 0.6133 | 0.0655 |
| GRB2 | P62993 | 8min | Rep1 | 5 | 0.2737 | 0.6874 | 0.0389 | 5 | 0.2382 | 0.6889 | 0.0730 |
| GRB2 | P62993 | 8min | Rep2 | 4 | 0.3421 | 0.6579 |  | 4 | 0.2805 | 0.7195 |  |
| GRB2 | P62993 | 8min | Rep3 | 6 | 0.2719 | 0.6407 | 0.0873 | 6 | 0.2205 | 0.6603 | 0.1192 |
| GRB2 | P62993 | 8min | Rep4 | 5 | 0.4068 | 0.4971 | 0.0961 | 5 | 0.2908 | 0.5602 | 0.1490 |
| GRB2 | P62993 | 20min | Rep1 | 4 | 0.4164 | 0.5836 |  | 4 | 0.3325 | 0.6675 |  |
| GRB2 | P62993 | 20min | Rep2 | 5 | 0.3830 | 0.5534 | 0.0636 | 5 | 0.3124 | 0.5954 | 0.0922 |
| GRB2 | P62993 | 20min | Rep3 | 5 | 0.4156 | 0.4486 | 0.1358 | 5 | 0.2851 | 0.5357 | 0.1792 |
| GRB2 | P62993 | 20min | Rep4 | 4 | 0.5515 | 0.4485 |  | 4 | 0.4449 | 0.5551 |  |
| GRB2 | P62993 | 90min | Rep1 | 4 | 0.7286 | 0.2714 |  | 4 | 0.6649 | 0.3351 |  |
| GRB2 | P62993 | 90min | Rep2 | 4 | 0.6497 | 0.3503 |  | 4 | 0.5789 | 0.4211 |  |
| GRB2 | P62993 | 90min | Rep3 | 5 | 0.6466 | 0.3122 | 0.0411 | 5 | 0.5720 | 0.3700 | 0.0580 |
| GRB2 | P62993 | 90min | Rep4 | 4 | 0.6613 | 0.3387 |  | 4 | 0.5897 | 0.4103 |  |
| SHC1 | P29353 | CTRL | Rep1 | 5 | 0.8608 | 0.1226 | 0.0166 | 5 | 0.8094 | 0.1572 | 0.0334 |
| SHC1 | P29353 | CTRL | Rep2 | 5 | 0.7554 | 0.2068 | 0.0377 | 5 | 0.6992 | 0.2377 | 0.0631 |
| SHC1 | P29353 | CTRL | Rep3 | 5 | 0.7654 | 0.2202 | 0.0145 | 5 | 0.6571 | 0.3162 | 0.0267 |
| SHC1 | P29353 | CTRL | Rep4 | 6 | 0.5502 | 0.3213 | 0.1285 | 6 | 0.5215 | 0.3424 | 0.1360 |
| SHC1 | P29353 | 2min | Rep1 | 6 | 0.5168 | 0.4328 | 0.0504 | 6 | 0.4280 | 0.4809 | 0.0911 |
| SHC1 | P29353 | 2min | Rep2 | 6 | 0.4367 | 0.4907 | 0.0726 | 6 | 0.3792 | 0.5134 | 0.1074 |
| SHC1 | P29353 | 2min | Rep3 | 6 | 0.4906 | 0.4305 | 0.0790 | 6 | 0.4079 | 0.4695 | 0.1226 |
| SHC1 | P29353 | 2min | Rep4 | 6 | 0.4332 | 0.5081 | 0.0587 | 6 | 0.3745 | 0.5455 | 0.0800 |
| SHC1 | P29353 | 8min | Rep1 | 6 | 0.4570 | 0.4697 | 0.0733 | 6 | 0.4061 | 0.4858 | 0.1080 |
| SHC1 | P29353 | 8min | Rep2 | 6 | 0.4456 | 0.4796 | 0.0748 | 6 | 0.3689 | 0.5265 | 0.1045 |
| SHC1 | P29353 | 8min | Rep3 | 6 | 0.5001 | 0.4141 | 0.0858 | 6 | 0.4212 | 0.4579 | 0.1209 |
| SHC1 | P29353 | 8min | Rep4 | 6 | 0.5143 | 0.3631 | 0.1226 | 6 | 0.3835 | 0.4288 | 0.1878 |
| SHC1 | P29353 | 20min | Rep1 | 6 | 0.6292 | 0.3051 | 0.0657 | 6 | 0.5170 | 0.3745 | 0.1085 |
| SHC1 | P29353 | 20min | Rep2 | 6 | 0.5409 | 0.3804 | 0.0788 | 6 | 0.4676 | 0.4136 | 0.1187 |
| SHC1 | P29353 | 20min | Rep3 | 6 | 0.5684 | 0.2392 | 0.1925 | 6 | 0.4429 | 0.2636 | 0.2935 |
| SHC1 | P29353 | 20min | Rep4 | 6 | 0.6628 | 0.2709 | 0.0663 | 6 | 0.5498 | 0.3365 | 0.1137 |
| SHC1 | P29353 | 90min | Rep1 | 5 | 0.6395 | 0.3174 | 0.0431 | 5 | 0.5482 | 0.3751 | 0.0767 |
| SHC1 | P29353 | 90min | Rep2 | 6 | 0.6726 | 0.2587 | 0.0687 | 6 | 0.5891 | 0.3114 | 0.0995 |
| SHC1 | P29353 | 90min | Rep3 | 6 | 0.6604 | 0.2554 | 0.0842 | 6 | 0.5818 | 0.3024 | 0.1158 |
| SHC1 | P29353 | 90min | Rep4 | 5 | 0.7001 | 0.2667 | 0.0332 | 5 | 0.6483 | 0.3072 | 0.0446 |
| CBL | P22681 | CTRL | Rep1 | 3 | 0.7521 | 0.2479 |  | 3 | 0.6365 | 0.3635 |  |
| CBL | P22681 | CTRL | Rep2 | 3 | 0.8998 | 0.1002 |  | 3 | 0.8319 | 0.1681 |  |
| CBL | P22681 | CTRL | Rep3 | 3 | 0.8203 | 0.1797 |  | 3 | 0.7001 | 0.2999 |  |
| CBL | P22681 | CTRL | Rep4 | 4 | 0.0077 | 0.9923 |  | 4 | 0.0079 | 0.9921 |  |
| CBL | P22681 | 2min | Rep1 | 3 | 0.1496 | 0.5956 | 0.2548 | 3 | 0.0969 | 0.4893 | 0.4138 |
| CBL | P22681 | 2min | Rep2 | 4 | 0.1303 | 0.8697 |  | 4 | 0.0969 | 0.9031 |  |
| CBL | P22681 | 2min | Rep3 | 4 | 0.1137 | 0.6563 | 0.2300 | 4 | 0.0888 | 0.6361 | 0.2751 |
| CBL | P22681 | 2min | Rep4 | 3 | 0.0954 | 0.6716 | 0.2329 | 3 | 0.0785 | 0.6422 | 0.2793 |
| CBL | P22681 | 8min | Rep1 | 3 |  | 0.8201 | 0.1799 | 3 |  | 0.7474 | 0.2526 |
| CBL | P22681 | 8min | Rep2 | 3 | 0.0868 | 0.9132 |  | 3 | 0.0554 | 0.9446 |  |
| CBL | P22681 | 8min | Rep3 | 3 | 0.2454 | 0.6468 | 0.1077 | 3 | 0.1870 | 0.6326 | 0.1803 |
| CBL | P22681 | 8min | Rep4 | 5 | 0.1300 | 0.8042 | 0.0658 | 5 | 0.0723 | 0.8478 | 0.0799 |
| CBL | P22681 | 20min | Rep1 | 2 | 0.1521 | 0.8479 |  | 2 | 0.1178 | 0.8822 |  |
| CBL | P22681 | 20min | Rep2 | 4 | 0.0321 | 0.8950 | 0.0729 | 4 | 0.0169 | 0.9098 | 0.0732 |
| CBL | P22681 | 20min | Rep3 | 3 | 0.0691 | 0.7198 | 0.2111 | 3 | 0.0326 | 0.7615 | 0.2059 |
| CBL | P22681 | 20min | Rep4 | 2 | 0.2657 | 0.7343 |  | 2 | 0.2021 | 0.7979 |  |
| CBL | P22681 | 90min | Rep1 | 2 | 0.9356 | 0.0644 |  | 2 | 0.8800 | 0.1200 |  |
| CBL | P22681 | 90min | Rep2 | 3 | 0.0577 | 0.9423 |  | 3 | 0.0386 | 0.9614 |  |
| CBL | P22681 | 90min | Rep3 | 1 | 1.0000 |  |  | 1 | 1.0000 |  |  |
| CBL | P22681 | 90min | Rep4 | 1 | 1.0000 |  |  | 1 | 1.0000 |  |  |
| EGFR | P00533 | CTRL | Rep1 | 5 | 0.0313 | 0.7565 | 0.2122 | 5 | 0.0199 | 0.6736 | 0.3065 |
| EGFR | P00533 | CTRL | Rep2 | 5 | 0.0372 | 0.6428 | 0.3200 | 5 | 0.0263 | 0.5473 | 0.4264 |
| EGFR | P00533 | CTRL | Rep3 | 4 |  | 0.6698 | 0.3302 | 4 |  | 0.5608 | 0.4392 |
| EGFR | P00533 | CTRL | Rep4 | 3 |  | 0.8373 | 0.1627 | 3 |  | 0.8288 | 0.1712 |
| EGFR | P00533 | 2min | Rep1 | 3 |  | 0.8079 | 0.1921 | 3 |  | 0.6914 | 0.3086 |
| EGFR | P00533 | 2min | Rep2 | 4 |  | 0.6607 | 0.3393 | 4 |  | 0.5586 | 0.4414 |
| EGFR | P00533 | 2min | Rep3 | 4 |  | 0.5817 | 0.4183 | 4 |  | 0.4760 | 0.5240 |
| EGFR | P00533 | 2min | Rep4 | 4 |  | 0.5622 | 0.4378 | 4 |  | 0.4678 | 0.5322 |
| EGFR | P00533 | 8min | Rep1 | 4 |  | 0.6858 | 0.3142 | 4 |  | 0.5918 | 0.4082 |
| EGFR | P00533 | 8min | Rep2 | 4 |  | 0.7672 | 0.2328 | 4 |  | 0.6689 | 0.3311 |
| EGFR | P00533 | 8min | Rep3 | 4 |  | 0.7575 | 0.2425 | 4 |  | 0.7204 | 0.2796 |
| EGFR | P00533 | 8min | Rep4 | 4 |  | 0.5070 | 0.4930 | 4 |  | 0.4075 | 0.5925 |
| EGFR | P00533 | 20min | Rep1 | 4 |  | 0.9012 | 0.0988 | 4 |  | 0.8761 | 0.1239 |
| EGFR | P00533 | 20min | Rep2 | 5 | 0.0084 | 0.5586 | 0.4330 | 5 | 0.0069 | 0.4828 | 0.5104 |
| EGFR | P00533 | 20min | Rep3 | 4 | 0.0044 | 0.4854 | 0.5102 | 4 | 0.0029 | 0.4115 | 0.5857 |
| EGFR | P00533 | 20min | Rep4 | 4 |  | 0.9127 | 0.0873 | 4 |  | 0.8803 | 0.1197 |
| EGFR | P00533 | 90min | Rep1 | 3 |  | 0.9563 | 0.0437 | 3 |  | 0.9297 | 0.0703 |
| EGFR | P00533 | 90min | Rep2 | 4 |  | 0.8706 | 0.1294 | 4 |  | 0.8374 | 0.1626 |
| EGFR | P00533 | 90min | Rep3 | 4 |  | 0.8830 | 0.1170 | 4 |  | 0.8667 | 0.1333 |
| EGFR | P00533 | 90min | Rep4 | 4 |  | 0.8746 | 0.1254 | 4 |  | 0.8377 | 0.1623 |

## ⑤ 附表说明

- `01_groups.tsv`：每个 gene×protein_group 一行，4 行；列 `gene, protein_group, genes_raw, n_runs_detected, is_primary`。
- `01_shares_secondary.tsv`：非 primary group 的份额，列同 `01_shares.tsv`，0 行。
- secondary group：0 个。

## ⑥ 跳过项 / 记录

- 解析失败 run：0；多 run 同三元组：0；找不到 a 的三元组：0。
- 值冲突的 (gene, group, tp, fraction, rep) 格：0。
- 有值 fraction 之和 = 0（三格留空）的行：0。
- 份额落在 [0,1] 外的行：0。
- 不复现 DAPAR 的 LOESS 归一化与 KNN 填补（方案 §1），只用 raw 与 centered 两版；不填补缺失。
- 代码冒烟测试：脚本在 scratchpad 合成夹具上跑通（夹具数据与结果均未写入本目录）。

### 自检（方案验收项）

| 检查 | 结果 | 值 |
|---|---|---|
| 01_groups.tsv 列 = 方案列 | PASS |  |
| 01_shares.tsv / 01_shares_secondary.tsv 列 = 方案列 | PASS |  |
| run_factors proteome 行 = 120，(tp, fr, rep) 唯一且 a 全为数值 | PASS | 120 行 / 有效 a 120 |
| 01_shares：Cyt+Mem+Nuc = 1 ± 1e-9（三格非空的行） | PASS | 80 行，max|偏差| = 1.00e-12 |
| 01_shares：有区室留空的行，非空区室之和 = 1 ± 1e-9 | PASS | 80 行 |
| 01_shares：raw 与 centered 行数相同 | PASS | 80 / 80 |
| 01_shares：n_fractions_present ∈ 1..6 | PASS | 1..6 |
| 01_shares_secondary：Cyt+Mem+Nuc = 1 ± 1e-9（三格非空的行） | 不适用 | 0 行 |
| 01_shares_secondary：raw 与 centered 行数相同 | PASS | 0 / 0 |
| 01_shares_secondary：n_fractions_present ∈ 1..6 | 不适用 | 0 行 |
| 0 行输入：所有表只有表头 | 不适用 | 输入 344 行 |
| 01_shares.tsv 行数 = 160（数据齐全时） | PASS | 160 |
| 每个目标基因在 01_groups.tsv 至少一行 | PASS |  |
| 每个测到的基因恰一个 primary | PASS |  |

