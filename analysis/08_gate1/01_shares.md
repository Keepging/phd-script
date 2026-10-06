# 08 任务 1：GRB2 / SHC1 / CBL / EGFR 区室份额表

脚本：`analysis/08_gate1/01_shares.py`（`python3`，任意目录可复跑，输出确定）。只报计数和数值。

## ① 输入

- `data/pilot/gate1_proteins.tsv`：0 行（不含表头），列 `run, genes, protein_group, quantity`；读法 `pd.read_csv(sep='\t', dtype=str, keep_default_na=False)`，空串 = 缺失。
- `analysis/07_pilot2/01_run_factors.tsv`：240 行，layer=proteome 120 行，有效 (timepoint, fraction, rep) → a：120 个。
- distinct run：0；解析成功 0/120；解析失败 0；覆盖的 (timepoint, fraction, rep) 三元组 0/120；多个 run 映射到同一三元组：0；在 run_factors 中找不到 a 的三元组：0。
- 匹配目标基因的行：0（按基因：GRB2 0、SHC1 0、CBL 0、EGFR 0）；其中 run 解析失败 0、quantity 缺失（空串/非数）0、quantity ≤ 0 0。
- **四个蛋白各在多少 run 测到（/120，任一 protein_group 有有效 quantity 的三元组数）**：GRB2 0/120，SHC1 0/120，CBL 0/120，EGFR 0/120。

group 表（`01_groups.tsv`）：

缺数据：`gate1_proteins.tsv` 0 行，跳过（`01_groups.tsv` 只有表头）。

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

缺数据：`gate1_proteins.tsv` 0 行，跳过。

## ④ 份额表（primary，raw 与 centered 并列）

缺数据：`gate1_proteins.tsv` 0 行，跳过（`01_shares.tsv` 只有表头）。

## ⑤ 附表说明

缺数据：`gate1_proteins.tsv` 0 行，跳过（`01_shares_secondary.tsv`、`01_groups.tsv` 只有表头）。

## ⑥ 跳过项 / 记录

- `data/pilot/gate1_proteins.tsv` 0 行数据：③ 缺 fraction 清单、④ 份额表、⑤ 附表全部跳过；三张 tsv 照常生成，只含表头。数据到位后直接复跑本脚本即可。
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
| 01_shares：Cyt+Mem+Nuc = 1 ± 1e-9（三格非空的行） | 不适用 | 0 行 |
| 01_shares：raw 与 centered 行数相同 | PASS | 0 / 0 |
| 01_shares：n_fractions_present ∈ 1..6 | 不适用 | 0 行 |
| 01_shares_secondary：Cyt+Mem+Nuc = 1 ± 1e-9（三格非空的行） | 不适用 | 0 行 |
| 01_shares_secondary：raw 与 centered 行数相同 | PASS | 0 / 0 |
| 01_shares_secondary：n_fractions_present ∈ 1..6 | 不适用 | 0 行 |
| 0 行输入：01_groups / 01_shares / 01_shares_secondary 只有表头 | PASS | 0 / 0 / 0 |

