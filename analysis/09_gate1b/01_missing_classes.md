# 09_gate1b 任务 1：缺失核查（Q2b）

脚本：`analysis/09_gate1b/01_missing_check.py`（`python3 analysis/09_gate1b/01_missing_check.py` 可复跑，输出确定）。pandas 2.2.3。逐行结果：`analysis/09_gate1b/01_missing_classes.tsv`（480 行）。

## ① 输入

| 文件 | 用途 | 行数 / 说明 |
|---|---|---|
| `data/pilot/gate1_missing_check.tsv` | n_precursors、min_PG_qvalue、min_precursor_qvalue | 480 行 = 4 蛋白 × 120 run；n_precursors == 0 的行 113（GRB2 27 / SHC1 4 / CBL 49 / EGFR 33），这些行两列 q 值均为空串 |
| `analysis/gate1_proteins.tsv` | A 类判定（有 quantity） | 344 行；GRB2 90（P62993） / SHC1 115（P29353） / CBL 59（P22681） / EGFR 80（P00533） |
| `data/pilot/proteome_pg_matrix_120.tsv` | 交叉核对 | 10000 行 × 122 列；四基因各 1 行 |
| `code/Protein contour/Zhihan/Test/07_main_analysis.py:32-38` | run 名解析 | `parse_design` 三个正则 `_(2min\|8min\|20min\|90min\|CTRL)_`（re.I）、`_(FR\d)_`、`_(Rep\d)`；脚本照抄，并核对源文件该段文本一致 |

- run 解析：120 个 run 全部可解析，解析失败 0；(timepoint, fraction, rep) 共 120 种 = 5 × 6 × 4，与 run 一一对应。
- 与 pg_matrix 交叉核对（四基因行的非空 run 集合 vs 本表 A 类 run 集合；数值为 gate1_proteins quantity 与 pg_matrix 值的最大相对差）：

| gene | pg_matrix Protein.Group | pg_matrix 非空 run 数 | A 类 run 数 | 集合一致 | 最大相对差 |
|---|---|---|---|---|---|
| GRB2 | P62993 | 90 | 90 | 是 | 0.00e+00 |
| SHC1 | P29353 | 115 | 115 | 是 | 0.00e+00 |
| CBL | P22681 | 59 | 59 | 是 | 0.00e+00 |
| EGFR | P00533 | 80 | 80 | 是 | 0.00e+00 |

- `data/pilot/gate1_proteins.tsv`：现为 344 行，与 `analysis/gate1_proteins.tsv` 逐字节相同；00_design §0 记为 0 行，与现状不符。本轮按任务书只用 `analysis/gate1_proteins.tsv`。

## ② 定义

对每个 (gene ∈ {GRB2, SHC1, CBL, EGFR}, run ∈ 120)：

- **A**：`analysis/gate1_proteins.tsv` 中该 (gene, run) 的 `quantity` 非空（gene 匹配：`genes` 按 `;` 拆开，任一 token 等于目标基因）。
- **B**：非 A 且 `gate1_missing_check.tsv` 的 `n_precursors > 0`。
- **C**：非 A 且 `n_precursors == 0`。
- B 类按 `min_PG_qvalue` 拆分：≤ 0.01 / > 0.01（`min_PG_qvalue` 转 float 比较；恰等于 0.01 的 B 行 0 个）。
- 读表 `sep='\t', dtype=str, keep_default_na=False`，空串 = 缺失；tsv 中 `n_precursors`、q 值、`pg_quantity` 均为原字符串（非 A 的 `pg_quantity` 为空串）。
- tsv 行序：gene（GRB2, SHC1, CBL, EGFR）→ timepoint（CTRL, 2min, 8min, 20min, 90min）→ fraction → rep。

## ③ 每蛋白 A/B/C 计数

| gene | A | B | C | 合计 | advisor 参考 A/B/C | 一致 |
|---|---|---|---|---|---|---|
| GRB2 | 90 | 3 | 27 | 120 | 90/3/27 | 是 |
| SHC1 | 115 | 1 | 4 | 120 | 115/1/4 | 是 |
| CBL | 59 | 12 | 49 | 120 | 59/12/49 | 是 |
| EGFR | 80 | 7 | 33 | 120 | 80/7/33 | 是 |
| 合计 | 344 | 23 | 113 | 480 |  |  |

- A 类中 `n_precursors == 0` 的行数：0（A 类 `n_precursors` 最小值 1）。

## ④ 每蛋白 × fraction 的 A/B/C

**GRB2**

| fraction | A | B | C | 合计 |
|---|---|---|---|---|
| FR1 | 20 | 0 | 0 | 20 |
| FR2 | 20 | 0 | 0 | 20 |
| FR3 | 20 | 0 | 0 | 20 |
| FR4 | 19 | 0 | 1 | 20 |
| FR5 | 2 | 0 | 18 | 20 |
| FR6 | 9 | 3 | 8 | 20 |
| 合计 | 90 | 3 | 27 | 120 |

**SHC1**

| fraction | A | B | C | 合计 |
|---|---|---|---|---|
| FR1 | 20 | 0 | 0 | 20 |
| FR2 | 20 | 0 | 0 | 20 |
| FR3 | 20 | 0 | 0 | 20 |
| FR4 | 20 | 0 | 0 | 20 |
| FR5 | 16 | 1 | 3 | 20 |
| FR6 | 19 | 0 | 1 | 20 |
| 合计 | 115 | 1 | 4 | 120 |

**CBL**

| fraction | A | B | C | 合计 |
|---|---|---|---|---|
| FR1 | 18 | 2 | 0 | 20 |
| FR2 | 9 | 2 | 9 | 20 |
| FR3 | 11 | 3 | 6 | 20 |
| FR4 | 12 | 0 | 8 | 20 |
| FR5 | 3 | 1 | 16 | 20 |
| FR6 | 6 | 4 | 10 | 20 |
| 合计 | 59 | 12 | 49 | 120 |

**EGFR**

| fraction | A | B | C | 合计 |
|---|---|---|---|---|
| FR1 | 0 | 1 | 19 | 20 |
| FR2 | 4 | 4 | 12 | 20 |
| FR3 | 19 | 0 | 1 | 20 |
| FR4 | 20 | 0 | 0 | 20 |
| FR5 | 17 | 2 | 1 | 20 |
| FR6 | 20 | 0 | 0 | 20 |
| 合计 | 80 | 7 | 33 | 120 |

## ⑤ B 类清单

B 类共 23 行。按 `min_PG_qvalue` 拆分：

| gene | B | min_PG_qvalue ≤ 0.01 | min_PG_qvalue > 0.01 |
|---|---|---|---|
| GRB2 | 3 | 0 | 3 |
| SHC1 | 1 | 0 | 1 |
| CBL | 12 | 1 | 11 |
| EGFR | 7 | 0 | 7 |
| 合计 | 23 | 1 | 22 |

清单（行序同 tsv）：

| gene | run | n_precursors | min_PG_qvalue | min_precursor_qvalue |
|---|---|---|---|---|
| GRB2 | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_8min_FR6_Rep2 | 1 | 0.01452969666570425 | 0.0019097927724942565 |
| GRB2 | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_20min_FR6_Rep1 | 1 | 0.01370906364172697 | 0.0017660752637311816 |
| GRB2 | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_20min_FR6_Rep4 | 1 | 0.015033695846796036 | 0.0021566024515777826 |
| SHC1 | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_90min_FR5_Rep1 | 1 | 0.012748222798109055 | 0.0005662870244123042 |
| CBL | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_2min_FR2_Rep1 | 1 | 0.013234269805252552 | 0.0006000599823892117 |
| CBL | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_2min_FR6_Rep2 | 1 | 0.028565557673573494 | 0.0017698535230010748 |
| CBL | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_2min_FR6_Rep4 | 1 | 0.07257365435361862 | 0.005842665210366249 |
| CBL | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_8min_FR1_Rep1 | 1 | 0.009316770359873772 | 0.0005066996673122048 |
| CBL | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_8min_FR1_Rep2 | 1 | 0.013557483442127705 | 0.0007145997951738536 |
| CBL | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_8min_FR5_Rep3 | 1 | 0.011372251436114311 | 0.006838476285338402 |
| CBL | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_8min_FR6_Rep2 | 1 | 0.029719853773713112 | 0.0017582597211003304 |
| CBL | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_20min_FR2_Rep1 | 1 | 0.011085973121225834 | 0.000495154585223645 |
| CBL | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_20min_FR3_Rep1 | 1 | 0.060092806816101074 | 0.00609272625297308 |
| CBL | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_20min_FR3_Rep4 | 1 | 0.06696730852127075 | 0.006113321986049414 |
| CBL | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_90min_FR3_Rep3 | 1 | 0.06800000369548798 | 0.006688389461487532 |
| CBL | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_90min_FR6_Rep3 | 1 | 0.03382616490125656 | 0.0021903770975768566 |
| EGFR | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR1_Rep4 | 1 | 0.04904306307435036 | 0.003141406225040555 |
| EGFR | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR2_Rep3 | 1 | 0.011505047790706158 | 0.0021681119687855244 |
| EGFR | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_2min_FR2_Rep4 | 1 | 0.027561837807297707 | 0.0012558638118207455 |
| EGFR | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_2min_FR5_Rep1 | 1 | 0.012588943354785442 | 0.003258973825722933 |
| EGFR | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_8min_FR2_Rep1 | 1 | 0.015246188268065453 | 0.0006119451718404889 |
| EGFR | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_20min_FR2_Rep1 | 1 | 0.018033867701888084 | 0.00389629858545959 |
| EGFR | 20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_90min_FR5_Rep1 | 1 | 0.02005251869559288 | 0.00230617169290781 |

## ⑥ 跳过项

- 无跳过项：三个输入文件均存在且结构符合方案。
- 自检（脚本内全部执行）：

| 检查 | 结果 | 细节 |
|---|---|---|
| 07_main_analysis.py:32-38 为 parse_design，三个正则与本脚本照抄一致 | 通过 | tp = re.search(r"_(2min\|8min\|20min\|90min\|CTRL)_", run, re.I); fr = re.search(r"_(FR\d)_", run); rep = re.search(r"_(Rep\d)", run) |
| gate1_missing_check.tsv 列名 | 通过 | ['run', 'gene', 'n_precursors', 'min_PG_qvalue', 'min_precursor_qvalue'] |
| gate1_proteins.tsv 列名 | 通过 | ['run', 'genes', 'protein_group', 'quantity'] |
| gate1_missing_check.tsv 480 行 | 通过 | 480 |
| analysis/gate1_proteins.tsv 344 行 | 通过 | 344 |
| missing_check 120 个 run 且全部可解析 | 通过 | n_run=120, 解析失败=0 |
| run 解析为 5 tp × 6 FR × 4 rep 且一一对应 | 通过 | 不同 (tp,FR,rep) = 120 |
| missing_check 的 gene 集合 = 四蛋白，每基因 120 行 | 通过 | {'GRB2': 120, 'SHC1': 120, 'CBL': 120, 'EGFR': 120} |
| missing_check (run, gene) 无重复 | 通过 | 0 |
| n_precursors 全为非负整数 | 通过 |  |
| n_precursors == 0 的行两列 q 值均为空串；> 0 的行均非空 | 通过 | n_precursors==0 行数 113 |
| gate1_proteins 每行都匹配到四蛋白之一，(gene, run) 无重复 | 通过 | 未匹配 0 行，重复 0 个 |
| gate1_proteins 的 run 都在 missing_check 的 120 run 内 | 通过 | 多出 0 |
| gate1_proteins quantity 全部非空且可转为正数 | 通过 |  |
| 所有 gate1_proteins 的 (gene, run) 都落入 A | 通过 | A=344, gate1_proteins=344 |
| 输出 480 行、列顺序正确 | 通过 | 480 |
| GRB2 A+B+C = 120 | 通过 | 120 |
| GRB2 A/B/C 与 advisor 参考一致 | 通过 | 本轮 (90, 3, 27) vs 参考 (90, 3, 27) |
| SHC1 A+B+C = 120 | 通过 | 120 |
| SHC1 A/B/C 与 advisor 参考一致 | 通过 | 本轮 (115, 1, 4) vs 参考 (115, 1, 4) |
| CBL A+B+C = 120 | 通过 | 120 |
| CBL A/B/C 与 advisor 参考一致 | 通过 | 本轮 (59, 12, 49) vs 参考 (59, 12, 49) |
| EGFR A+B+C = 120 | 通过 | 120 |
| EGFR A/B/C 与 advisor 参考一致 | 通过 | 本轮 (80, 7, 33) vs 参考 (80, 7, 33) |
| A 类中 n_precursors == 0 的行数 = 0 | 通过 | 0 |
| B 类 min_PG_qvalue ≤ 0.01 计数与 advisor 参考一致（CBL 1，其余 0） | 通过 | {'GRB2': 0, 'SHC1': 0, 'CBL': 1, 'EGFR': 0} |
| 每 (蛋白, fraction) 行合计 = 20 | 通过 | [20] |
| pg_matrix 的 120 个 run 列与 missing_check 的 run 集合一致 | 通过 | pg_matrix run 列 120 |
| pg_matrix GRB2 行（P62993）非空 run 集合 = A 类 run 集合 | 通过 | pg_matrix 非空 90，A 90，对称差 0 |
| pg_matrix GRB2 行与 gate1_proteins 的 protein_group 一致 | 通过 | ['P62993'] vs P62993 |
| pg_matrix SHC1 行（P29353）非空 run 集合 = A 类 run 集合 | 通过 | pg_matrix 非空 115，A 115，对称差 0 |
| pg_matrix SHC1 行与 gate1_proteins 的 protein_group 一致 | 通过 | ['P29353'] vs P29353 |
| pg_matrix CBL 行（P22681）非空 run 集合 = A 类 run 集合 | 通过 | pg_matrix 非空 59，A 59，对称差 0 |
| pg_matrix CBL 行与 gate1_proteins 的 protein_group 一致 | 通过 | ['P22681'] vs P22681 |
| pg_matrix EGFR 行（P00533）非空 run 集合 = A 类 run 集合 | 通过 | pg_matrix 非空 80，A 80，对称差 0 |
| pg_matrix EGFR 行与 gate1_proteins 的 protein_group 一致 | 通过 | ['P00533'] vs P00533 |

