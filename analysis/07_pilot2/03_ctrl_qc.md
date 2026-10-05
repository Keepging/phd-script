# 03 CTRL 质检（run 级稳健 z）

脚本：`analysis/07_pilot2/03_ctrl_qc.py`（确定性，可复跑）。只报计数和数值。

## ① 输入

- `data/pilot/run_medians.tsv`：240 行（不含表头），列 `run, n_precursors, median_log2`；读法 `pd.read_csv(sep='\t', dtype=str, keep_default_na=False)`，空串 = 缺失（本表空串 0 个）。
- 解析：`layer` = `run` 第一个 `/` 前的前缀（phospho | proteome）；timepoint / fraction / rep 用 `code/Protein contour/Zhihan/Test/07_main_analysis.py` `parse_design`（`:32-38`）的三个正则 `_(2min|8min|20min|90min|CTRL)_`（re.I）、`_(FR\d)_`、`_(Rep\d)`。
- 解析失败：0 行。layer 计数：phospho 120、proteome 120；CTRL run：48（phospho 24、proteome 24）。
- **与方案不符**：方案写 24 个 CTRL run（phospho 12 + proteome 12）；输入中 CTRL run 为 48 个（2 layer × 6 fraction × 4 rep）。本表输出全部 48 个，未做子集挑选。方案验收参考值（25 条、phospho 9 / proteome 16、4 个 z）均在这 48 个 run 上复现，见自检。

## ② 定义

- 组（cell）= 同 layer × fraction 的 20 个 run（5 timepoint × 4 rep，**含非 CTRL**）；共 2 × 6 = 12 个 cell。
- 对 `n_precursors` 和 `median_log2` 分别在每个 cell 的 20 个值上算：`med = median(20 值)`，`MAD = median(|x − med|)`，`z = (x − med) / (1.4826 × MAD)`；中位数用 `statistics.median`（偶数个取中间两值均值）。
- MAD = 0 时 z 留空（本数据 MAD = 0 的 (cell, 指标) 数见 ⑥）。
- 只输出 CTRL run（本数据 48 个：phospho 24、proteome 24）。`flag_n` / `flag_log2`：|z| > 2 → yes，否则 no（判定用未取整的 z；z 为空时 flag 留空）。
- 数值列保留 4 位小数；`n_precursors`、`median_log2` 两列为输入原字符串。

### cell 统计（`03_ctrl_qc_cells.tsv`，12 行）

| layer | fraction | n_runs | median_n | mad_n | median_log2 | mad_log2 |
|---|---|---|---|---|---|---|
| phospho | FR1 | 20 | 5490.5000 | 500.0000 | 20.4646 | 0.3901 |
| phospho | FR2 | 20 | 6273.0000 | 553.0000 | 19.5388 | 0.1462 |
| phospho | FR3 | 20 | 3301.0000 | 486.5000 | 17.2626 | 0.2881 |
| phospho | FR4 | 20 | 4854.0000 | 607.0000 | 18.6335 | 0.2013 |
| phospho | FR5 | 20 | 5883.0000 | 485.5000 | 19.2894 | 0.2918 |
| phospho | FR6 | 20 | 4319.0000 | 324.0000 | 18.0560 | 0.2340 |
| proteome | FR1 | 20 | 23237.0000 | 654.5000 | 23.0663 | 0.0357 |
| proteome | FR2 | 20 | 33255.5000 | 1714.5000 | 22.9138 | 0.1166 |
| proteome | FR3 | 20 | 23118.0000 | 1763.0000 | 21.9967 | 0.0805 |
| proteome | FR4 | 20 | 29762.0000 | 737.0000 | 22.7379 | 0.0602 |
| proteome | FR5 | 20 | 34965.5000 | 1316.0000 | 22.4800 | 0.0785 |
| proteome | FR6 | 20 | 28825.5000 | 963.5000 | 21.9562 | 0.0946 |

## ③ 48 个 CTRL run（`03_ctrl_qc.tsv`，48 行）

| layer | fraction | rep | run | n_precursors | median_log2 | cell_median_n | cell_mad_n | z_n_precursors | cell_median_log2 | cell_mad_log2 | z_median_log2 | flag_n | flag_log2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| phospho | FR1 | Rep1 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR1_Rep1/` | 6446 | 19.2291 | 5490.5000 | 500.0000 | 1.2890 | 20.4646 | 0.3901 | -2.1366 | no | yes |
| phospho | FR1 | Rep2 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR1_Rep2/` | 4593 | 20.2204 | 5490.5000 | 500.0000 | -1.2107 | 20.4646 | 0.3901 | -0.4224 | no | no |
| phospho | FR1 | Rep3 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR1_Rep3/` | 5382 | 20.9367 | 5490.5000 | 500.0000 | -0.1464 | 20.4646 | 0.3901 | 0.8163 | no | no |
| phospho | FR1 | Rep4 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR1_Rep4/` | 5182 | 18.8867 | 5490.5000 | 500.0000 | -0.4162 | 20.4646 | 0.3901 | -2.7287 | no | yes |
| phospho | FR2 | Rep1 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR2_Rep1/` | 7483 | 19.2925 | 6273.0000 | 553.0000 | 1.4758 | 19.5388 | 0.1462 | -1.1363 | no | no |
| phospho | FR2 | Rep2 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR2_Rep2/` | 6346 | 19.8852 | 6273.0000 | 553.0000 | 0.0890 | 19.5388 | 0.1462 | 1.5981 | no | no |
| phospho | FR2 | Rep3 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR2_Rep3/` | 6539 | 19.8185 | 6273.0000 | 553.0000 | 0.3244 | 19.5388 | 0.1462 | 1.2904 | no | no |
| phospho | FR2 | Rep4 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR2_Rep4/` | 5613 | 18.4218 | 6273.0000 | 553.0000 | -0.8050 | 19.5388 | 0.1462 | -5.1533 | no | yes |
| phospho | FR3 | Rep1 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR3_Rep1/` | 4156 | 16.8545 | 3301.0000 | 486.5000 | 1.1854 | 17.2626 | 0.2881 | -0.9556 | no | no |
| phospho | FR3 | Rep2 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR3_Rep2/` | 2699 | 17.0776 | 3301.0000 | 486.5000 | -0.8346 | 17.2626 | 0.2881 | -0.4332 | no | no |
| phospho | FR3 | Rep3 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR3_Rep3/` | 3423 | 17.5524 | 3301.0000 | 486.5000 | 0.1691 | 17.2626 | 0.2881 | 0.6786 | no | no |
| phospho | FR3 | Rep4 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR3_Rep4/` | 4419 | 18.0977 | 3301.0000 | 486.5000 | 1.5500 | 17.2626 | 0.2881 | 1.9554 | no | no |
| phospho | FR4 | Rep1 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR4_Rep1/` | 6053 | 18.7845 | 4854.0000 | 607.0000 | 1.3323 | 18.6335 | 0.2013 | 0.5057 | no | no |
| phospho | FR4 | Rep2 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR4_Rep2/` | 3057 | 17.6090 | 4854.0000 | 607.0000 | -1.9968 | 18.6335 | 0.2013 | -3.4321 | no | yes |
| phospho | FR4 | Rep3 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR4_Rep3/` | 5420 | 18.9520 | 4854.0000 | 607.0000 | 0.6289 | 18.6335 | 0.2013 | 1.0668 | no | no |
| phospho | FR4 | Rep4 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR4_Rep4/` | 5894 | 19.2360 | 4854.0000 | 607.0000 | 1.1556 | 18.6335 | 0.2013 | 2.0181 | no | yes |
| phospho | FR5 | Rep1 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR5_Rep1/` | 5745 | 18.1787 | 5883.0000 | 485.5000 | -0.1917 | 19.2894 | 0.2918 | -2.5679 | no | yes |
| phospho | FR5 | Rep2 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR5_Rep2/` | 5029 | 18.7601 | 5883.0000 | 485.5000 | -1.1864 | 19.2894 | 0.2918 | -1.2238 | no | no |
| phospho | FR5 | Rep3 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR5_Rep3/` | 6229 | 19.2126 | 5883.0000 | 485.5000 | 0.4807 | 19.2894 | 0.2918 | -0.1777 | no | no |
| phospho | FR5 | Rep4 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR5_Rep4/` | 6422 | 19.8152 | 5883.0000 | 485.5000 | 0.7488 | 19.2894 | 0.2918 | 1.2155 | no | no |
| phospho | FR6 | Rep1 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR6_Rep1/` | 5208 | 18.0271 | 4319.0000 | 324.0000 | 1.8507 | 18.0560 | 0.2340 | -0.0833 | no | no |
| phospho | FR6 | Rep2 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR6_Rep2/` | 2777 | 17.7050 | 4319.0000 | 324.0000 | -3.2101 | 18.0560 | 0.2340 | -1.0120 | yes | no |
| phospho | FR6 | Rep3 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR6_Rep3/` | 4193 | 17.9405 | 4319.0000 | 324.0000 | -0.2623 | 18.0560 | 0.2340 | -0.3330 | no | no |
| phospho | FR6 | Rep4 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR6_Rep4/` | 5373 | 19.8674 | 4319.0000 | 324.0000 | 2.1942 | 18.0560 | 0.2340 | 5.2224 | yes | yes |
| proteome | FR1 | Rep1 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR1_Rep1/` | 26214 | 22.9766 | 23237.0000 | 654.5000 | 3.0679 | 23.0663 | 0.0357 | -1.6924 | yes | no |
| proteome | FR1 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR1_Rep2/` | 25209 | 23.0913 | 23237.0000 | 654.5000 | 2.0322 | 23.0663 | 0.0357 | 0.4717 | yes | no |
| proteome | FR1 | Rep3 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR1_Rep3/` | 22484 | 23.0867 | 23237.0000 | 654.5000 | -0.7760 | 23.0663 | 0.0357 | 0.3849 | no | no |
| proteome | FR1 | Rep4 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR1_Rep4/` | 24401 | 23.0873 | 23237.0000 | 654.5000 | 1.1996 | 23.0663 | 0.0357 | 0.3962 | no | no |
| proteome | FR2 | Rep1 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR2_Rep1/` | 36702 | 23.1255 | 33255.5000 | 1714.5000 | 1.3559 | 22.9138 | 0.1166 | 1.2244 | no | no |
| proteome | FR2 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR2_Rep2/` | 34323 | 23.1644 | 33255.5000 | 1714.5000 | 0.4200 | 22.9138 | 0.1166 | 1.4493 | no | no |
| proteome | FR2 | Rep3 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR2_Rep3/` | 32371 | 22.9254 | 33255.5000 | 1714.5000 | -0.3480 | 22.9138 | 0.1166 | 0.0674 | no | no |
| proteome | FR2 | Rep4 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR2_Rep4/` | 36205 | 22.5531 | 33255.5000 | 1714.5000 | 1.1603 | 22.9138 | 0.1166 | -2.0853 | no | yes |
| proteome | FR3 | Rep1 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR3_Rep1/` | 28973 | 22.2016 | 23118.0000 | 1763.0000 | 2.2400 | 21.9967 | 0.0805 | 1.7164 | yes | no |
| proteome | FR3 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR3_Rep2/` | 29118 | 22.2407 | 23118.0000 | 1763.0000 | 2.2955 | 21.9967 | 0.0805 | 2.0440 | yes | yes |
| proteome | FR3 | Rep3 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR3_Rep3/` | 22992 | 22.0532 | 23118.0000 | 1763.0000 | -0.0482 | 21.9967 | 0.0805 | 0.4730 | no | no |
| proteome | FR3 | Rep4 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR3_Rep4/` | 25268 | 22.6767 | 23118.0000 | 1763.0000 | 0.8225 | 21.9967 | 0.0805 | 5.6971 | no | yes |
| proteome | FR4 | Rep1 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR4_Rep1/` | 31083 | 22.7523 | 29762.0000 | 737.0000 | 1.2090 | 22.7379 | 0.0602 | 0.1612 | no | no |
| proteome | FR4 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR4_Rep2/` | 34488 | 23.0124 | 29762.0000 | 737.0000 | 4.3252 | 22.7379 | 0.0602 | 3.0730 | yes | yes |
| proteome | FR4 | Rep3 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR4_Rep3/` | 30239 | 22.6961 | 29762.0000 | 737.0000 | 0.4365 | 22.7379 | 0.0602 | -0.4679 | no | no |
| proteome | FR4 | Rep4 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR4_Rep4/` | 34688 | 22.8910 | 29762.0000 | 737.0000 | 4.5082 | 22.7379 | 0.0602 | 1.7139 | yes | no |
| proteome | FR5 | Rep1 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR5_Rep1/` | 33815 | 22.3637 | 34965.5000 | 1316.0000 | -0.5897 | 22.4800 | 0.0785 | -0.9993 | no | no |
| proteome | FR5 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR5_Rep2/` | 38516 | 22.7336 | 34965.5000 | 1316.0000 | 1.8197 | 22.4800 | 0.0785 | 2.1790 | no | yes |
| proteome | FR5 | Rep3 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR5_Rep3/` | 36064 | 22.5628 | 34965.5000 | 1316.0000 | 0.5630 | 22.4800 | 0.0785 | 0.7114 | no | no |
| proteome | FR5 | Rep4 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR5_Rep4/` | 39993 | 22.8495 | 34965.5000 | 1316.0000 | 2.5767 | 22.4800 | 0.0785 | 3.1748 | yes | yes |
| proteome | FR6 | Rep1 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR6_Rep1/` | 28766 | 21.9134 | 28825.5000 | 963.5000 | -0.0417 | 21.9562 | 0.0946 | -0.3046 | no | no |
| proteome | FR6 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR6_Rep2/` | 32435 | 22.2527 | 28825.5000 | 963.5000 | 2.5268 | 21.9562 | 0.0946 | 2.1133 | yes | yes |
| proteome | FR6 | Rep3 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR6_Rep3/` | 29201 | 21.9457 | 28825.5000 | 963.5000 | 0.2629 | 21.9562 | 0.0946 | -0.0745 | no | no |
| proteome | FR6 | Rep4 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR6_Rep4/` | 31497 | 22.7412 | 28825.5000 | 963.5000 | 1.8702 | 21.9562 | 0.0946 | 5.5944 | no | yes |

## ④ |z| > 2 清单（25 条 (run, 指标)）

| # | layer | fraction | rep | run | metric | value | z |
|---|---|---|---|---|---|---|---|
| 1 | phospho | FR1 | Rep1 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR1_Rep1/` | median_log2 | 19.2291 | -2.1366 |
| 2 | phospho | FR1 | Rep4 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR1_Rep4/` | median_log2 | 18.8867 | -2.7287 |
| 3 | phospho | FR2 | Rep4 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR2_Rep4/` | median_log2 | 18.4218 | -5.1533 |
| 4 | phospho | FR4 | Rep2 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR4_Rep2/` | median_log2 | 17.6090 | -3.4321 |
| 5 | phospho | FR4 | Rep4 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR4_Rep4/` | median_log2 | 19.2360 | 2.0181 |
| 6 | phospho | FR5 | Rep1 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR5_Rep1/` | median_log2 | 18.1787 | -2.5679 |
| 7 | phospho | FR6 | Rep2 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR6_Rep2/` | n_precursors | 2777 | -3.2101 |
| 8 | phospho | FR6 | Rep4 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR6_Rep4/` | n_precursors | 5373 | 2.1942 |
| 9 | phospho | FR6 | Rep4 | `phospho/20200214_EXPL2_Evo1_AMV_SA_21m_DIA_Phos_HeLa_SUBCELL_EGF_CTRL_FR6_Rep4/` | median_log2 | 19.8674 | 5.2224 |
| 10 | proteome | FR1 | Rep1 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR1_Rep1/` | n_precursors | 26214 | 3.0679 |
| 11 | proteome | FR1 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR1_Rep2/` | n_precursors | 25209 | 2.0322 |
| 12 | proteome | FR2 | Rep4 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR2_Rep4/` | median_log2 | 22.5531 | -2.0853 |
| 13 | proteome | FR3 | Rep1 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR3_Rep1/` | n_precursors | 28973 | 2.2400 |
| 14 | proteome | FR3 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR3_Rep2/` | n_precursors | 29118 | 2.2955 |
| 15 | proteome | FR3 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR3_Rep2/` | median_log2 | 22.2407 | 2.0440 |
| 16 | proteome | FR3 | Rep4 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR3_Rep4/` | median_log2 | 22.6767 | 5.6971 |
| 17 | proteome | FR4 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR4_Rep2/` | n_precursors | 34488 | 4.3252 |
| 18 | proteome | FR4 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR4_Rep2/` | median_log2 | 23.0124 | 3.0730 |
| 19 | proteome | FR4 | Rep4 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR4_Rep4/` | n_precursors | 34688 | 4.5082 |
| 20 | proteome | FR5 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR5_Rep2/` | median_log2 | 22.7336 | 2.1790 |
| 21 | proteome | FR5 | Rep4 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR5_Rep4/` | n_precursors | 39993 | 2.5767 |
| 22 | proteome | FR5 | Rep4 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR5_Rep4/` | median_log2 | 22.8495 | 3.1748 |
| 23 | proteome | FR6 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR6_Rep2/` | n_precursors | 32435 | 2.5268 |
| 24 | proteome | FR6 | Rep2 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR6_Rep2/` | median_log2 | 22.2527 | 2.1133 |
| 25 | proteome | FR6 | Rep4 | `proteome/20200217_EXPL2_Evo1_AMV_SA_21m_DIA_Prot_HeLa_SUBCELL_EGF_CTRL_FR6_Rep4/` | median_log2 | 22.7412 | 5.5944 |

## ⑤ 计数

各 layer 的 CTRL run 中 |z| > 2 的 run 数（n = `n_precursors`，log2 = `median_log2`）：

| layer | n_ctrl_runs | n_flag_n | n_flag_log2 | n_flag_either | n_flag_both | n_entries |
|---|---|---|---|---|---|---|
| phospho | 24 | 2 | 7 | 8 | 1 | 9 |
| proteome | 24 | 8 | 8 | 12 | 4 | 16 |
| 合计 | 48 | 10 | 15 | 20 | 5 | 25 |

列说明：`n_flag_n` = n 上 |z|>2 的 CTRL run 数；`n_flag_log2` = log2 上 |z|>2 的 CTRL run 数；`n_flag_either` = 任一指标；`n_flag_both` = 两指标同时；`n_entries` = (run, 指标) 条目数。

按 layer × fraction 的条目数：

| layer | fraction | n_entries_n | n_entries_log2 |
|---|---|---|---|
| phospho | FR1 | 0 | 2 |
| phospho | FR2 | 0 | 1 |
| phospho | FR3 | 0 | 0 |
| phospho | FR4 | 0 | 2 |
| phospho | FR5 | 0 | 1 |
| phospho | FR6 | 2 | 1 |
| proteome | FR1 | 2 | 0 |
| proteome | FR2 | 0 | 1 |
| proteome | FR3 | 2 | 2 |
| proteome | FR4 | 2 | 1 |
| proteome | FR5 | 1 | 2 |
| proteome | FR6 | 1 | 2 |

## ⑥ 跳过项 / 记录

- MAD = 0 的 (cell, 指标)：0 个；因此 z 留空的格：0 个。
- 解析失败行：0。
- 输入缺失值（空串）：0 个 cell 有缺失。
- |z| 与阈值 2 相差 < 0.0005（4 位取整可能影响判定）的格：0 个。
- 方案「24 个 CTRL run / 24 行」与输入不符（实为 48），未按 24 截取；无其他跳过项（本任务只依赖 `run_medians.tsv`，已存在）。

## 自检

| 检查 | 结果 | 值 |
|---|---|---|
| 输入 240 行 | PASS | 240 |
| 解析失败 0 行 | PASS | 0 |
| 03_ctrl_qc.tsv 24 行（方案原文） | FAIL | 48 |
| phospho / proteome CTRL 各 12（方案原文） | FAIL | [24, 24] |
| 03_ctrl_qc.tsv 行数 = 输入中 CTRL run 数 = 2x6x4 | PASS | 48 / 48 |
| 每个 layer x fraction 的 CTRL run = 4（Rep1-4） | PASS |  |
| 03_ctrl_qc_cells.tsv 12 行 | PASS | 12 |
| 每个 cell n_runs = 20 | PASS | [20] |
| 每个 cell 5 tp x 4 rep 各 1 | PASS |  |
| 48 个 CTRL run 的 n_precursors / median_log2 与输入原值一致 | PASS |  |
| |z|>2 条目共 25 | PASS | 25 |
| phospho 9 条 / proteome 16 条 | PASS | 9 / 16 |
| phospho FR2 Rep4 z_median_log2 ≈ -5.153 | PASS | -5.1533 |
| phospho FR6 Rep4 z_median_log2 ≈ 5.222 | PASS | 5.2224 |
| proteome FR3 Rep4 z_median_log2 ≈ 5.697 | PASS | 5.6971 |
| proteome FR6 Rep4 z_median_log2 ≈ 5.594 | PASS | 5.5944 |
| MAD=0 的 (cell, 指标) | PASS | 0 |
