# 01 统一蛋白分母：旧/新分类比较

## ① 输入

- `data/pilot/occupancy_variants_long.tsv`：306035 行（不含表头）；取 `site_id, timepoint, fraction, rep, phos_intensity, denom_residue, denom_protein`。`phos_intensity` 空 0 行；`denom_residue` 空 219111 行；`denom_protein` 空 36694 行；(site_id, fraction) 39813 个，与 site_traj 集合相同。
- `data/pilot/run_medians.tsv`：240 行（phospho 120、proteome 120）；layer = `run` 路径前缀；(timepoint, fraction, rep) 用 07 `parse_design` 解析，240/240 成功。
- `data/pilot/site_traj.tsv`：39813 行、16 列；旧分类与步骤 0 校准靶。
- 读表：`pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)`，空串 = 缺失。
- 脚本：`analysis/07_pilot2/01_build_v2.py`；输出：`01_run_factors.tsv`、`01_site_traj_v2.tsv`、`01_class_compare.tsv`、`01_class_compare.md`（均在 `analysis/07_pilot2/`）。

## ② 定义

### 07 occupancy 规则（07 `:238-306`，逐行核对后照做）
1. run 因子（log2 加性形式 `a_run`）：`num = phos_intensity × 2^a_phospho_run`，`den = denom × 2^a_proteome_run`。
2. 轨迹 = (site_id, fraction)；观测 = 长表该轨迹下全部 (timepoint, rep) 行（`phos_intensity` 均非空）。
3. 旧档位：`n_res` = `denom_residue` 非空的观测数；`n_res >= 0.5 * len(obs)` → residue 档（d = den_res），否则 protein 档（d = den_prot）；d 缺失的观测不进 occ（07 `if d:`）。
4. `inten[tp] += log2(num)`；`occ[tp] += log2(num / d)`。
5. responder：`len(CTRL) < 2` → 非响应；`base = median(CTRL)`；按 2min→8min→20min→90min，跳过 `len < 2` 的 tp，取第一个满足 `abs(median − base) >= 1.0` 且所有重复 `(x − base) * fc > 0` 的 tp；`peak_tp = o_tp or i_tp`（occupancy 达标时取其 tp，否则取 intensity 的 tp）。
6. class：o 与 i 皆响应 → both；仅 o → occupancy_only；仅 i → intensity_only；否则 none。`log2occ_tp` / `log2int_tp` = 该 tp 值数 ≥ 2 时 `round(median, 3)`（写出格式同 07 csv.writer，即 `str(round(x, 3))`），否则空。
- `median()` 用 07 `:117-118` 原函数（奇数取中间值，偶数取两中间值均值；与 `statistics.median` 同义）。

### 步骤 0 校准方法
- 初值 `a0_run = anchor_layer − median_log2_run`，`anchor_layer = statistics.median(该 layer 120 个 median_log2)`：phospho anchor = 18.84585，proteome anchor = 22.67100。
- phospho：取 site_traj `log2int_tp` 非空且长表该 (site, fraction, tp) 恰 2 个重复 {i, j} 的组，`diff = log2int_tp(07) − median(log2(phos_i) + a0_i, log2(phos_j) + a0_j)`，模型 `diff = (δ_i + δ_j)/2`；每个 (fraction, tp) cell 4 个未知数（Rep1–Rep4），`numpy.linalg.lstsq` 最小二乘；`a = a0 + δ`。
- proteome：按规则 3 定每条观测的分母，取 site_traj `log2occ_tp` 非空且恰 2 个重复有分母的组，`diff_occ = log2occ_tp(07) − median(log2(phos_i) + a_i^phos − log2(den_i) − a0_i^prot, …)`，模型 `diff_occ = −(δ_i^prot + δ_j^prot)/2`，同法解（`a^phos` 用已校准值）。
- 每个 cell 的设计矩阵秩：rank < 4 的 cell 数 = 0（60 个 cell 全部 rank = 4），因此未使用 3 重复组补充，也没有 run 保持 a0。
- `01_run_factors.tsv` 的 `n_equations`、`max_residual` 为该 run 所在 cell 的值（同 cell 4 个 run 相同）；`max_residual = max|diff − 模型值|`。a0/delta/a 写 6 位小数，脚本内部用全精度。
- 与 site_traj 比对：log2 值比较 `round(本脚本值, 3)` 与 07 表值，差以 0.001 为单位（整数千分位）计。

### 步骤 1（新）
- 规则 2、4、5、6 不变；规则 3 改为所有观测 `d = denom_protein × 2^a_prot`（a 为校准后值），`denom_protein` 空的观测不进 occ；`denom_tier` 固定写 `protein`。`protein` 列 = site_id 去掉 `_残基位置` 后缀（与 site_traj 逐行相同）。

### 步骤 2（比较）
- 旧 = `data/pilot/site_traj.tsv`；旧复现 = 步骤 0 校准后重算；新 = `01_site_traj_v2.tsv`。
- 位点级 occupancy 响应 = 该 site_id 任一 fraction 行 class ∈ {both, occupancy_only}。

## ③ 校准结果

### 与 site_traj 的一致性
| 因子 | class 一致 | denom_tier 一致 | peak_tp 一致 | 10 个 log2 列有无全一致的行 | 有无不一致单元格 | 比较的非空单元格 | 最大差 | 差 = 0.001 的单元格 | 差 > 0.001 的单元格 |
|---|---|---|---|---|---|---|---|---|---|
| a0（run_medians） | 38458/39813 | 39813/39813 | 38541/39813 | 39813/39813 | 0 | 150846 | 0.223 | 1283 | 148939 |
| a0 + δ^phos（proteome 仍 a0） | 39775/39813 | 39813/39813 | 39763/39813 | 39813/39813 | 0 | 150846 | 0.020 | 10647 | 56205 |
| a = a0 + δ（phospho + proteome） | 39813/39813 | 39813/39813 | 39812/39813 | 39813/39813 | 0 | 150846 | 0.001 | 3041 | 0 |

校准后 class/denom_tier/peak_tp 任一不一致的行：1 行。
| site_id | fraction | class 07 | class 本脚本 | tier 07 | tier 本脚本 | peak_tp 07 | peak_tp 本脚本 | 本脚本 occ fc 2min | fc 8min | fc 20min | fc 90min |
|---|---|---|---|---|---|---|---|---|---|---|---|
| O60232_S103 | FR1 | both | both | residue | residue | 8min | 2min | 1.000023 | 1.095450 | -0.373760 | -0.751468 |

（fc = median(occ_tp) − median(occ_CTRL)，用校准后 a，未取整；阈值 FC_THRESH = 1.0。）

### δ 范围
| layer | δ min | δ max | cell 数 | 每 cell 方程数 min | 每 cell 方程数 max | 每 run 方程数 min | max_residual 最大 |
|---|---|---|---|---|---|---|---|
| phospho | -0.222410 | 0.180227 | 30 | 339 | 888 | 30 | 0.000562 |
| proteome | -0.014775 | 0.019467 | 30 | 298 | 809 | 50 | 0.000574 |

### 每个 cell 的方程数、秩与最大残差
| fraction | timepoint | phospho n_eq | rank | max_res | δ Rep1/2/3/4 | proteome n_eq | rank | max_res | δ Rep1/2/3/4 |
|---|---|---|---|---|---|---|---|---|---|
| FR1 | CTRL | 699 | 4 | 0.000516 | -0.16246 / 0.09065 / 0.15709 / 0.11857 | 638 | 4 | 0.000517 | 0.00025 / 0.00602 / 0.00742 / 0.01782 |
| FR1 | 2min | 625 | 4 | 0.000533 | 0.17839 / 0.17703 / 0.13772 / 0.14664 | 574 | 4 | 0.000574 | 0.00081 / 0.01775 / -0.00346 / -0.00210 |
| FR1 | 8min | 772 | 4 | 0.000527 | 0.07953 / 0.14602 / 0.18023 / 0.16762 | 642 | 4 | 0.000506 | -0.00834 / 0.00013 / -0.00938 / 0.00255 |
| FR1 | 20min | 750 | 4 | 0.000521 | 0.15105 / 0.07405 / -0.08514 / 0.09376 | 675 | 4 | 0.000520 | 0.00484 / 0.01018 / -0.00011 / 0.00345 |
| FR1 | 90min | 749 | 4 | 0.000530 | 0.15624 / 0.13429 / 0.06291 / 0.12982 | 683 | 4 | 0.000526 | -0.00168 / 0.01289 / 0.00222 / 0.01458 |
| FR2 | CTRL | 779 | 4 | 0.000529 | -0.10541 / 0.06359 / 0.09596 / 0.07824 | 762 | 4 | 0.000543 | -0.01417 / -0.01179 / 0.00573 / -0.00462 |
| FR2 | 2min | 746 | 4 | 0.000518 | 0.03981 / 0.03873 / 0.11549 / 0.07672 | 714 | 4 | 0.000514 | -0.00301 / -0.00160 / -0.00715 / 0.00037 |
| FR2 | 8min | 781 | 4 | 0.000531 | -0.04214 / 0.10295 / 0.16538 / 0.03149 | 755 | 4 | 0.000524 | 0.00531 / 0.00338 / 0.00094 / -0.00724 |
| FR2 | 20min | 705 | 4 | 0.000533 | 0.07129 / 0.13088 / 0.04214 / 0.03325 | 676 | 4 | 0.000512 | -0.00831 / -0.00664 / -0.00510 / -0.00691 |
| FR2 | 90min | 772 | 4 | 0.000516 | 0.05554 / 0.08749 / 0.13100 / 0.14365 | 759 | 4 | 0.000516 | -0.01399 / -0.01054 / -0.00400 / -0.00499 |
| FR3 | CTRL | 514 | 4 | 0.000533 | -0.22241 / -0.01444 / 0.05239 / -0.01423 | 427 | 4 | 0.000526 | -0.00329 / 0.00917 / 0.00477 / 0.00787 |
| FR3 | 2min | 339 | 4 | 0.000562 | -0.08734 / -0.10378 / -0.06235 / 0.01459 | 298 | 4 | 0.000509 | 0.00677 / 0.00506 / -0.00375 / 0.00361 |
| FR3 | 8min | 390 | 4 | 0.000539 | -0.13642 / -0.00292 / -0.03107 / 0.00111 | 328 | 4 | 0.000538 | 0.00036 / 0.01342 / -0.00105 / 0.00593 |
| FR3 | 20min | 434 | 4 | 0.000538 | -0.04022 / -0.08255 / -0.06380 / 0.03708 | 358 | 4 | 0.000557 | -0.00359 / 0.01947 / 0.00088 / 0.00527 |
| FR3 | 90min | 505 | 4 | 0.000506 | -0.02229 / 0.03415 / -0.00872 / 0.08046 | 439 | 4 | 0.000540 | -0.00158 / 0.01177 / -0.00246 / 0.00964 |
| FR4 | CTRL | 888 | 4 | 0.000522 | -0.12084 / -0.10961 / -0.04157 / 0.07566 | 809 | 4 | 0.000520 | 0.00279 / -0.00221 / -0.00189 / -0.00151 |
| FR4 | 2min | 612 | 4 | 0.000511 | -0.07541 / 0.05218 / 0.05390 / -0.01161 | 551 | 4 | 0.000515 | -0.00183 / -0.00125 / -0.00612 / 0.01934 |
| FR4 | 8min | 589 | 4 | 0.000547 | -0.13175 / 0.01754 / 0.01447 / -0.05297 | 521 | 4 | 0.000517 | 0.00064 / 0.00421 / 0.00304 / 0.00148 |
| FR4 | 20min | 565 | 4 | 0.000532 | -0.02379 / -0.05606 / -0.06354 / -0.04944 | 507 | 4 | 0.000531 | -0.00159 / 0.00845 / 0.00012 / 0.00707 |
| FR4 | 90min | 565 | 4 | 0.000518 | 0.02111 / 0.09249 / 0.05363 / 0.02258 | 542 | 4 | 0.000518 | -0.00980 / 0.00153 / -0.00068 / 0.01561 |
| FR5 | CTRL | 679 | 4 | 0.000541 | -0.14560 / -0.00281 / 0.04153 / 0.05712 | 689 | 4 | 0.000511 | -0.01247 / -0.00236 / 0.00062 / -0.00448 |
| FR5 | 2min | 595 | 4 | 0.000540 | 0.03868 / 0.11959 / 0.12641 / 0.04766 | 588 | 4 | 0.000530 | -0.00111 / -0.00343 / -0.00841 / -0.01067 |
| FR5 | 8min | 644 | 4 | 0.000541 | 0.00620 / 0.09826 / 0.14428 / 0.10349 | 622 | 4 | 0.000539 | 0.01176 / -0.00631 / -0.00927 / -0.01478 |
| FR5 | 20min | 617 | 4 | 0.000504 | 0.07122 / 0.09653 / 0.04569 / 0.05748 | 603 | 4 | 0.000515 | -0.00425 / 0.00169 / -0.00189 / -0.00050 |
| FR5 | 90min | 600 | 4 | 0.000540 | 0.03923 / 0.08967 / 0.11343 / 0.10521 | 603 | 4 | 0.000510 | -0.00034 / -0.00120 / -0.00688 / -0.00013 |
| FR6 | CTRL | 619 | 4 | 0.000513 | -0.20422 / -0.08508 / -0.05375 / 0.08123 | 593 | 4 | 0.000530 | 0.00766 / -0.00269 / 0.01015 / 0.00153 |
| FR6 | 2min | 498 | 4 | 0.000517 | -0.14706 / -0.07711 / -0.08745 / -0.05070 | 494 | 4 | 0.000536 | 0.01159 / -0.01441 / 0.00380 / 0.01036 |
| FR6 | 8min | 457 | 4 | 0.000524 | -0.05755 / -0.05016 / -0.06804 / -0.12120 | 463 | 4 | 0.000524 | 0.00726 / -0.00526 / 0.00916 / 0.01049 |
| FR6 | 20min | 559 | 4 | 0.000541 | -0.01819 / 0.00024 / -0.02482 / -0.08079 | 529 | 4 | 0.000537 | 0.00097 / -0.00063 / -0.00026 / 0.01534 |
| FR6 | 90min | 570 | 4 | 0.000534 | -0.05751 / -0.00199 / 0.09837 / 0.05827 | 567 | 4 | 0.000543 | 0.00935 / 0.00445 / 0.00250 / 0.00981 |

## ④ 比较表

### C1 四分类计数
| class | n_old (site_traj) | n_old_replicated (步骤0, 校准后 a) | n_new |
|---|---|---|---|
| both | 3202 | 3202 | 3291 |
| occupancy_only | 2447 | 2447 | 2128 |
| intensity_only | 2072 | 2072 | 1983 |
| none | 32092 | 32092 | 32411 |
| 合计 | 39813 | 39813 | 39813 |

### C2 交叉表（行 = 旧 site_traj class，列 = 新 class）
| 旧 \ 新 | both | occupancy_only | intensity_only | none | 行合计 |
|---|---|---|---|---|---|
| both | 3064 | 0 | 138 | 0 | 3202 |
| occupancy_only | 0 | 1919 | 0 | 528 | 2447 |
| intensity_only | 227 | 0 | 1845 | 0 | 2072 |
| none | 0 | 209 | 0 | 31883 | 32092 |
| 列合计 | 3291 | 2128 | 1983 | 32411 | 39813 |

`n_changed_rows`（旧 class ≠ 新 class）= 1102

### C3 位点级 occupancy 响应（任一 fraction 为 both/occupancy_only）
| n_old | n_new | 两者皆 | 仅旧 | 仅新 |
|---|---|---|---|---|
| 3787 | 3711 | 3496 | 291 | 215 |

### C4 按 fraction 的 both + occupancy_only 行数
| fraction | 旧 (site_traj) | 新 |
|---|---|---|
| FR1 | 1130 | 1105 |
| FR2 | 1135 | 1053 |
| FR3 | 652 | 612 |
| FR4 | 848 | 828 |
| FR5 | 1132 | 1091 |
| FR6 | 752 | 730 |
| all | 5649 | 5419 |

### 补充计数
| 项目 | 口径 | 值 |
|---|---|---|
| rows_log2occ_CTRL_nonempty | old_site_traj | 13919 |
| rows_log2occ_CTRL_nonempty | old_replicated_step0 | 13919 |
| rows_log2occ_CTRL_nonempty | new | 14094 |
| obs_total | value | 306035 |
| obs_occ_missing_due_to_denom_protein_missing | new | 36694 |
| obs_entering_occ | old_replicated_step0 | 263503 |
| obs_entering_occ | new | 269341 |
| log2int_cells_nonempty_either | new_vs_site_traj | 80684 |
| log2int_cells_string_diff | new_vs_site_traj | 1653 |
| log2int_cells_presence_diff | new_vs_site_traj | 0 |
| log2int_cells_numeric_diff | new_vs_site_traj | 1653 |
| log2int_max_abs_diff | new_vs_site_traj | 0.001 |
| log2int_cells_string_diff | new_vs_old_replicated_step0 | 0 |
| class_count_both | new_with_a0_uncalibrated | 3255 |
| class_count_occupancy_only | new_with_a0_uncalibrated | 2094 |
| class_count_intensity_only | new_with_a0_uncalibrated | 2054 |
| class_count_none | new_with_a0_uncalibrated | 32410 |
| rows_class_diff | new_a0_vs_new_calibrated | 1337 |

### 自检
| 检查 | 结果 |
|---|---|
| v2 行数 39813、16 列 | 通过 |
| v2 (site_id, fraction) 集合与顺序同 site_traj | 通过 |
| v2 denom_tier 全为 protein | 通过 |
| log2int_* 与 site_traj 逐格相同（字符串） | 未通过 |
| log2int_* 与 步骤0 校准复现 逐格相同（字符串） | 通过 |
| 校准后 class 39813/39813 | 通过 |
| 校准后 denom_tier 39813/39813 | 通过 |
| 校准后 peak_tp 39813/39813 | 未通过 |
| 校准后 log2 有无一致（单元格） | 通过 |
| 校准后 非空 log2 差 ≤ 0.001 | 通过 |
| 所有 cell 残差最大值 ≤ 0.0015 | 通过 |
| C1 三列合计各 39813 | 通过 |
| C2 行边际 = n_old、列边际 = n_new | 通过 |
| C3 两者皆+仅旧 = n_old、两者皆+仅新 = n_new | 通过 |
| C4 六个 fraction 之和 = all | 通过 |
| 新 log2occ_CTRL 非空行数 ≥ 旧 | 通过 |

## ⑤ 跳过/缺失项

- 3 重复组补充：60 个 cell 均 rank = 4，未触发，未使用。
- 保持 a0 的 run：0 个。
- `data/pilot/sites.tsv`、`summary*.txt`：缺，本任务不依赖，跳过。
