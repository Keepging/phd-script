# 09 验收记录（advisor，pilot2）

方案：`00_design.md`。执行：四个 Opus subagent（任务 1 与 3 并行；任务 2、4 在任务 1 通过后并行）。验收方法同上一轮：advisor 先用独立实现算参考数（scratchpad，不入库），交回后逐格比对，再按方案验收标准逐条核。四个任务都是一轮通过，没有重派。

## 输入与分支说明
- 新输入 `occupancy_long.tsv`、`occupancy_variants_long.tsv`、`human_2026-08-19_UP2026_02.fasta` 原本只在分支 `claude/determined-goldberg-atmb8n`（commit 25bbda6），已合并进当前分支（merge commit c90097c），`data/` 未改。
- `sites.tsv`、`summary*.txt` 不在 `data/pilot/`：缺，跳过（本轮任务不依赖）。
- 关键事实：**`occupancy_long.tsv` 的 `prot_intensity` 是残基级分母**（05 `:250-252`），蛋白级分母在 `occupancy_variants_long.tsv` 的 `denom_protein`（06 `:222`）。任务 1 用的是后者。

## 任务 1 统一蛋白分母 — 合格
文件：`01_site_traj_v2.tsv`（39813 行 × 16 列，列、行序同 site_traj）、`01_class_compare.tsv/.md`、`01_run_factors.tsv`（240 行）、`01_build_v2.py`。
验收：
- 校准：用 `run_medians.tsv` 推出的因子只能复现 site_traj 的 class 38458/39813（advisor 预先发现，见方案 0 节）；按方案从"恰 2 个重复"的组最小二乘解出每 run 常数 δ（phospho δ ∈ [−0.222, 0.180]，proteome δ ∈ [−0.015, 0.019]，60 个 cell 全满秩，残差 ≤ 0.00057），校准后旧分类复现：class 39813/39813、denom_tier 39813/39813、peak_tp 39812/39813、全部非空 log2 值差 ≤ 0.001。唯一 peak_tp 不一致行 O60232_S103 FR1：2min 的 occupancy fc = 1.000023，离阈值 1.0 差 2.3e-5，class 两边都是 both。
- advisor 用 subagent 的因子独立重算：旧分类与 site_traj、新分类与 `01_site_traj_v2.tsv` 均逐行一致（class / tier / peak_tp / 10 列 log2 值）。
- 新旧四分类（C1）：both 3202→3291，occupancy_only 2447→2128，intensity_only 2072→1983，none 32092→32411。换类行 1102（C2 只在 both↔intensity_only、occupancy_only↔none 之间移动，强度列不变所致）。位点级 occupancy 响应 3787→3711，两者皆 3496，仅旧 291，仅新 215。有 occ 值的行 13919→14094。
方案两处未完全达成、advisor 接受：(1) "log2int 逐格相同"有 1653 格差 0.001，是 07 表 3 位小数取整与 δ 的 1e-5 精度共同造成；(2) peak_tp 差 1 行（上述边界值）。两者都不是规则差异。

## 任务 2 漏斗重跑 + 丰度 — 合格
文件：`02_funnel_v2.tsv`（91 行）、`02_funnel_v2.md`、`02_protein_abundance.tsv`（3356 个蛋白）、`02_funnel_v2.py`。
验收：70 个有参考值的格全部与 advisor 一致（含丰度中位数，容差 2e-4）；版本 B 的 step1–3 与版本 A 逐格相同；A step2 n_rows = 14963 与旧漏斗相同；两版单调不增；每步 6 个 fraction 之和 = all；抽查 3 个蛋白丰度一致。
all 口径 n_sites：版本 A ① 18268 → ② 6861 → ③ 6502 → ④ 3173 → ⑤ 2657 → ⑥ 1977；版本 B ③ 6502 → ⑤ 5228 → ⑥ 3894。旧漏斗（残基/蛋白混用）对应 ③ 6439 / ④ 3129 / ⑤ 2616 / ⑥ 1905。
丰度列：`protein_abundance_median`（蛋白所在 CTRL run 的 log2 `denom_protein` 中位数的中位数）各步 25.26–25.59；有丰度值的蛋白 3356/4996，因为长表只在该蛋白有位点被观测的 run 里才有 `denom_protein`。

## 任务 3 CTRL 质检 — 合格
文件：`03_ctrl_qc.tsv`（48 行）、`03_ctrl_qc_cells.tsv`（12 行）、`03_ctrl_qc.md`、`03_ctrl_qc.py`。
验收：|z|>2 的 (run, 指标) 共 25 条，与 advisor 逐条一致（容差 0.001）；flag 列与 z 一致；48 个 CTRL run 原值与输入一致；MAD=0 的 cell 0 个。
**任务原文"24 个 CTRL run（磷酸 12 + 蛋白 12）"与数据不符**：数据是 2 layer × 6 fraction × 4 rep = 48 个。subagent 按 48 报，advisor 已修正方案。
|z|>2 的 CTRL run：phospho 24 个里 8 个（n 上 2、log2 上 7），proteome 24 个里 12 个（n 上 8、log2 上 8）。最大 |z|：proteome FR3 Rep4 median_log2 5.697、proteome FR6 Rep4 5.594、phospho FR6 Rep4 5.222、phospho FR2 Rep4 −5.153。

## 任务 4 时间规律重跑 + ERK — 合格
文件：`04_time_rules_v2.md`、`04_time_rules_v2_rows.tsv`（5419 行）、`04_erk_sites.tsv`（18268 行 × 17 列）、`04_time_rules_v2.py`、`04_erk_score.py`。
验收：
- (a) 旧列与 06_pilot 完全相同；新列（R = 5419 行 / 3711 位点 / 1809 蛋白）A1、位点类型（1153/1628/483/447）、n_frac（2566/756/258/96/27/8）、n_comp（2854/712/145）与 advisor 参考一致。
- (b) venv 安装、导入 kinase-library 1.8.0 成功，批量接口 `PhosphoProteomics(...).predict(kin_type='ser_thr')`，耗时 50 s；18268 个窗口与 advisor 独立取的逐个相同，中心残基不一致 0；+1P 7148 与 advisor 一致；17676 个 S/T 位点打分成功，592 个 Y 位点库不打 ser_thr 分（分数空）。ERK 候选（ERK1 或 ERK2 percentile ≥ 95）1698 个（ERK1 1362、ERK2 1252、两者 916）。advisor 在 venv 内用逐条接口抽查 3 个位点：score 一致，percentile 差 ≤ 0.01。E1 全部 n 从 tsv 复算一致。
- 全局 pandas/numpy/matplotlib 前后均为 2.2.3 / 1.26.4 / 3.8.4，未改。
E1（位点口径，早峰 = 任一响应行 peak_tp_max ∈ {2min, 8min}）：

| | 旧 n | 旧早峰 | 旧比例 | 新 n | 新早峰 | 新比例 |
|---|---|---|---|---|---|---|
| 全部响应位点 | 3787 | 2179 | 0.5754 | 3711 | 2133 | 0.5748 |
| ERK 候选（定义 1） | 606 | 376 | 0.6205 | 588 | 369 | 0.6276 |
| +1P（定义 2） | 2300 | 1348 | 0.5861 | 2273 | 1325 | 0.5829 |
| 集合前 5%（对照） | 336 | 220 | 0.6548 | 325 | 216 | 0.6646 |

行口径见 `04_time_rules_v2.md` ⑤。

## 改过什么
1. 方案任务 3 的 CTRL run 数 24 → 48（数据事实）。
2. 方案任务 4 "前 5%"的基数：subagent 按"全部打分位点"（17676 个 S/T）取，方案原文写 18268；advisor 接受 subagent 的口径（Y 位点无 ser_thr 分数，无法排序），记录在此。
3. 无重派；四个任务输出未被 advisor 手改。

## 没解决 / 需要知道的问题
1. **`run_medians.tsv` 不是 07 内部的 run 中位数**：与 07 实际使用的因子相差每 run 一个常数（phospho 最大 0.22 log2）。本轮靠从 site_traj 反解 δ 才复现 07；若以后要脱离 site_traj 重算，需要 07 原始的 run 中位数（或重跑 07 导出）。
2. **蛋白档分母覆盖**：36694 个观测（12.0%）没有 `denom_protein`，这些观测在新规则下不进 occupancy；旧规则下进 occ 的观测 263503，新规则 269341。
3. **峰值口径**：07 的 `peak_tp` 是第一个达标时间点；本轮 3a/E1 的早峰用 |fc| 最大的 `peak_tp_max`。
4. **ERK 定义的百分位**是 kinase-library 自带的、相对其背景集的百分位，不是本数据集内的百分位；集合内前 5% 另列为对照。592 个 Y 位点不在 ERK 候选判断范围内（E1 的"全部响应位点"行含这些 Y 位点：旧 49、新 50）。
5. **丰度**只覆盖 3356/4996 个蛋白；按 fraction 报时丰度不分 fraction（蛋白级同一值）。
6. 任务 1 的 1 行 peak_tp 差异和 1653 格 0.001 差异（取整/精度）如上。
