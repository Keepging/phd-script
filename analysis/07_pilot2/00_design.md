# 07_pilot2 设计方案（advisor 写，subagent 执行）

沿用 `analysis/06_pilot/00_design.md` 的通用约定：读表 `pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)`，空串 = 缺失，`site_features.tsv` 的字面值 `NA` 也是缺失；中位数用 `statistics.median`（07 `median()` 同义），保留 4 位小数（与 07 表比对时按 07 的 3 位）；run 名解析用 07 `parse_design` 的三个正则；阈值 `FC_THRESH = 1.0`、`MIN_REPS = 2`；`EGF_TPS = [2min, 8min, 20min, 90min]`、`TPS = [CTRL] + EGF_TPS`、`FRS = FR1..FR6`。只读 `data/`、`code/`；输出全部写 `analysis/07_pilot2/`，每个任务保存可复跑脚本 `0N_*.py`；只报计数和数值；缺输入写"缺 X，跳过"；subagent 不做 git 提交。**kinase-library 只装在独立 venv 里，不碰全局 pandas/numpy/matplotlib。**

## 0. 输入现状（advisor 已核）
| 文件 | 行数 | 列 / 说明 |
|---|---|---|
| `data/pilot/occupancy_long.tsv` | 306035 | `site_id, timepoint, fraction, rep, phos_intensity, prot_intensity, occupancy_index`（05 产物）。**`prot_intensity` 是残基级分母**（05 `:250-252`：覆盖该残基的肽段强度和），不是蛋白级。空 219111 行。 |
| `data/pilot/occupancy_variants_long.tsv` | 306035 | `site_id, timepoint, fraction, rep, phos_intensity, denom_residue, denom_protein, denom_residue_repmean, occ_v1_strict, occ_vA_fallback, vA_type, occ_vB_repmean`（06 产物）。**`denom_protein` 才是蛋白级分母**（06 `:222`：该蛋白在该 run 的全部可定位肽段强度和，= 07 的 `den_prot`，07 `:230`）。`denom_residue` 与 `occupancy_long.prot_intensity` 逐行相等；`phos_intensity` 两表相等。`denom_protein` 空 36694 行；`vA_type` protein 182417 / residue 86924 / none 36694。 |
| `data/pilot/run_medians.tsv` | 240 | `run, n_precursors, median_log2`。 |
| `data/pilot/site_traj.tsv` | 39813 | 07 产物（旧分类）。列见 06_pilot 方案。 |
| `data/pilot/site_features.tsv` | 18268 | 同上一轮。 |
| `data/pilot/human_2026-08-19_UP2026_02.fasta` | 20431 条 | `>sp|ACC|NAME ...`，05 `load_fasta` 正则 `>\w+\|([^|]+)\|` 全部可解析；site_traj 的 4996 个蛋白全在。 |
| `data/pilot/sites.tsv`、`summary*.txt` | **缺** | 不在 `data/pilot/`，涉及处写"缺，跳过"（本轮任务不依赖它们）。 |
| `kinase_name_mapping.csv` | 247 | MAPK1→ERK2、MAPK3→ERK1。 |

**07 的 occupancy 规则（任务 1 必须逐条照做，07 `:238-290`）**
1. 每个 run 一个 log2 中心化因子：`anchor_layer = median(该 layer 全部 run 的 run 中位数)`，`factor_run = 2^(anchor − m_run)`；`num = phos_intensity × factor_phospho_run`，`den = denom × factor_proteome_run`。
2. 轨迹 = (site_id, fraction)，观测 = 该轨迹下所有 (timepoint, rep) 行（`phos_intensity` 存在即为观测）。
3. 旧分母档位：`n_res = 有残基分母的观测数`，`n_res ≥ 0.5 × 观测数` → residue 档（d = den_res），否则 protein 档（d = den_prot）；d 缺失的观测不进 occ。
4. `inten[tp] = log2(num)`；`occ[tp] = log2(num / d)`。
5. responder：CTRL ≥ MIN_REPS 个值，`base = median(CTRL)`；按 2min→8min→20min→90min 顺序，取第一个满足 `len ≥ MIN_REPS` 且 `|median − base| ≥ FC_THRESH` 且所有重复 `(x − base) × fc > 0` 的 tp；没有则非响应。`peak_tp` = 该第一个达标 tp。
6. class：both / occupancy_only / intensity_only / none；表里 `log2occ_tp`、`log2int_tp` = 该 tp ≥ MIN_REPS 时的中位数（3 位小数），否则空。

**advisor 已发现的事实（任务 1 必须处理）**：用 `run_medians.tsv` 的 `median_log2` 按规则 1 算因子，再按规则 2–6 重算旧分类，与 `site_traj.tsv` 比：`denom_tier` 39813/39813 一致，`log2occ_CTRL` 有无一致 39813/39813，但 `class` 只有 38458/39813 一致，`log2int_*`/`log2occ_*` 与 07 表的差最大 0.223。差值在"恰 2 个重复"的 (site, fraction, tp) 组内只取 2 个相邻值（±0.001 取整），说明 **07 实际用的每 run 因子与 run_medians 推出的因子相差一个 run 常数 δ_run**（`run_medians.tsv` 不是 07 内部的 run 中位数）。δ 可从数据精确解出（见任务 1 步骤 0）。

---

## 任务 1：统一蛋白分母 → `01_site_traj_v2.tsv` + `01_class_compare.md` + `01_class_compare.tsv` + `01_run_factors.tsv` + `01_build_v2.py`

### 输入
`occupancy_variants_long.tsv`（取 `phos_intensity`、`denom_residue`、`denom_protein`）、`run_medians.tsv`、`site_traj.tsv`（旧分类与校准靶）。

### 步骤 0：校准 07 的 run 因子（log2 加性形式 `a_run`，`num = phos × 2^a`）
- 初值：`a0_run = anchor_layer − median_log2_run`（规则 1）。
- 对 phospho run：取 site_traj 中 `log2int_tp` 非空且该 (site, fraction, tp) 在长表里**恰有 2 个重复** {i, j} 的组：`diff = log2int_tp(07) − median(log2(phos_i) + a0_i, log2(phos_j) + a0_j) = (δ_i + δ_j)/2`。同一 (fraction, tp) 的 4 个 run 构成 4 个未知数，各重复对给方程，最小二乘解 δ；报每个 cell 的方程数、残差最大值（预期 ≤ 0.0015，因为 07 表取 3 位小数）。`a_run = a0_run + δ_run`。
- 对 proteome run：先用旧档位规则（规则 3）确定每条观测用的分母，取 site_traj `log2occ_tp` 非空且**恰有 2 个重复有分母**的组，`diff_occ = log2occ_tp(07) − median(log2(phos_i/den_i) + a_i^phos − a0_i^prot, …) = −(δ_i^prot + δ_j^prot)/2`，同法解 δ^prot。
- 若某 cell 方程不足（<4 个独立对），用 3 重复组补（中位数 = 中间值，所用 rep 由排序确定），并在 md 记录；仍不足则该 run 保持 a0，记录。
- 输出 `01_run_factors.tsv`：`layer, run, timepoint, fraction, rep, median_log2, a0, delta, a, n_equations, max_residual`。
- **验收门槛**：用校准后的 a 重算旧分类（规则 2–6），与 site_traj 比：`class` 一致 39813/39813，`denom_tier` 一致 39813/39813，`peak_tp` 一致 39813/39813，所有非空 `log2occ_*`/`log2int_*` 差 ≤ 0.001，有无一致。达不到 100% 就把不一致行数、最大差值如实写进 md，仍继续后续步骤，不得为了凑数改规则。

### 步骤 1：蛋白档重算
规则 2、4、5、6 不变，规则 3 改为：**所有观测一律 d = den_prot**（`denom_protein × 2^a_prot`），缺 `denom_protein` 的观测不进 occ；`denom_tier` 列固定写 `protein`。输出 `01_site_traj_v2.tsv`，列与 site_traj 完全相同：`site_id, protein, fraction, denom_tier, class, peak_tp, log2occ_CTRL..log2occ_90min, log2int_CTRL..log2int_90min`，行数 39813，(site_id, fraction) 顺序同 site_traj。

### 步骤 2：新旧比较 → `01_class_compare.tsv` / `.md`
- 表 C1：四分类计数并列，三列：`class, n_old(site_traj), n_old_replicated(步骤0), n_new`，各列合计 39813。
- 表 C2：交叉表 旧(site_traj class) × 新，4×4 + 边际；`n_changed_rows`。
- 表 C3：位点级 occupancy 响应（任一 fraction 为 both/occupancy_only）：`n_old, n_new, 两者皆, 仅旧, 仅新`。
- 表 C4：按 fraction 的新旧 both+occupancy_only 行数。
- 另报：有 occ 值（log2occ_CTRL 非空）的行数 旧 / 新；`denom_protein` 缺失导致 occ 缺失的观测数。
- `log2int_*` 新旧必须逐格相同（分母不影响强度），报不同格数（预期 0）。

### 验收标准（advisor 参考值，用 run_medians 初值 a0 算的，校准后 old 应变为 100% 一致；new 的数会随校准略变，接受偏差但要说明）
1. v2 表 39813 行、16 列、(site_id, fraction) 集合与 site_traj 相同；`log2int_*` 与 site_traj 逐格相同（校准成功时）。
2. 校准复现 site_traj 的 class 一致率写明；目标 39813/39813。
3. C1 三列各合计 39813；C2 行边际 = n_old，列边际 = n_new；C3 "两者皆 + 仅旧" = n_old，"两者皆 + 仅新" = n_new。
4. 新分类有 occ 的行数 ≥ 旧（protein 档是 residue 档的超集，06 `:222-226`）。
5. 脚本可复跑、确定。

---

## 任务 3：CTRL 质检 → `03_ctrl_qc.tsv` + `03_ctrl_qc_cells.tsv` + `03_ctrl_qc.md` + `03_ctrl_qc.py`
（与任务 1 并行，不依赖它。）

### 输入
`run_medians.tsv`。解析同 06_pilot 任务 2。

### 定义
- 组 = 同 layer × fraction 的 20 个 run（5 tp × 4 rep）。
- 对 `n_precursors` 和 `median_log2` 各算：`med = median(20 值)`，`MAD = median(|x − med|)`，`z = (x − med) / (1.4826 × MAD)`；MAD = 0 时 z 留空并记录。
- 报全部 CTRL run。**advisor 修正**：数据里 CTRL run 是 48 个（2 layer × 6 fraction × 4 rep），任务原文写的 24（磷酸 12 + 蛋白 12）与数据不符，按 48 报。

### 输出
- `03_ctrl_qc.tsv`（48 行；原稿写 24，见上）：`layer, fraction, rep, run, n_precursors, median_log2, cell_median_n, cell_mad_n, z_n_precursors, cell_median_log2, cell_mad_log2, z_median_log2, flag_n(|z|>2 yes/no), flag_log2(|z|>2 yes/no)`。
- `03_ctrl_qc_cells.tsv`（12 行）：`layer, fraction, n_runs, median_n, mad_n, median_log2, mad_log2`。
- `03_ctrl_qc.md`：① 输入 ② 定义 ③ 24 行表 ④ |z|>2 清单（run、指标、值、z）⑤ 计数：phospho / proteome 各有几个 CTRL run 在 n、在 log2 上 |z|>2 ⑥ 跳过项。

### 验收标准
1. 48 行（原稿 24，advisor 修正）；每个 cell n_runs = 20；z 用 20 个 run（含非 CTRL）算，不是只用 CTRL。
2. advisor 参考：|z|>2 的 (run, 指标) 共 25 条，其中 phospho 9 条（含 FR2 Rep4 median_log2 z = −5.153、FR6 Rep4 median_log2 z = 5.222）、proteome 16 条（含 FR3 Rep4 median_log2 z = 5.697、FR6 Rep4 median_log2 z = 5.594）。逐条对上（容差 0.001）。
3. 48 个 CTRL run 的 n_precursors、median_log2 与 run_medians.tsv 原值一致。

---

## 任务 2：漏斗重跑 + 丰度 → `02_funnel_v2.tsv` + `02_funnel_v2.md` + `02_protein_abundance.tsv` + `02_funnel_v2.py`
（**依赖任务 1 的 `01_site_traj_v2.tsv`**，任务 1 验收通过后再派。）

### 输入
`analysis/07_pilot2/01_site_traj_v2.tsv`（新分类）、`data/pilot/occupancy_variants_long.tsv`（算蛋白丰度）、`data/pilot/site_traj.tsv`（蛋白归属）。

### 蛋白丰度
`protein_abundance = median over CTRL runs of log2(denom_protein)`：对每个蛋白，取长表中 `timepoint == CTRL` 且 `denom_protein` 非空的 (fraction, rep) 组合（最多 24 个 run），**按 (protein, fraction, rep) 去重**（同一蛋白多个位点的 `denom_protein` 是同一个数），取 `log2(denom_protein)` 的中位数；用原值，不做 run 中心化。输出 `02_protein_abundance.tsv`：`protein, n_ctrl_runs_with_denom, abundance_median_log2`。每层漏斗的 `protein_abundance_median` = 剩余位点所在蛋白（去重）的 `abundance_median_log2` 的中位数；无值的蛋白不参与，并报参与蛋白数。

### 漏斗
定义同 06_pilot 任务 1（c2、c3、c4、c5、c6、S3、响应/非响应），但全部基于 v2 表：c3 = `log2occ_CTRL` 非空且某 EGF `log2occ` 非空（v2 里分母只有 protein 档）；对照行 ctrl_a/ctrl_b 不再适用（都是 protein 档），改报：ctrl_c（任一 log2occ 非空）、ctrl_d（10 列全空）。
- 版本 A（含④）：① → ② → ③ → ④ → ⑤ → ⑥。
- 版本 B（不含④）：① → ② → ③ → ⑤ → ⑥，⑤ 在 ③ 的剩余集合上取交。
- 两版都按 fraction 再报一次（口径同 06_pilot）。

### 输出列（`02_funnel_v2.tsv`）
`version, step, condition, scope, n_rows, n_sites, n_proteins, n_proteins_with_abundance, protein_abundance_median, site_log2int_CTRL_median`。version ∈ {A, B}；A 有 step 1–6 + ctrl_c、ctrl_d；B 有 1,2,3,5,6（ctrl 行不重复）。scope ∈ {all, FR1..FR6}。md 渲染两版表并附定义；另附一张"06_pilot 旧漏斗 all 口径 vs 本轮版本 A all 口径"的并列表（旧数直接抄 `analysis/06_pilot/01_funnel.tsv`）。

### 验收标准
1. 版本 A step1/all = 39813 / 18268 / 4996；版本 B step1–3 与 A 相同；两版 step5、6 的数各自单调不增；版本 B 的 step5 n_sites ≥ 版本 A 的 step5 n_sites。
2. 每步 6 个 fraction 的 n_rows 之和 = all。
3. 丰度：`02_protein_abundance.tsv` 的蛋白数 ≤ 4996；advisor 抽查 3 个蛋白手算一致。
4. `log2int` 相关的格（step2 的 n_rows 等）与 06_pilot 旧漏斗相同（强度不随分母变）：step2/all n_rows 必须 = 14963。

---

## 任务 4：时间规律重跑 + ERK → `04_time_rules_v2.md` + `04_time_rules_v2_rows.tsv` + `04_erk_sites.tsv` + `04_time_rules_v2.py` + `04_erk_score.py`
（**(a) 依赖任务 1**；(b) 的打分不依赖，但比例计算依赖，整体在任务 1 通过后派。）

### (a) 3a / 3b 重跑，新旧并列
- 旧 = `data/pilot/site_traj.tsv`，新 = `analysis/07_pilot2/01_site_traj_v2.tsv`。各自取 R = class ∈ {both, occupancy_only}。
- 3a：同 06_pilot 任务 3（`fc_tp`、`peak_tp_max`、`type_90`、表 A1、位点级类型）；3b：`n_frac_resp` 1/2/3+、`n_comp_resp` 1/2/3。每张表旧、新两列并排；不做 backbone 对照（上一轮已做），不做检验。
- `04_time_rules_v2_rows.tsv`：新 R 的逐行表，列同 06_pilot `03_time_pilots_rows.tsv`。

### (b) ERK 阳性对照
- venv：`python3 -m venv <scratchpad>/venv_kl && <venv>/bin/pip install <scratchpad>/kl/kinase_library-1.8.0-py3-none-any.whl`（wheel 已在 `/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad/kl/`）。打分脚本 `04_erk_score.py` 用 venv 的 python 跑；其余脚本用系统 python。先 `python3 -c "import pandas;print(pandas.__version__)"` 记录全局版本，结束后再记录一次，必须不变。
- 序列窗口：按 05 `load_fasta` 读 FASTA；按 `kl_core.cut_window(acc, pos, half=7)` 的逻辑取 15-mer（两端 `_` 补齐）；位点 = site_traj 的 18268 个（site_id 解析为 acc、残基、位置）。核对窗口中心残基 == site_id 的残基，不一致的计数并排除。
- 打分：对每个窗口用 kinase-library 的 ser_thr 打分（优先批量接口 `kl.PhosphoProteomics`，否则逐条 `kl.Substrate(seq).predict(kin_type='ser_thr')`；记录实际用的接口和版本），取 `ERK1`、`ERK2` 的 score 和 percentile（kinase-library 自带的百分位）。**ERK 候选定义 1**：`max(percentile_ERK1, percentile_ERK2) ≥ 95`（即"百分位前 5%"）；同时记 ERK1、ERK2 各自 ≥95 的数。**定义 2（粗代理）**：+1 位残基 == `P`。另报集合内相对口径：按 `max(score_ERK1, score_ERK2)` 在 18268 位点中排前 5% 的数（只作对照）。
- 比例：对旧 R 和新 R 各算（位点口径：位点任一响应行 `peak_tp_max ∈ {2min, 8min}`；行口径：响应行）：全部响应位点早峰比例；ERK 候选（定义 1）中早峰比例；+1 P 位点中早峰比例。表 E1：行 = {全部响应位点, ERK 候选(定义1), +1P(定义2), 集合前5%(对照)} × 列 = {旧 n, 旧 早峰 n, 旧 比例, 新 n, 新 早峰 n, 新 比例}，位点口径和行口径各一张。
- `04_erk_sites.tsv`（18268 行）：`site_id, protein, residue, position, window15, center_ok, plus1_P, erk1_score, erk1_percentile, erk2_score, erk2_percentile, erk_candidate_p95, erk_top5pct_set, resp_old, resp_new, early_old, early_new`。

### 验收标准
1. (a) 旧列的 A1、位点类型、3b 计数与 06_pilot `03_time_pilots.md` 完全相同（704/255/316…；1114/1683/538/452；2547/814/426；2874/738/175）。新列合计 = 新 R 的行数 / 位点数（任务 1 的 C3 n_new）。
2. (b) 18268 位点全部有窗口；center_ok 不一致数写明；打分成功数写明；ERK 候选数、+1P 数写明；E1 的 n 列能和 04_erk_sites.tsv 过滤复算一致。
3. 全局 pandas/numpy/matplotlib 版本前后相同。
