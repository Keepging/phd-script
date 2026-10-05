# 06_pilot 设计方案（advisor 写，subagent 执行）

## 0. 通用约定（三个任务共用）

**输入**（只读，不改）：`data/pilot/`
| 文件 | 行数(不含表头) | 列 |
|---|---|---|
| `site_traj.tsv` | 39813 | `site_id, protein, fraction, denom_tier, class, peak_tp, log2occ_CTRL, log2occ_2min, log2occ_8min, log2occ_20min, log2occ_90min, log2int_CTRL, log2int_2min, log2int_8min, log2int_20min, log2int_90min`。一行 = 一个 位点×fraction 轨迹。`log2occ_*`/`log2int_*` 为该时间点 ≥2 个重复时的中位数，否则为空串（07 `:294-295`）。`class` ∈ {both, occupancy_only, intensity_only, none}（07 `:288-290`）。`peak_tp` = 07 responder 返回的**第一个**达标时间点，不是峰值（07 `:278-284`）。`denom_tier` ∈ {residue, protein}（07 `:263`）。 |
| `site_profiles.tsv` | 26740 | `site_id, protein, proteome_tier, timepoint, phos_shift_vs_ctrl, prot_shift_vs_ctrl, site_vs_prot_divergence, divergence_delta_vs_ctrl`（本轮任务不用）。 |
| `run_medians.tsv` | 240 | `run, n_precursors, median_log2`。**没有单独的 layer 列**，layer 是 `run` 的路径前缀 `phospho/` 或 `proteome/`。 |
| `site_features.tsv` | 18268 | `site_id, accession, residue, position, single_{backbone,sidechain,disoMine,earlyFolding,helix,sheet,coil,ppII}, win5_{同上}`。**无序列列**。site_id 与 site_traj 完全重合（18268/18268）。**特征列缺失值是字面字符串 `NA`**（single/win5 各列 35 行，disoMine 46 行；响应位点里 backbone 为 NA 的 5 个），必须当缺失处理。 |
| `kinase_name_mapping.csv` | 247 | `our_name, atlas_name, kin_type, status`。MAPK1→ERK2、MAPK3→ERK1。 |
| `kinase_pssm/` | 不存在 | — |

**脚本里已有的阈值，不改**：`FC_THRESH = 1.0`、`MIN_REPS = 2`（07 `:22-23`）、`PCT = 90`（08 `:45`）。`EGF_TPS = [2min, 8min, 20min, 90min]`，`FRS = FR1..FR6`。

**读表**：`pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)`，空串 = 缺失；需要数值时显式 `float()`。不要让 pandas 把 "none"/"NA" 当成 NaN。

**输出**：全部写到 `analysis/06_pilot/`。每个任务同时保存执行脚本 `0N_*.py`（可复跑）。md 文件固定结构：① 输入（路径、行数）② 定义（照抄本方案）③ 结果表 ④ 跳过/缺失项。只报计数和数值，不解释生物学，不给建议。缺输入写"缺 X，跳过"，不造数据。subagent 不做 git 提交。

**基线数（advisor 已核，subagent 结果必须对上）**：site_traj 39813 行、18268 位点、4996 蛋白；class 计数 none 32092 / both 3202 / occupancy_only 2447 / intensity_only 2072；denom_tier protein 29687 / residue 10126；occupancy 响应行（both+occupancy_only）5649 行、3787 位点、1837 蛋白，其中 `log2occ_90min` 为空的 801 行、`log2occ_CTRL` 为空的 0 行。

---

## 任务 1：漏斗表 → `01_funnel.md` + `01_funnel.tsv` + `01_funnel.py`

### 输入
`data/pilot/site_traj.tsv`（全部列）。蛋白级丰度：三张表都**没有**蛋白级强度列（site_traj 只有位点级 `log2int_*`），该列固定填 `无`；另加一列位点级对照值（见下）。

### 行级条件（在一个 位点×fraction 行上判断）
- `c2`（②）：`log2int_CTRL != ''` 且 `log2int_{tp} != ''` 对某个 tp ∈ EGF_TPS 成立。（非空 ⇔ 该时间点 ≥2 重复，07 `:295`。）
- `c3`（③"有蛋白档分母"）：`log2occ_CTRL != ''` 且 `log2occ_{tp} != ''` 对某个 tp ∈ EGF_TPS 成立。解释：occupancy 值存在 ⇔ 07 为该轨迹找到了 proteome 分母（residue 档或 protein 档均算，07 `:270-272`；residue 档存在时 protein 档必然存在，07 `:230-234`）。按 `denom_tier` 拆开报 protein / residue 两个子计数，residue 子计数即"残基档分母覆盖"对照行。
- 对照行（宽口径）：任一 `log2occ_*` 非空（不要求 CTRL+EGF 配对）。

### 位点级 / 蛋白级条件
- 记 `S3` = 在 ≥1 个 fraction 上满足 `c2 & c3` 的位点集合（③ 的位点级口径）。
- `c4`（④）：位点满足 `c2 & c3` 的 fraction 数 ≥ 2。
- `c5`（⑤）：位点所在蛋白在 `S3` 中有 ≥2 个位点。
- 位点"occupancy 响应" := 该位点任一行 `class ∈ {both, occupancy_only}`；"非响应" := 位点 ∈ S3 且无此类行。
- `c6`（⑥）：位点所在蛋白在 `S3` 中既有 ≥1 个响应位点又有 ≥1 个非响应位点。

### 漏斗（累积：每步在上一步剩余集合上再加条件）
| step | 剩余行 | 剩余位点 |
|---|---|---|
| ① 全部 | 全表 | 全表 |
| ② | 行满足 c2 | 有 ≥1 行满足 c2 的位点 |
| ③ | 行满足 c2 & c3 | S3 |
| ④ | 行满足 c2&c3 且位点满足 c4 | S3 ∩ c4 |
| ⑤ | 行满足 c2&c3 且位点 ∈ (④剩余 ∩ c5) | ④剩余 ∩ c5 |
| ⑥ | 行满足 c2&c3 且位点 ∈ (⑤剩余 ∩ c6) | ⑤剩余 ∩ c6 |
每步报：`n_rows`、`n_sites`（distinct site_id）、`n_proteins`（distinct protein）、`protein_abundance_median` = `无`、`site_log2int_CTRL_median` = 剩余行中非空 `log2int_CTRL` 的中位数（标注"位点级磷酸肽强度，非蛋白丰度"）。
**按 fraction 再报一次**：scope = FR1..FR6，行 = 该 fraction 内满足该步行级条件且位点属于该步剩余集合的行；n_sites/n_proteins 在这些行上数。
**对照行**：(a) ③ 的 protein 档子计数；(b) ③ 的 residue 档子计数；(c) 宽口径"任一 log2occ 非空"的行/位点/蛋白数；(d) 全表中 10 个数值列全空的行数。

### 输出列（`01_funnel.tsv`，制表符）
`step, condition, scope, n_rows, n_sites, n_proteins, protein_abundance_median, site_log2int_CTRL_median`
step ∈ {1..6, ctrl_a, ctrl_b, ctrl_c, ctrl_d}；scope ∈ {all, FR1..FR6}。`01_funnel.md` 内渲染同一张表并附定义。

### 验收标准
1. step1/all：n_rows 39813、n_sites 18268、n_proteins 4996。
2. step2/all n_rows = 14963；step3/all n_rows = 13106，ctrl_a = 9037、ctrl_b = 4069（两者之和 = 13106）；ctrl_d = 16854。
3. all 口径下 n_rows、n_sites、n_proteins 随 step 1→6 单调不增。
4. 每一步 6 个 fraction 的 n_rows 之和 = 该步 all 的 n_rows（行级口径下必须相等）。
5. tsv 与 md 数字一致；无 NaN；每步每个 scope 各恰好一行（10 steps × 7 scopes = 70 行）。
6. 脚本可复跑，输出确定。

---

## 任务 2：归一化趋势 → `02_run_medians.md` + `02_run_medians.tsv` + `02_run_medians_summary.tsv` + `02_run_medians.png`（另存 `02_run_medians_phospho.png`、`02_run_medians_proteome.png`）+ `02_run_medians.py`

### 输入
`data/pilot/run_medians.tsv`（`run, n_precursors, median_log2`）。

### 解析（照 07 `parse_design`，`:32-38`）
- `layer` = `run` 在第一个 `/` 之前的部分（phospho | proteome）。
- `timepoint` = 正则 `_(2min|8min|20min|90min|CTRL)_`（忽略大小写）；`fraction` = `_(FR\d)_`；`rep` = `_(Rep\d)`。解析失败的行计数并列出（预期 0）。
- 预期 2 layer × 5 tp × 6 FR × 4 rep = 240，每个 cell 恰 4 行。

### 偏移定义（两种都算，都只是计数）
- `anchor_{layer}` = 该 layer 120 个 `median_log2` 的中位数（07 `factors()` `:239-242` 的 anchor）。`offset_anchor` = `median_log2 − anchor_layer`。
- `offset_ctrl` = `median_log2 − median_log2(同 layer、同 fraction、同 rep 的 CTRL run)`；CTRL 行本身留空。
- "方向一致" := 同一 layer×fraction×timepoint 的 4 个重复的 offset 符号全同（全 >0 或全 <0；出现 0 记为不一致，并单独计数）。

### 输出
- `02_run_medians.tsv`（240 行）：`layer, run, timepoint, fraction, rep, n_precursors, median_log2, anchor_layer, offset_anchor, offset_ctrl`
- `02_run_medians_summary.tsv`（60 行 = 2×6×5）：`layer, fraction, timepoint, n_reps, mean_median_log2, sd_median_log2, min_median_log2, max_median_log2, mean_n_precursors, n_pos_anchor, n_neg_anchor, n_zero_anchor, consistent_anchor, n_pos_ctrl, n_neg_ctrl, n_zero_ctrl, consistent_ctrl`。sd 用样本标准差（ddof=1）。CTRL 行的 `*_ctrl` 列留空。
- `02_run_medians.md`：① 输入 ② 定义 ③ anchor 值（两个 layer）④ summary 表全文（60 行）⑤ 一致性计数：每个 layer×fraction 中 consistent_anchor=yes 的时间点数（/5）、consistent_ctrl=yes 的时间点数（/4）；总计 /60 和 /48 ⑥ 解析失败列表（应为空）。
- 图：`02_run_medians.png` 2 行（上 phospho、下 proteome）× 6 列（FR1..FR6）；每子图横轴 `CTRL,2min,8min,20min,90min`（按此顺序等距），纵轴 `median_log2`，4 条线 = Rep1..Rep4，图例一次；同一 layer 的 6 个子图共享 y 范围；标题含 layer 和 fraction；dpi ≥150。另存两张单 layer 图（1×6）。画图前先调用 `dataviz` skill 取配色与规范。

### 验收标准
1. 解析失败 0 行；240 行全部有 layer/timepoint/fraction/rep；60 个 cell 每个 `n_reps = 4`。
2. advisor 抽查 2 个 cell：手算 mean/sd 与 summary 一致（容差 1e-6）。
3. `offset_ctrl` 在 CTRL 行为空，在 EGF 行 = 该行 median_log2 − 对应 CTRL 行值（抽查 2 行）。
4. 三张 png 存在、非空、能打开；子图数正确。
5. md 中一致性计数与 summary.tsv 的 consistent 列逐 cell 复核一致。

---

## 任务 3：时间规律小实验 → `03_time_pilots.md` + `03_time_pilots_rows.tsv` + `03_time_pilots_sites.tsv` + `03_time_pilots.py`

### 输入
`data/pilot/site_traj.tsv`；`data/pilot/site_features.tsv`（存在）；`data/pilot/kinase_name_mapping.csv`；kinase-library wheel 已下载在 `/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad/kl/kinase_library-1.8.0-py3-none-any.whl`。

### 基集
`R` = site_traj 中 `class ∈ {both, occupancy_only}` 的行（预期 5649 行、3787 位点、1837 蛋白）。

### 3a 形状
每行：`fc_{tp} = float(log2occ_{tp}) − float(log2occ_CTRL)`，tp ∈ EGF_TPS 且该列非空。
- `peak_tp_max` = |fc| 最大的 tp（并列取最早）。同时保留 07 的 `peak_tp`（第一个达标 tp）。
- `fc_90` = `fc_90min`；`type_90`：`log2occ_90min` 为空 → `90min缺失`；`|fc_90| < FC_THRESH(1.0)` → `transient`（暂时型）；否则 `sustained`（持续型）。
- 报：表 A1 = `peak_tp_max` × `type_90` 计数（4×3 + 合计）；表 A2 = 07 `peak_tp` × `type_90`；表 A3 = `peak_tp_max` × 07 `peak_tp` 交叉表。
- 位点级类型（一个位点可有多个响应 fraction）：忽略 `90min缺失` 行后，全 transient → `transient`，全 sustained → `sustained`，两者都有 → `mixed`，无可判行 → `untyped`。报四类位点数。
- backbone 对照（site_features 存在）：用 `single_backbone`（另用 `win5_backbone` 重复一遍）。找同时有 ≥1 个 transient 位点和 ≥1 个 sustained 位点的蛋白（mixed/untyped 位点不参与）；每个蛋白算 `d = median(backbone of transient sites) − median(backbone of sustained sites)`。报：蛋白数、d 的中位数、d>0 / d<0 / d=0 的蛋白数、缺 backbone 的位点数（`NA` 字面值视为缺失；响应位点里预期 5）。不做检验。

### 3b 空间特异性
每个位点：`n_frac_resp` = 其 `class ∈ {both, occupancy_only}` 的 fraction 数；`n_frac_evaluable` = 其 `log2occ_CTRL` 非空且任一 EGF `log2occ` 非空的 fraction 数（含非响应行，在全表上数）。
- 报 `n_frac_resp` = 1 / 2 / 3+ 的位点数（3+ 再拆 3/4/5/6），合计 3787。
- 报 `n_frac_resp × n_frac_evaluable` 交叉表（行 1..6，列 1..6）。
- 区室合并：FR1,FR2→Cyt；FR3,FR4→Mem；FR5,FR6→Nuc。`n_comp_resp` = 响应 fraction 覆盖的不同区室数；报 1/2/3 的位点数，合计 3787；再报 `n_comp_resp × n_comp_evaluable` 交叉表。

### 3c ERK 阳性对照
- 序列来源检查：列出 site_features.tsv 的全部列名、site_traj/site_profiles 列名，确认无序列/窗口列；`find` 仓库内 `*.fasta`、`biophys_json`（预期均无）。无序列 ⇒ kinase-library 打分和"+1 位是 P"代理**都无法执行**，写"缺序列，跳过"。
- 仍需记录：`pip install <上面的 wheel 路径>`（限时 10 分钟）是否成功、`import kinase_library` 是否成功、`kl.get_kinase_list(kin_type='ser_thr')` 中是否含 `ERK1`、`ERK2`（名称映射用 kinase_name_mapping.csv 的 MAPK1→ERK2、MAPK3→ERK1，调用方式参考 `code/Protein contour/Phospho/martinez_network_check/kl_core.py:7-8,35-42`）。装不上也只记录，不重试超过一次。
- 报基线：全部响应行中 `peak_tp_max ∈ {2min, 8min}` 的比例（行口径和位点口径各一个，位点口径用位点的任一响应行），ERK 候选比例一栏写"缺序列，跳过"。

### 输出
- `03_time_pilots_rows.tsv`（5649 行）：`site_id, protein, fraction, class, denom_tier, peak_tp_07, peak_tp_max, fc_2min, fc_8min, fc_20min, fc_90min, fc_at_peak, type_90`
- `03_time_pilots_sites.tsv`（3787 行）：`site_id, protein, n_frac_resp, n_frac_evaluable, n_comp_resp, n_comp_evaluable, site_type_90, single_backbone, win5_backbone, peak_early_any`（peak_early_any = 任一响应行 peak_tp_max ∈ {2min,8min}）
- `03_time_pilots.md`：① 输入 ② 定义 ③ 3a 表 A1/A2/A3、位点类型计数、backbone 对照（single 与 win5 两段）④ 3b 四张表 ⑤ 3c 序列检查结果、pip/import/kinase 列表检查结果、基线比例 ⑥ 跳过项。

### 验收标准
1. R = 5649 行 / 3787 位点 / 1837 蛋白；A1 合计 5649；`90min缺失` = 801；A1 与 A2 行合计相同。
2. 位点类型四类之和 = 3787；3b 的 1/2/3+ 之和 = 3787；区室 1/2/3 之和 = 3787；交叉表行列边际与上述一致。
3. `peak_tp_max = 90min` 的行里 `transient` 必须为 0（峰在 90 且 |fc_90| ≥ 该行所有 |fc| ≥ 1）。
4. backbone 对照：d>0 + d<0 + d=0 = 蛋白数；缺 backbone 位点数 = 5（single 与 win5 相同；advisor 修正，原稿写 0 是因为未识别 `NA` 字面值）。
5. 3c 明确写出"缺序列，跳过"，且 pip/import 结果有记录。
6. rows.tsv 5649 行、sites.tsv 3787 行，无 NaN（允许空串表示缺失）。
