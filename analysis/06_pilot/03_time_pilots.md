# 03 时间规律小实验（任务 3）

## ① 输入

- `data/pilot/site_traj.tsv`：39813 行（不含表头），18268 个 site_id，4996 个 protein。class：none 32092，both 3202，occupancy_only 2447，intensity_only 2072；denom_tier：protein 29687，residue 10126。
- `data/pilot/site_features.tsv`：18268 行，18268 个 site_id，20 列；site_traj 位点在其中的覆盖：18268/18268。`single_backbone` 为字面 "NA" 的行 35，`win5_backbone` 为 "NA" 的行 35（全表）。
- `data/pilot/site_profiles.tsv`：26740 行（仅在 3c 列出列名）。
- `data/pilot/kinase_name_mapping.csv`：247 行（逗号分隔，`pd.read_csv(..., sep=',', dtype=str, keep_default_na=False)`）。
- kinase-library wheel：`/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad/kl/kinase_library-1.8.0-py3-none-any.whl`（存在）。
- 读表：`pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)`，空串 = 缺失，数值处显式 `float()`。
- 基集 R（class ∈ {both, occupancy_only}）：5649 行、3787 位点、1837 蛋白；其中 `log2occ_90min` 为空 801 行，`log2occ_CTRL` 为空 0 行。
- 阈值（脚本已有，未改）：FC_THRESH = 1.0，MIN_REPS = 2，PCT = 90（任务 3 无用到 PCT 之处）。
- 脚本：`analysis/06_pilot/03_time_pilots.py`；输出：`analysis/06_pilot/03_time_pilots_rows.tsv`、`analysis/06_pilot/03_time_pilots_sites.tsv`、`analysis/06_pilot/03_time_pilots.md`。

## ② 定义（照抄 00_design.md 任务 3 的 基集/3a/3b/3c，原文如下）

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

**执行补充（协调方更正，现已并入方案）**：site_features.tsv 的特征列缺失值为字面字符串 `"NA"`（非空串）。`"NA"` 视为缺失：该位点不参与 backbone 中位数计算，sites.tsv 中该位点的 `single_backbone`/`win5_backbone` 留空串；“缺 backbone 的位点数”预期为 5（方案初稿写 0）。

**实现细节**：`peak_tp_max` 按 EGF_TPS 顺序（2min→8min→20min→90min）遍历，仅当 |fc| 严格大于当前最大值时替换，即并列取最早。rows.tsv 中 fc 列以 3 位小数写出（输入 `log2occ_*` 均 ≤3 位小数）；缺该时间点则为空串。`peak_tp_07` = site_traj 的 `peak_tp` 原值。`peak_early_any` 写作 `True`/`False`。区室：FR1/FR2→Cyt，FR3/FR4→Mem，FR5/FR6→Nuc。backbone 对照中，蛋白的入选（≥1 transient 且 ≥1 sustained）按剔除 backbone 缺失位点后的位点判断；中位数用 `statistics.median`。

## ③ 3a 形状

### 表 A1：peak_tp_max × type_90（行数）

| peak_tp_max \ type_90 | transient | sustained | 90min缺失 | 合计 |
|---|---|---|---|---|
| 2min | 704 | 255 | 316 | 1275 |
| 8min | 768 | 428 | 261 | 1457 |
| 20min | 513 | 298 | 224 | 1035 |
| 90min | 0 | 1882 | 0 | 1882 |
| 合计 | 1985 | 2863 | 801 | 5649 |

### 表 A2：07 peak_tp × type_90（行数）

| peak_tp_07 \ type_90 | transient | sustained | 90min缺失 | 合计 |
|---|---|---|---|---|
| 2min | 914 | 737 | 411 | 2062 |
| 8min | 675 | 596 | 228 | 1499 |
| 20min | 396 | 329 | 162 | 887 |
| 90min | 0 | 1201 | 0 | 1201 |
| 合计 | 1985 | 2863 | 801 | 5649 |

### 表 A3：peak_tp_max × 07 peak_tp（行数）

| peak_tp_max \ peak_tp_07 | 2min | 8min | 20min | 90min | 合计 |
|---|---|---|---|---|---|
| 2min | 1184 | 47 | 24 | 20 | 1275 |
| 8min | 394 | 1009 | 31 | 23 | 1457 |
| 20min | 193 | 186 | 628 | 28 | 1035 |
| 90min | 291 | 257 | 204 | 1130 | 1882 |
| 合计 | 2062 | 1499 | 887 | 1201 | 5649 |

- 对角线（peak_tp_max = peak_tp_07）行数：3951；非对角线：1698。
- 每行 max|fc| < FC_THRESH 的响应行数：0。
- 精度核对（同一定义改用十进制精确减法重算）：peak_tp_max 不一致 0 行，type_90 不一致 0 行；十进制下 |fc_90| 恰为 1.000 的行 1（按 `<1.0` 规则归 sustained）；max|fc| 出现并列的行：浮点 5，十进制 6（两种算法下 peak_tp_max 结果均一致）。

### 位点级类型（忽略 90min缺失 行）

| site_type_90 | 位点数 |
|---|---|
| transient | 1114 |
| sustained | 1683 |
| mixed | 538 |
| untyped | 452 |
| 合计 | 3787 |

### backbone 对照：`single_backbone`

| 项 | 值 |
|---|---|
| 同时有 ≥1 transient 与 ≥1 sustained 位点（均有 backbone 值）的蛋白数 | 347 |
| 参与计算的 transient 位点数 / sustained 位点数 | 546 / 684 |
| d 的中位数 | 0.0185 |
| d > 0 蛋白数 | 198 |
| d < 0 蛋白数 | 148 |
| d = 0 蛋白数（精确相等） | 1 |
| abs(d) < 1e-9 蛋白数（容差核对） | 1 |
| 缺 backbone 的位点数（R 的 3787 位点中；字面 "NA"） | 5（"NA" 5） |
| 缺 backbone 位点的 site_type_90 分布 | transient 3，sustained 1，mixed 0，untyped 1 |
| 缺 backbone 位点 site_id | P08195_S33，P78312_S933，Q14669_S1066，Q14669_S354，Q8WX92_S605 |
| 若不剔除缺失位点、按位点类型判断入选的蛋白数 | 348 |
| 因剔除缺失位点而落选的蛋白 | Q14669 |

### backbone 对照：`win5_backbone`

| 项 | 值 |
|---|---|
| 同时有 ≥1 transient 与 ≥1 sustained 位点（均有 backbone 值）的蛋白数 | 347 |
| 参与计算的 transient 位点数 / sustained 位点数 | 546 / 684 |
| d 的中位数 | 0.00872727 |
| d > 0 蛋白数 | 200 |
| d < 0 蛋白数 | 147 |
| d = 0 蛋白数（精确相等） | 0 |
| abs(d) < 1e-9 蛋白数（容差核对） | 0 |
| 缺 backbone 的位点数（R 的 3787 位点中；字面 "NA"） | 5（"NA" 5） |
| 缺 backbone 位点的 site_type_90 分布 | transient 3，sustained 1，mixed 0，untyped 1 |
| 缺 backbone 位点 site_id | P08195_S33，P78312_S933，Q14669_S1066，Q14669_S354，Q8WX92_S605 |
| 若不剔除缺失位点、按位点类型判断入选的蛋白数 | 348 |
| 因剔除缺失位点而落选的蛋白 | Q14669 |

## ④ 3b 空间特异性

- 响应行中不满足 evaluable 条件的行数：0。

### 表 B1：n_frac_resp 分布（位点数）

| n_frac_resp | 位点数 |
|---|---|
| 1 | 2547 |
| 2 | 814 |
| 3+ | 426 |
|   其中 3 | 283 |
|   其中 4 | 99 |
|   其中 5 | 35 |
|   其中 6 | 9 |
| 合计 | 3787 |

### 表 B2：n_frac_resp × n_frac_evaluable（位点数）

| n_frac_resp \ n_frac_evaluable | 1 | 2 | 3 | 4 | 5 | 6 | 合计 |
|---|---|---|---|---|---|---|---|
| 1 | 1406 | 625 | 260 | 117 | 77 | 62 | 2547 |
| 2 | 0 | 331 | 227 | 114 | 75 | 67 | 814 |
| 3 | 0 | 0 | 91 | 87 | 59 | 46 | 283 |
| 4 | 0 | 0 | 0 | 32 | 27 | 40 | 99 |
| 5 | 0 | 0 | 0 | 0 | 14 | 21 | 35 |
| 6 | 0 | 0 | 0 | 0 | 0 | 9 | 9 |
| 合计 | 1406 | 956 | 578 | 350 | 252 | 245 | 3787 |

### 表 B3：n_comp_resp 分布（位点数；FR1-2 Cyt，FR3-4 Mem，FR5-6 Nuc）

| n_comp_resp | 位点数 |
|---|---|
| 1 | 2874 |
| 2 | 738 |
| 3 | 175 |
| 合计 | 3787 |

### 表 B4：n_comp_resp × n_comp_evaluable（位点数）

| n_comp_resp \ n_comp_evaluable | 1 | 2 | 3 | 合计 |
|---|---|---|---|---|
| 1 | 1972 | 612 | 290 | 2874 |
| 2 | 0 | 441 | 297 | 738 |
| 3 | 0 | 0 | 175 | 175 |
| 合计 | 1972 | 1053 | 762 | 3787 |

## ⑤ 3c ERK 阳性对照

### 序列来源检查

- `data/pilot/site_features.tsv` 列名（20）：`site_id`, `accession`, `residue`, `position`, `single_backbone`, `single_sidechain`, `single_disoMine`, `single_earlyFolding`, `single_helix`, `single_sheet`, `single_coil`, `single_ppII`, `win5_backbone`, `win5_sidechain`, `win5_disoMine`, `win5_earlyFolding`, `win5_helix`, `win5_sheet`, `win5_coil`, `win5_ppII`
- `data/pilot/site_traj.tsv` 列名（16）：`site_id`, `protein`, `fraction`, `denom_tier`, `class`, `peak_tp`, `log2occ_CTRL`, `log2occ_2min`, `log2occ_8min`, `log2occ_20min`, `log2occ_90min`, `log2int_CTRL`, `log2int_2min`, `log2int_8min`, `log2int_20min`, `log2int_90min`
- `data/pilot/site_profiles.tsv` 列名（8）：`site_id`, `protein`, `proteome_tier`, `timepoint`, `phos_shift_vs_ctrl`, `prot_shift_vs_ctrl`, `site_vs_prot_divergence`, `divergence_delta_vs_ctrl`
- 列名匹配 `seq|window|flank|motif|peptide|fasta|aa_`（不分大小写）：`data/pilot/site_features.tsv` 0 个；`data/pilot/site_traj.tsv` 0 个；`data/pilot/site_profiles.tsv` 0 个
- 仓库内（排除 .git）`find -name '*.fasta'`：0 个
- 仓库内（排除 .git）名称含 `biophys_json` 的文件/目录：0 个
- 结论：无序列 / 窗口来源。site_id（如 `P08195_S33`）只含残基与位置，不含侧翼序列。
- kinase-library 打分：**缺序列，跳过**。
- “+1 位是 P” 代理：**缺序列，跳过**。

### pip / import / kinase 列表检查

- `python3 -m pip install <wheel>` 第 1 次：成功，returncode=0，耗时 37.9 s，输出末行：`      Successfully uninstalled matplotlib-3.11.2` / `Successfully installed adjustText-1.3.0 alabaster-1.0.0 babel-2.18.0 beautifulsoup4-4.15.0 biopython-1.88 bleach-6.4.0 cloudpickle-3.1.2 defusedxml-0.7.1 docutils-0.21.2 et-xmlfile-2.0.0 fastjsonschem` / `WARNING: Running pip as the 'root' user can result in broken permissions and conflicting behaviour with the system package manager. It is recommended to use a virtual environment instead: https://pip.`
- `import kinase_library`：成功（version 1.8.0）；耗时 2.3 s。
- `kl.get_kinase_list(kin_type='ser_thr')`：311 个激酶；含 `ERK1`：True；含 `ERK2`：True。
- 调用方式参考：`code/Protein contour/Phospho/martinez_network_check/kl_core.py:7-8,35-42`。
- pip 前后环境版本（`importlib.metadata`）：pandas 3.0.6 → 2.2.3；numpy 2.4.6 → 1.26.4；matplotlib 3.11.2 → 3.8.4。本脚本进程内的 pandas 为 pip 前已加载的 3.0.6。
- kinase_name_mapping.csv 中 ERK 映射：MAPK1→ERK2（ser_thr，alias）；MAPK3→ERK1（ser_thr，alias）

### 基线比例

| 口径 | peak_tp_max ∈ {2min, 8min} | 分母 | 比例 | ERK 候选比例 |
|---|---|---|---|---|
| 行（全部响应行） | 2732 | 5649 | 0.4836 | 缺序列，跳过 |
| 位点（任一响应行） | 2179 | 3787 | 0.5754 | 缺序列，跳过 |

## ⑥ 跳过项

- 3c kinase-library 打分（ERK1/ERK2）：缺序列，跳过。
- 3c “+1 位是 P” 代理：缺序列，跳过。
- 3c ERK 候选比例（行口径、位点口径）：缺序列，跳过。
- 未做任何统计检验（按方案）。
- backbone 对照中缺 backbone（字面 "NA"）的位点：single 5 个、win5 5 个，不参与中位数计算。

### 附：验收自检（脚本自动核对）

| # | 项 | 实际 | 期望 | 通过 |
|---|---|---|---|---|
| 1 | R 行数 | 5649 | 5649 | 是 |
| 1 | R 位点数 | 3787 | 3787 | 是 |
| 1 | R 蛋白数 | 1837 | 1837 | 是 |
| 1 | A1 合计 | 5649 | 5649 | 是 |
| 1 | A1 中 90min缺失 列合计 | 801 | 801 | 是 |
| 1 | R 中 log2occ_90min 为空行数 | 801 | 801 | 是 |
| 1 | A1 合计 = A2 合计 | True | True | 是 |
| 1 | A1 与 A2 的 type_90 列合计逐列相同 | True | True | 是 |
| 1 | A3 合计 | 5649 | 5649 | 是 |
| 2 | 位点类型四类之和 | 3787 | 3787 | 是 |
| 2 | n_frac_resp 1/2/3+ 之和 | 3787 | 3787 | 是 |
| 2 | n_comp_resp 1/2/3 之和 | 3787 | 3787 | 是 |
| 2 | n_frac 交叉表行边际 = n_frac_resp 分布 | [2547, 814, 283, 99, 35, 9] | [2547, 814, 283, 99, 35, 9] | 是 |
| 2 | n_frac 交叉表列边际 = n_frac_evaluable 分布 | [1406, 956, 578, 350, 252, 245] | [1406, 956, 578, 350, 252, 245] | 是 |
| 2 | n_comp 交叉表行边际 = n_comp_resp 分布 | [2874, 738, 175] | [2874, 738, 175] | 是 |
| 2 | n_comp 交叉表列边际 = n_comp_evaluable 分布 | [1972, 1053, 762] | [1972, 1053, 762] | 是 |
| 2 | 交叉表合计（frac / comp） | (3787, 3787) | (3787, 3787) | 是 |
| 3 | peak_tp_max=90min 且 transient 的行数 | 0 | 0 | 是 |
| 4 | single_backbone: d>0 + d<0 + d=0 = 蛋白数 | 347 | 347 | 是 |
| 4 | single_backbone: 缺 backbone 位点数（方案验收 4：5；初稿写 0） | 5 | 5 | 是 |
| 4 | win5_backbone: d>0 + d<0 + d=0 = 蛋白数 | 347 | 347 | 是 |
| 4 | win5_backbone: 缺 backbone 位点数（方案验收 4：5；初稿写 0） | 5 | 5 | 是 |
| 5 | 3c 写明“缺序列，跳过”且有 pip/import 记录 | True | True | 是 |
| 6 | rows.tsv 数据行数 | 5649 | 5649 | 是 |
| 6 | sites.tsv 数据行数 | 3787 | 3787 | 是 |
| 6 | rows.tsv / sites.tsv 中 nan/NaN/None/NA 文本个数 | (0, 0) | (0, 0) | 是 |
| 6 | rows.tsv / sites.tsv 列名与方案一致 | (True, True) | (True, True) | 是 |
