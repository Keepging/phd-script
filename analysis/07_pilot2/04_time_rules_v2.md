# 04 时间规律重跑（旧/新分类并列）+ ERK 阳性对照（任务 4）

## ① 输入与环境

### 输入

| 文件 | 行数（不含表头） | 说明 |
|---|---|---|
| `data/pilot/site_traj.tsv` | 39813 | 旧分类；18268 位点、4996 蛋白 |
| `analysis/07_pilot2/01_site_traj_v2.tsv` | 39813 | 新分类（任务 1，蛋白档分母）；(site_id, fraction) 顺序与旧表相同：True |
| `data/pilot/human_2026-08-19_UP2026_02.fasta` | 20431 条序列 | 05 `load_fasta`（`05_build_occupancy.py:40-57`）读法 |
| `data/pilot/kinase_name_mapping.csv` | 247 | ERK 映射：MAPK1→ERK2（ser_thr，alias）；MAPK3→ERK1（ser_thr，alias） |
| `analysis/06_pilot/03_time_pilots_rows.tsv` | 5649 | 仅用于自检（旧逐行表须逐格相同）及 E1 行口径旧列复算 |

- 读表：`pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)`，空串 = 缺失。阈值 FC_THRESH = 1.0，MIN_REPS = 2（未改）。
- 旧 R（class ∈ {both, occupancy_only}）：5649 行、3787 位点、1837 蛋白；class：both 3202，occupancy_only 2447；`log2occ_90min` 为空 801 行。
- 新 R：5419 行、3711 位点、1809 蛋白；class：both 3291，occupancy_only 2128；`log2occ_90min` 为空 760 行。

### 环境

| 项 | 记录 |
|---|---|
| 全局版本检查命令 | `python3 -c "import pandas,numpy,matplotlib;print(pandas.__version__,numpy.__version__,matplotlib.__version__)"`（`python3` = `/usr/local/bin/python3`） |
| 开始前（任务开始、建 venv 之前） | `2.2.3 1.26.4 3.8.4`（2026-10-05T22:33:33+00:00） |
| 结束后（本脚本全部计算与 tsv 写出之后，子进程执行） | `2.2.3 1.26.4 3.8.4`（2026-10-05T22:41:34+00:00） |
| 前后是否相同 | 相同 |
| venv 路径 | `/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad/venv_kl`（`python3 -m venv` 创建；pyvenv.cfg：`home = /usr/local/bin`；`include-system-site-packages = false`；`version = 3.11.15`；`executable = /usr/bin/python3.11`；`command = /usr/local/bin/python3 -m venv /tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad/venv_kl`） |
| venv 安装命令 | `/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad/venv_kl/bin/pip install /tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad/kl/kinase_library-1.8.0-py3-none-any.whl`（联网装依赖，限时 900 s） |
| venv 安装结果 | pip_rc=0；pip_end 2026-10-05T22:34:32+00:00 elapsed=53s；venv 创建 venv_create_start 2026-10-05T22:33:39+00:00 |
| venv 内相关包 | `kinase-library-1.8.0`、`matplotlib-3.8.4`、`numpy-1.26.4`、`pandas-2.2.3` |
| kinase-library 版本 / 导入位置 | 1.8.0 / `/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad/venv_kl/lib/python3.11/site-packages/kinase_library/__init__.py` |
| 打分进程 | `/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad/venv_kl/bin/python`（sys.prefix = venv，base_prefix = `/usr`；venv 内 pandas 2.2.3） |
| 本脚本进程 | `/usr/local/bin/python3`（系统 python3，pandas 2.2.3） |
| 全局环境备注 | 全局 site-packages 中原已有 kinase-library 1.8.0（上一轮 `analysis/06_pilot/03_time_pilots.py` 用全局 pip 安装，见 06 的 md 第 186-190 行）；本轮未调用、未改动它 |
| 安装/打分日志 | `/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad/t4_venv_install.log`、`/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad/t4_erk_score_log.json` |

## ② 定义

- **R**：class ∈ {both, occupancy_only} 的行；旧 = `data/pilot/site_traj.tsv`，新 = `analysis/07_pilot2/01_site_traj_v2.tsv`。
- **3a**（沿用 06_pilot 任务 3）：`fc_tp = float(log2occ_tp) − float(log2occ_CTRL)`（tp ∈ 2min/8min/20min/90min 且非空）；`peak_tp_max` = |fc| 最大的 tp（按 2→8→20→90 遍历、严格大于才替换，即并列取最早）；`type_90`：`log2occ_90min` 空 → `90min缺失`；|fc_90| < 1.0 → `transient`；否则 `sustained`。表 A1 = peak_tp_max × type_90 行数。位点级类型：忽略 `90min缺失` 行后，全 transient → transient，全 sustained → sustained，两者都有 → mixed，无可判行 → untyped。
- **3b**：`n_frac_resp` = 位点 class ∈ {both, occupancy_only} 的 fraction 数；区室 FR1/FR2→Cyt、FR3/FR4→Mem、FR5/FR6→Nuc，`n_comp_resp` = 响应 fraction 覆盖的区室数。
- **窗口**：site_id 按 `^(.+)_([A-Z])(\d+)$` 解析为 acc、残基、位置；FASTA 按 05 `load_fasta`；15-mer 按 `kl_core.cut_window(acc, pos, half=7)`（中心 ±7，两端 `_` 补齐）；`center_ok` = 窗口第 8 位 == site_id 残基（不一致的不送打分）；`plus1_P` = 窗口第 9 位（+1 位）== `P`（C 端补齐的 `_` 记 False）。
- **打分**：kinase-library ser_thr，实际接口与参数见 ④；`erk1_*`/`erk2_*` = 库输出的 `ERK1_score`、`ERK1_percentile`、`ERK2_score`、`ERK2_percentile`（库默认取整：score 3 位、percentile 2 位，tsv 原样写出）。
- **定义 1（ERK 候选）**：`max(erk1_percentile, erk2_percentile) ≥ 95`；未打分位点该列为空串。另报 ERK1、ERK2 各自 ≥ 95 的数。
- **定义 2（+1P）**：`plus1_P == True`。
- **对照（集合前 5%）**：在全部打分成功的位点中按 `max(erk1_score, erk2_score)` 降序，k = floor(0.05 × 打分成功数)，取 `max_score ≥ 第 k 名的值`（边界并列全纳入）；未打分位点该列为空串。
- **早峰**：响应行 `peak_tp_max ∈ {2min, 8min}`。位点口径：位点任一响应行早峰 → `early_* = True`；非响应位点 `early_*` 为空串。`resp_old`/`resp_new` = 位点在旧/新 R 中（True/False）。
- **表 E1**：位点口径 = 在 `resp_* == True` 的位点上按行定义过滤（全部 / 定义1 True / plus1_P True / 前5% True），早峰 n = 其中 `early_* == True` 数；行口径 = 旧/新 R 的响应行按 site_id 带上位点属性后同样过滤，早峰 n = 该行 peak_tp_max ∈ {2min, 8min} 的行数。比例 = 早峰 n / n（4 位小数）。未打分位点只出现在“全部响应位点”行。
- 只报计数与数值，不做检验。

## ③ (a) 3a / 3b 旧新并列

### 表 A1：peak_tp_max × type_90（行数）

| peak_tp_max | 旧 transient | 旧 sustained | 旧 90min缺失 | 旧 合计 | 新 transient | 新 sustained | 新 90min缺失 | 新 合计 |
|---|---|---|---|---|---|---|---|---|
| 2min | 704 | 255 | 316 | 1275 | 672 | 220 | 298 | 1190 |
| 8min | 768 | 428 | 261 | 1457 | 781 | 424 | 249 | 1454 |
| 20min | 513 | 298 | 224 | 1035 | 496 | 286 | 213 | 995 |
| 90min | 0 | 1882 | 0 | 1882 | 0 | 1780 | 0 | 1780 |
| 合计 | 1985 | 2863 | 801 | 5649 | 1949 | 2710 | 760 | 5419 |

### 位点级类型（位点数）

| site_type_90 | 旧 | 新 |
|---|---|---|
| transient | 1114 | 1153 |
| sustained | 1683 | 1628 |
| mixed | 538 | 483 |
| untyped | 452 | 447 |
| 合计 | 3787 | 3711 |

### 3b：n_frac_resp（位点数）

| n_frac_resp | 旧 | 新 |
|---|---|---|
| 1 | 2547 | 2566 |
| 2 | 814 | 756 |
| 3+ | 426 | 389 |
|   其中 3 | 283 | 258 |
|   其中 4 | 99 | 96 |
|   其中 5 | 35 | 27 |
|   其中 6 | 9 | 8 |
| 合计 | 3787 | 3711 |

### 3b：n_comp_resp（位点数）

| n_comp_resp | 旧 | 新 |
|---|---|---|
| 1 | 2874 | 2854 |
| 2 | 738 | 712 |
| 3 | 175 | 145 |
| 合计 | 3787 | 3711 |

- 新逐行表 `analysis/07_pilot2/04_time_rules_v2_rows.tsv`：5419 行，列同 06_pilot `03_time_pilots_rows.tsv`（`denom_tier` 全为 protein）。
- 旧列由同一算法对 site_traj 重算，逐行表与 `analysis/06_pilot/03_time_pilots_rows.tsv` 逐格相同：True。

## ④ (b) 窗口与打分统计

| 项 | 值 |
|---|---|
| site_traj 位点数 | 18268 |
| 取到 15-mer 窗口 | 18268（长度 15：18268；含 `_` 补齐 540） |
| acc 不在 FASTA / 位置越界 | 0 / 0 |
| 中心残基一致 / 不一致 | 18268 / 0 |
| 残基构成（全部） | S 14435，T 3241，Y 592 |
| +1 位是 P（定义 2） | 7148（+1 位为 `_` 补齐：15） |
| 送打分（center_ok） | 18268 |
| 实际接口 | `kl.PhosphoProteomics(df, seq_col='window15').predict(kin_type='ser_thr')` |
| 接口参数 | predict 全部取库默认值（kinases=None → 全部 ser_thr 激酶；st_fav=True；score_round_digits=3；percentile_round_digits=2）；PhosphoProteomics 默认 pp=False、drop_invalid_subs=True |
| 库内过滤 | 无效序列剔除 0；进入 ser_thr 集合 17676；中心为 Y、归入 tyrosine 集合（ser_thr 不打分）592 |
| ser_thr 激酶数（预测表） | 311（`get_kinase_list('ser_thr')` 311；含 ERK1 True、ERK2 True） |
| 耗时 | PhosphoProteomics 初始化 0.05 s；predict 50.18 s；批量接口合计 50.32 s；04_erk_score.py 全程 59.41 s |
| 打分成功（4 个 ERK 值均非空） | 17676（S 14435，T 3241） |
| 未打分 | 592（Y 592） |
| 逐条接口抽查 | 前 20 个打分位点用 kl.Substrate(window15, kin_type='ser_thr').predict()（库默认 score_round_digits=4、percentile_round_digits=2） 重算：abs(Δscore) 最大 0.0006，abs(Δpercentile) 最大 0.02 |
| ERK1 percentile ≥ 95 | 1362 |
| ERK2 percentile ≥ 95 | 1252 |
| 两者皆 ≥ 95 / 仅 ERK1 / 仅 ERK2 | 916 / 446 / 336 |
| ERK 候选（定义 1，并集） | 1698 |
| 集合前 5%（对照） | 883（k = 883，阈值 max_score ≥ 5.566，阈值处并列 1） |
| 交叠：定义1∩前5% / 定义1∩+1P / 前5%∩+1P | 883 / 1695 / 883 |
| +1P 且打分成功 | 7113 |

按响应集合（位点数）：

|  | 旧 R 位点 | 新 R 位点 |
|---|---|---|
| 响应位点 | 3787 | 3711 |
| 其中未打分（Y 中心） | 49（Y 49） | 50（Y 50） |
| ERK1 ≥ 95 | 491 | 470 |
| ERK2 ≥ 95 | 445 | 435 |
| ERK 候选（定义 1） | 606 | 588 |
| +1P（定义 2） | 2300 | 2273 |
| 集合前 5%（对照） | 336 | 325 |

## ⑤ 表 E1

### E1-位点口径（位点数；早峰 = 位点任一响应行 peak_tp_max ∈ {2min, 8min}）

|  | 旧 n | 旧 早峰 n | 旧 比例 | 新 n | 新 早峰 n | 新 比例 |
|---|---|---|---|---|---|---|
| 全部响应位点 | 3787 | 2179 | 0.5754 | 3711 | 2133 | 0.5748 |
| ERK 候选(定义1) | 606 | 376 | 0.6205 | 588 | 369 | 0.6276 |
| +1P(定义2) | 2300 | 1348 | 0.5861 | 2273 | 1325 | 0.5829 |
| 集合前5%(对照) | 336 | 220 | 0.6548 | 325 | 216 | 0.6646 |

### E1-行口径（响应行数；早峰 = 该行 peak_tp_max ∈ {2min, 8min}）

|  | 旧 n | 旧 早峰 n | 旧 比例 | 新 n | 新 早峰 n | 新 比例 |
|---|---|---|---|---|---|---|
| 全部响应位点 | 5649 | 2732 | 0.4836 | 5419 | 2644 | 0.4879 |
| ERK 候选(定义1) | 974 | 486 | 0.4990 | 932 | 475 | 0.5097 |
| +1P(定义2) | 3575 | 1711 | 0.4786 | 3463 | 1668 | 0.4817 |
| 集合前5%(对照) | 553 | 287 | 0.5190 | 518 | 282 | 0.5444 |

- 复算：位点口径 n 列由写出的 `analysis/07_pilot2/04_erk_sites.tsv` 过滤复算一致：True；行口径由 `analysis/07_pilot2/04_time_rules_v2_rows.tsv`（新）/ `analysis/06_pilot/03_time_pilots_rows.tsv`（旧）按 site_id 连 `analysis/07_pilot2/04_erk_sites.tsv` 复算一致：True。

## ⑥ 跳过 / 缺失项

- 中心为 Y 的 592 个位点：kinase-library `PhosphoProteomics` 把它们归入 tyrosine 集合，`predict(kin_type='ser_thr')` 不对其打分；这些位点 `erk1_*`/`erk2_*`/`erk_candidate_p95`/`erk_top5pct_set` 为空串，只计入 E1 的“全部响应位点”行。
- 中心残基不一致：0 个（无排除）。取不到窗口：0 个。
- 方案写的逐条调用 `kl.Substrate(seq).predict(kin_type='ser_thr')` 与 1.8.0 的签名不符（`kin_type` 在 `Substrate(...)` 构造参数里，`predict()` 无此参数）；本轮用批量接口，未走逐条回退；抽查用 `kl.Substrate(seq, kin_type='ser_thr').predict()`。
- 不做 backbone 对照、A2/A3、3b 交叉表（任务 4 未要求）；不做统计检验。
- `sites.tsv`、`summary*.txt`：本任务不依赖，未读。

### 附：验收自检（脚本自动核对）

| 项 | 实际 | 期望 | 通过 |
|---|---|---|---|
| 旧 R 行/位点/蛋白 | (5649, 3787, 1837) | (5649, 3787, 1837) | 是 |
| 旧 A1 2min 行（transient/sustained/90min缺失）= 上一轮 | [704, 255, 316] | [704, 255, 316] | 是 |
| 旧 A1 8min 行（transient/sustained/90min缺失）= 上一轮 | [768, 428, 261] | [768, 428, 261] | 是 |
| 旧 A1 20min 行（transient/sustained/90min缺失）= 上一轮 | [513, 298, 224] | [513, 298, 224] | 是 |
| 旧 A1 90min 行（transient/sustained/90min缺失）= 上一轮 | [0, 1882, 0] | [0, 1882, 0] | 是 |
| 旧 位点类型 transient/sustained/mixed/untyped = 上一轮 | [1114, 1683, 538, 452] | [1114, 1683, 538, 452] | 是 |
| 旧 n_frac_resp 1/2/3+ = 上一轮 | [2547, 814, 426] | [2547, 814, 426] | 是 |
| 旧 n_frac_resp 3/4/5/6 = 上一轮 | [283, 99, 35, 9] | [283, 99, 35, 9] | 是 |
| 旧 n_comp_resp 1/2/3 = 上一轮 | [2874, 738, 175] | [2874, 738, 175] | 是 |
| 旧逐行表与 06_pilot/03_time_pilots_rows.tsv 逐格相同 | True | True | 是 |
| 新旧表 (site_id, fraction) 顺序相同 | True | True | 是 |
| 新 R 行数 / 位点数 = 任务 1（C4 all / C3 n_new） | (5419, 3711) | (5419, 3711) | 是 |
| 新 A1 合计 = 新 R 行数 | 5419 | 5419 | 是 |
| 新 A1 90min缺失 列合计 = 新 R 中 log2occ_90min 为空行数 | 760 | 760 | 是 |
| 新 位点类型四类合计 = 新 R 位点数 | 3711 | 3711 | 是 |
| 新 n_frac_resp 1/2/3+ 合计 = 新 R 位点数 | 3711 | 3711 | 是 |
| 新 n_comp_resp 1/2/3 合计 = 新 R 位点数 | 3711 | 3711 | 是 |
| 新 peak_tp_max=90min 且 transient 行数 | 0 | 0 | 是 |
| 04_time_rules_v2_rows.tsv 行数 / 列名 | (5419, True) | (5419, True) | 是 |
| 04_time_rules_v2_rows.tsv 中 nan/NaN/None/NA 文本个数 | 0 | 0 | 是 |
| 04_erk_sites.tsv 行数 / 列名严格 | (18268, True) | (18268, True) | 是 |
| 04_erk_sites.tsv site_id 集合 = site_traj 位点集合 | True | True | 是 |
| 18268 位点全有 15-mer 窗口 | (18268, 18268) | (18268, 18268) | 是 |
| 中心残基不一致数（advisor 参考 0） | 0 | 0 | 是 |
| +1 位是 P 位点数（advisor 参考 7148） | 7148 | 7148 | 是 |
| 旧响应位点中 +1P 数（advisor 参考 2300） | 2300 | 2300 | 是 |
| resp_old / resp_new True 数 = 旧/新 R 位点数 | (3787, 3711) | (3787, 3711) | 是 |
| 04_erk_sites.tsv 布尔列只含 True/False/空串 | True | True | 是 |
| 04_erk_sites.tsv 中 nan/NaN/None/NA 文本个数 | 0 | 0 | 是 |
| 定义1 候选数 = max(pct_ERK1, pct_ERK2) ≥ 95 复算 | 1698 | 1698 | 是 |
| 打分日志与 tsv 一致：成功数 / ERK1≥95 / ERK2≥95 / 并集 / 前5% | (17676, 1362, 1252, 1698, 883) | (17676, 1362, 1252, 1698, 883) | 是 |
| E1 位点口径 n 列由 04_erk_sites.tsv 过滤复算一致 | True | True | 是 |
| E1 行口径 n 列由逐行表 + 04_erk_sites.tsv 复算一致 | True | True | 是 |
| E1 位点口径'全部响应位点' n = 旧/新 R 位点数 | (3787, 3711) | (3787, 3711) | 是 |
| E1 行口径'全部响应位点' n = 旧/新 R 行数 | (5649, 5419) | (5649, 5419) | 是 |
| E1 旧口径早峰数 = 上一轮基线（行 2732 / 位点 2179） | (2732, 2179) | (2732, 2179) | 是 |
| 全局 pandas/numpy/matplotlib 版本 开始前 = 结束后 | True | True | 是 |
