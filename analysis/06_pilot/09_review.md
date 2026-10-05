# 09 验收记录（advisor）

方案：`00_design.md`。执行：三个 Opus subagent，各自只拿到本任务的方案节和输入路径。验收方法：advisor 在分派前用独立实现（scratchpad，不入库）算出一套参考数，交回后逐格比对；再按方案的验收标准逐条核。三个任务都是一轮通过，没有重派。

## 输入与分支说明
- `data/pilot/` 原本只在分支 `claude/determined-goldberg-atmb8n`（commit 856e6c5）上，本分支没有。已把该分支合并进当前分支（merge commit 556bdf1），`data/` 内容未改。`analysis/01–04` 也随合并进入本分支。
- `kinase_pssm/` 不存在；`site_features.tsv` 存在但无序列列。

## 任务 1 漏斗表 — 合格
文件：`01_funnel.tsv`（70 行 × 8 列）、`01_funnel.md`、`01_funnel.py`。
验收：
- 70 个 step×scope 格的 n_rows / n_sites / n_proteins / site_log2int_CTRL_median 与 advisor 参考值全部一致（0 处不符）。
- step1/all = 39813/18268/4996；step2 n_rows 14963；step3 n_rows 13106 = ctrl_a 9037 + ctrl_b 4069；ctrl_d 16854：全部对上。
- all 口径 6 步单调不增；每步 6 个 fraction 的 n_rows 之和 = all；`protein_abundance_median` 全部为"无"（输入无蛋白级强度列）。
- 复跑（在降级后的 pandas 2.2.3 下）tsv/md md5 不变。
subagent 自定的点（advisor 认可，与方案意图一致）：c5/c6 在整个 S3 上计数再与上一步剩余取交；ctrl_d 的中位数格填"无"；对照行按 fraction 拆分用"该 fraction 内满足条件的行"。
需注意的事实：ctrl_a 与 ctrl_b 的 n_sites 不可相加（4704+2397 > 6439），同一位点在不同 fraction 可属不同 denom_tier。

## 任务 2 归一化趋势 — 合格
文件：`02_run_medians.tsv`（240 行）、`02_run_medians_summary.tsv`（60 行 × 17 列）、`02_run_medians.md`、`02_run_medians.png`（2×6）、`02_run_medians_phospho.png`、`02_run_medians_proteome.png`（各 1×6）、`02_run_medians.py`。
验收：
- 解析失败 0；60 个 cell 每个 n_reps=4；anchor phospho 18.84585 / proteome 22.671 与参考一致。
- 抽查 4 个 cell（FR3 × CTRL/90min × 两 layer）mean/sd 与参考一致（容差 1e-6）；offset_ctrl 在 48 个 CTRL 行为空、192 个 EGF 行非空，抽查 2 行数值正确。
- 12 个 layer×fraction 的 consistent_anchor / consistent_ctrl 计数与参考逐一相同；总计 42/60、6/48。
- 三张 png 存在、可打开；合图目视确认 2 行 × 6 列、每图 4 条线、图例一次、同 layer 共享 y。
- 复跑 tsv/md md5 不变；png 字节因 matplotlib 版本变化而不同（已保留首次审过的版本）。
subagent 报告的偏离：dataviz skill 的两两色差校验中 Rep2↔Rep4（orange↔yellow）ΔE 13.7 < 15 判 FAIL；处理为保留 4 条线并给每个重复加不同 marker（o/s/^/D），未换色板、未减线数；仅出浅色版。advisor 接受，不影响数值。

## 任务 3 时间规律 — 合格
文件：`03_time_pilots_rows.tsv`（5649 行 × 13 列）、`03_time_pilots_sites.tsv`（3787 行 × 10 列）、`03_time_pilots.md`、`03_time_pilots.py`。
验收：
- R = 5649 行 / 3787 位点 / 1837 蛋白；表 A1 12 格全部与参考一致；90min缺失 801；peak_tp_max=90min 且 transient 的行 0。
- 位点类型 transient 1114 / sustained 1683 / mixed 538 / untyped 452（和 3787）；n_frac_resp 2547/814/283/99/35/9；n_comp_resp 2874/738/175；n_frac_evaluable 边际 1406/956/578/350/252/245：全部与参考一致。
- backbone 对照（single）：347 蛋白，d 中位数 0.0185，d>0 198 / d<0 148 / d=0 1；（win5）347，0.0087，200/147/0：与参考一致。缺 backbone 的响应位点 5 个。
- 3c：三张表无序列列，仓库无 fasta / biophys_json，kinase-library 打分与"+1 位是 P"代理均记"缺序列，跳过"；pip install 成功、`import kinase_library` 1.8.0 成功、ser_thr 列表 311 个含 ERK1/ERK2；基线比例行口径 2732/5649 = 0.4836、位点口径 2179/3787 = 0.5754。
- rows/sites 两表无 NaN 文本；`fc_at_peak` 与 `fc_{peak_tp_max}` 逐行一致。

## 改过什么
1. 方案 `00_design.md` 在执行中改了一处事实：`site_features.tsv` 的缺失值是字面字符串 `NA`（single/win5 各列 35 行，disoMine 46 行；响应位点里 backbone 为 NA 的 5 个），原稿写"预期 0"。已通知任务 3 subagent 按缺失处理，方案和验收标准同步改为 5。
2. 无重派。三个任务的输出未被 advisor 手改。

## 没解决 / 需要知道的问题
1. **蛋白丰度**：三张输入表都没有蛋白级强度列，漏斗表该列全为"无"；只给了位点级 `log2int_CTRL` 中位数作对照，不是蛋白丰度。
2. **3c 无法执行**：没有序列来源（无 FASTA、无 biophys_json、site_features 无窗口列），ERK 候选比例无法算；只有基线比例。
3. **环境副作用**：任务 3 的 `pip install kinase-library` 把全局 Python 包降级：pandas 3.0.6→2.2.3、numpy 2.4.6→1.26.4、matplotlib 3.11.2→3.8.4，并装入 scikit-learn、statsmodels 等。任务 1/2/3 的 tsv/md 在新版本下复跑结果字节相同；只有 png 字节不同。
4. **c3 的口径**：方案把"有蛋白档分母"定义为"occupancy 值存在（residue 或 protein 档均算）"，并用 ctrl_a/ctrl_b 拆开两档。若原意是仅 protein 档，step3 应取 ctrl_a 的数（9037 行 / 4704 位点 / 2115 蛋白），后续步骤需重算。
5. **peak 的两种口径**：07 的 `peak_tp` 是第一个达标时间点，不是峰值；3a 以 |fc| 最大的 `peak_tp_max` 为主，07 口径另列为表 A2，两者交叉表 A3 对角线 3951 / 非对角 1698 行。
6. **A1 与 A2 的验收措辞**：方案写"A1 与 A2 行合计相同"不准确（两表行变量不同），实际核的是总计与 type_90 列合计相同。
7. **PCT=90** 在本轮三个任务中都未用到。
8. 任务 2 的 png 只有浅色版；配色两两校验有一对 ΔE 低于 dataviz 门槛（见任务 2）。
