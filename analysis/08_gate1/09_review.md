# 09 验收记录（advisor，gate1）

方案：`00_design.md`。执行：三个 Opus subagent（任务 1 与 3 并行，任务 2 在任务 1 通过后派）。三个任务一轮通过，没有重派。

## 核心事实：输入为空
`data/pilot/gate1_proteins.tsv`（分支 `claude/determined-goldberg-atmb8n` commit bee801a，已合并进当前分支，`analysis/gate1_proteins.tsv` 为同一文件副本）**只有表头、0 行数据**。因此：
- GRB2 / SHC1 / CBL / EGFR 各在 **0/120** 个 run 测到。
- 任务 1、2、3 的全部数值结果：**缺数据，跳过**。三个任务的表只有表头，`02_events.png` 为 18 个写着 "no data" 的空图框。
- 三个脚本（`01_shares.py`、`02_events.py`、`03_egfr_total.py`）已按方案写好并对着真实文件跑过（退出码 0，复跑逐字节相同）；数据到位后按 01 → 03 → 02 的顺序直接复跑即可，不需改动。

## 验收方法
真实输入无数可对，advisor 改为验证脚本逻辑：在 scratchpad 建了一个明显合成的夹具（4 基因 × 120 run，含故意缺失的 fraction、一个基因两个 protein_group），用独立实现算出参考值，再用三个脚本的可选参数 `--proteins/--shares/--factors/--outdir` 在 scratchpad 跑夹具比对。夹具与其任何数字都不进入 `analysis/`。

## 任务 1 份额表 — 合格（脚本逻辑）
文件：`01_shares.tsv`、`01_shares_secondary.tsv`、`01_groups.tsv`（均只有表头）、`01_shares.md`、`01_shares.py`。
夹具比对：160 行（4 基因 × 5 tp × 4 rep × raw/centered）Cyt/Mem/Nuc 与 advisor 参考逐格一致（0 处不符，含留空格）；`n_fractions_present` 一致；三格非空行 Cyt+Mem+Nuc 与 1 的最大偏差 1e-12；primary/secondary 分流正确（run 数多者为主）；centered = quantity × 2^a 正确。真实输入：md ① 写明 0/120、解析 0/120。

## 任务 2 三个易位事件 — 合格（脚本逻辑）
文件：`02_events.md`、`02_events.tsv`（只有表头）、`02_events.png`（2×9 空图框）、`02_events.py`。
夹具比对：E-A 90 格 × (mean, sd, Δ) 全部一致（1e-5）；E-B 12 行 × 4 个计数一致（`n/4` 格式，分母为两边四格都有值的 rep 数）；E-C 24 行 MS_author（作者定义：fraction 级线性重复均值 → 份额 → 6 个 |Δ| 的均值）、MS_max1、MS_max2、ge_0.1 一致；E-D 24 行 Welch t、p 一致。夹具图目视：2 行 × 9 列、每图 4 条重复线、缺份额处断开、图例一次。
subagent 记录的口径选择（advisor 认可）：阈值按方案用 ≥ 0.1（作者第 181 行是 > 0.1，只在恰等于 0.1 时不同）；某 fraction 在 tp 或 CTRL 的 4 个 rep 全缺时 E-C 该行留空并列出；`n_mem_up`/`n_cyt_down` 为单项计数，排他计数可由减法得到；作者的 `pval_combi_FDR < 0.05` 条件未复现。

## 任务 3 EGFR 总量 — 合格（脚本逻辑）
文件：`03_egfr_total.tsv`、`03_egfr_ratio.tsv`（只有表头）、`03_egfr_total.md`、`03_egfr_total.py`。
夹具比对：40 行 total（raw/centered）与参考最大差 5e-10；`n_fractions_present` 一致；8 个按 rep 配对的 90min/CTRL 比值全部一致。`02_events.md` ⑦ 节已抄入其 ③④ 节。
subagent 口径：`n_fractions_present = 0` 的 (tp, rep) 也出一行（total 留空），与任务 1"一个 fraction 都没有则不出行"不一致，数据到位后两表行数会差这一类行；其余一致。

## 改过什么
- 方案 §1 引用作者代码第 181 行时漏了行尾 `+`（ggplot 续行符），subagent 已在 md 注明；方案文本未再改。
- 无重派，输出未被 advisor 手改。

## 没解决 / 需要知道的问题
1. **输入为空是本轮唯一阻塞项**：`gate1_proteins.tsv` 需要重新导出（每个 proteome run 的 `report.pg_matrix.tsv` 最后一列）。导出后复跑三个脚本即可得到全部结果。
2. **centered 的定义**：按 `quantity × 2^a`（`a` 为 `01_run_factors.tsv` 中 07 实际使用的 log2 因子，= anchor − m_07）。任务原文写"除以 2^该 run 的 factor"，若 factor 指 m_07 − anchor，两者等价；若指别的量，需要改一行。
3. **Movement Score 口径**（方案 §1）：作者只对 2 min vs CTRL 算，且输入是 DAPAR 按 fraction LOESS 归一化并填补后的值；本轮对四个 EGF 时间点都算，且 raw/centered 两版都不复现 DAPAR 归一化与填补，缺失 rep 直接从均值里剔除。
4. 份额的归一化基数是"有值的 fraction 之和"，缺 fraction 的 (tp, rep) 在 md 单列；没有填补。
5. 行数口径小差异：任务 3 对 `n_fractions_present = 0` 的格出空行，任务 1 不出行。
