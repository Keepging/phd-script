# 09 验收记录（advisor，gate1）— 2026-10-06 重跑版

## 本轮经过
1. 首轮（commit 157284f）：`data/pilot/gate1_proteins.tsv` 只有表头、0 行，三个任务由三个 Opus subagent 写好脚本并在 scratchpad 合成夹具上验证逻辑（advisor 用独立实现逐格比对，全部一致），真实结果全部"缺数据，跳过"。
2. 重新导出的数据以 `analysis/gate1_proteins.tsv`（344 行，commit 9af8707 "Add files via upload"）进入当前分支；`data/pilot/gate1_proteins.tsv` 仍为 0 行，`data/` 未改。三个脚本通过 `--proteins analysis/gate1_proteins.tsv`（03 为 `--gate1`）按 01 → 03 → 02 重跑，输出覆盖 `analysis/08_gate1/`。
3. 两处定义按确认执行：centered = quantity × 2^a（a 为 07 实际施加的 log2 因子，`01_run_factors.tsv` 的 `a` 列）；Movement Score 阈值改为作者 :181 的严格 `> 0.1`（`02_events.py` 与 `00_design.md` 已改，列名 `gt_0.1`）。
4. 验收：advisor 用独立实现在真实数据上重算份额（160 行）、EGFR 总量（40 行）与比值（8 行）、E-A/E-B/E-C/E-D（328 项），与脚本输出 0 处不符（份额 1e-9、总量 3e-8、MS 1e-6、p 1e-5 容差）。

## 输入事实
- 344 行，119 个 distinct run（120 个设计格中 119 个有行），解析失败 0；run 名无 `proteome/` 前缀，按 (timepoint, fraction, rep) 匹配因子。
- 每个基因只有一个 protein_group，无 secondary：GRB2 P62993、SHC1 P29353、CBL P22681、EGFR P00533。
- **四个蛋白测到的 run 数 /120**：GRB2 90、SHC1 115、CBL 59、EGFR 80。
- 每个基因 20 个 (timepoint, rep) 格都有 ≥1 个 fraction；**fraction 不全（<6）的格**：SHC1 5/20、GRB2 18/20、CBL 20/20、EGFR 20/20。份额按有值 fraction 之和归一，缺的在 `01_shares.md` ③ 逐行列出，未填补。

## 任务 1 份额表 — 合格
`01_shares.tsv` 160 行（4 基因 × 5 tp × 4 rep × raw/centered），三格非空行 Cyt+Mem+Nuc = 1（1e-12 内）；`01_shares_secondary.tsv` 只有表头；`01_groups.tsv` 4 行。

## 任务 2 三个易位事件 — 合格
- E-B 方向计数（Mem 上升且 Cyt 下降，按 rep 配对 vs CTRL；raw 与 centered 相同）：
  | gene | 2 min | 8 min |
  |---|---|---|
  | GRB2 | 4/4 | 4/4 |
  | SHC1 | 4/4 | 4/4 |
  | CBL | 3/4 | 2/3 |
- E-C Movement Score（作者定义）：可算的只有 SHC1 全部 8 格和 GRB2 8min 两格；GRB2 的 2/20/90 min 与 CBL 全部时间点因某 fraction 在 tp 或 CTRL 的 4 个 rep 全缺而留空（方案规定不在少于 6 个 fraction 上算）。**MS_author > 0.1 的 3 个**：(GRB2, raw, 8min) 0.1298；(GRB2, centered, 8min) 0.1288；(SHC1, centered, 2min) 0.1013。SHC1 raw 2min = 0.0981（低于阈值）。MS_max1 在上述格均为 FR4，MS_max2 为 FR1。
- E-A 均值/标准差/Δ份额 90 行、E-D Welch p 值 24 行见 `02_events.md`（n = 4/4 或实际 n，只作参考）。
- 图 `02_events.png`：2 × 9，缺份额处断开。

## 任务 3 EGFR 总量 — 合格
EGFR 90 min / CTRL 总量比值（按 rep 配对）：
| rep | raw | centered |
|---|---|---|
| Rep1 | 0.471 | 0.409 |
| Rep2 | 0.247 | 0.272 |
| Rep3 | 0.397 | 0.360 |
| Rep4 | 1.468 | 1.694 |
EGFR 的 20 个 (tp, rep) 格全部 fraction 不全，总量是有值 fraction 之和，详见 `03_egfr_total.md`。

## 改过什么
- 阈值 ≥ 0.1 → > 0.1（`00_design.md`、`02_events.py`，列名 `ge_0.1` → `gt_0.1`）。
- 输入路径：脚本默认仍指向 `data/pilot/gate1_proteins.tsv`（0 行），本轮用参数指向 `analysis/gate1_proteins.tsv`；若以后把导出放回 `data/pilot/`，不加参数即可。
- 本轮重跑由 advisor 直接执行已验收的脚本并独立核对，未再派 subagent。

## 没解决 / 需要知道的问题
1. **覆盖不全是主要限制**：CBL 只在 59/120 个 run 测到，所有 (tp, rep) 格都缺 fraction，Movement Score 全部无法按作者定义计算；GRB2 只有 8 min 可算。份额与方向计数在缺 fraction 的格上是按"有值 fraction 之和"归一的，分母随缺失情况变化。
2. 运行 `FR?`：120 个设计格中 1 个没有任何行（见 `01_shares.md` ①）。
3. 作者的 `pval_combi_FDR < 0.05`（limma + sumlog + BH）条件未复现；DAPAR 的 LOESS 归一化与填补未复现。
4. 任务 3 对 `n_fractions_present = 0` 的格出空行、任务 1 不出行；本数据下每格都 ≥1 个 fraction，两表行数一致。
