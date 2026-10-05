#!/usr/bin/env python3
"""02_run_medians.py — 任务 2：归一化趋势（run 级 median_log2）。

按 analysis/06_pilot/00_design.md「0. 通用约定」与「任务 2」执行。
输入（只读）: data/pilot/run_medians.tsv
输出: analysis/06_pilot/02_run_medians{.tsv,_summary.tsv,.md,.png,_phospho.png,_proteome.png}

用法: python3 analysis/06_pilot/02_run_medians.py   （在任意目录下可复跑）
"""
import re
import statistics
from collections import defaultdict
from pathlib import Path

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parents[2]
IN_PATH = REPO / "data" / "pilot" / "run_medians.tsv"
OUT_DIR = REPO / "analysis" / "06_pilot"
OUT_LONG = OUT_DIR / "02_run_medians.tsv"
OUT_SUM = OUT_DIR / "02_run_medians_summary.tsv"
OUT_MD = OUT_DIR / "02_run_medians.md"
OUT_PNG = OUT_DIR / "02_run_medians.png"
OUT_PNG_LAYER = {
    "phospho": OUT_DIR / "02_run_medians_phospho.png",
    "proteome": OUT_DIR / "02_run_medians_proteome.png",
}

LAYERS = ["phospho", "proteome"]
TPS = ["CTRL", "2min", "8min", "20min", "90min"]          # 07 :24
EGF_TPS = ["2min", "8min", "20min", "90min"]              # 07 :25
FRS = [f"FR{i}" for i in range(1, 7)]                     # 07 :26
REPS = [f"Rep{i}" for i in range(1, 5)]
DPI = 200

LONG_COLS = ["layer", "run", "timepoint", "fraction", "rep", "n_precursors",
             "median_log2", "anchor_layer", "offset_anchor", "offset_ctrl"]
SUM_COLS = ["layer", "fraction", "timepoint", "n_reps", "mean_median_log2",
            "sd_median_log2", "min_median_log2", "max_median_log2",
            "mean_n_precursors", "n_pos_anchor", "n_neg_anchor", "n_zero_anchor",
            "consistent_anchor", "n_pos_ctrl", "n_neg_ctrl", "n_zero_ctrl",
            "consistent_ctrl"]


# ---- 照抄 07_main_analysis.py parse_design (:32-38) ----
def parse_design(run):
    tp = re.search(r"_(2min|8min|20min|90min|CTRL)_", run, re.I)
    fr = re.search(r"_(FR\d)_", run)
    rep = re.search(r"_(Rep\d)", run)
    if not (tp and fr and rep):
        return None
    return (tp.group(1), fr.group(1), rep.group(1))


# ---- 照抄 07_main_analysis.py median (:100-101) ----
def median(v):
    return sorted(v)[len(v) // 2] if len(v) % 2 else sum(sorted(v)[len(v) // 2 - 1:len(v) // 2 + 1]) / 2


def fmt(x, nd=6):
    """数值输出：保留 nd 位小数；None -> 空串。"""
    if x is None:
        return ""
    return f"{x:.{nd}f}"


def sign_counts(vals):
    pos = sum(1 for v in vals if v > 0)
    neg = sum(1 for v in vals if v < 0)
    zero = sum(1 for v in vals if v == 0)
    return pos, neg, zero


def consistent(pos, neg, zero, n):
    # 全 >0 或全 <0 -> yes；出现 0 或符号混合 -> no
    return "yes" if n > 0 and (pos == n or neg == n) else "no"


# ======================================================================
# 1. 读表 + 解析
# ======================================================================
df = pd.read_csv(IN_PATH, sep="\t", dtype=str, keep_default_na=False)
n_input = len(df)

rows, failures = [], []
for i, r in df.iterrows():
    run = r["run"]
    layer = run.split("/", 1)[0] if "/" in run else ""
    parsed = parse_design(run)
    reasons = []
    if layer not in LAYERS:
        reasons.append(f"layer='{layer}'")
    if parsed is None:
        reasons.append("parse_design=None")
    try:
        med = float(r["median_log2"])
        npc = float(r["n_precursors"])
    except ValueError:
        reasons.append("non-numeric median_log2/n_precursors")
    if reasons:
        failures.append((i + 2, run, "; ".join(reasons)))   # +2: 表头 + 1-based 行号
        continue
    tp, fr, rep = parsed
    rows.append(dict(layer=layer, run=run, timepoint=tp, fraction=fr, rep=rep,
                     n_precursors_str=r["n_precursors"], n_precursors=npc,
                     median_log2_str=r["median_log2"], median_log2=med))

# ======================================================================
# 2. 偏移
# ======================================================================
anchor, n_anchor_runs = {}, {}
for layer in LAYERS:
    vals = [x["median_log2"] for x in rows if x["layer"] == layer]
    n_anchor_runs[layer] = len(vals)
    anchor[layer] = median(vals) if vals else None

ctrl_val = {}
dup_ctrl = []
for x in rows:
    if x["timepoint"] == "CTRL":
        key = (x["layer"], x["fraction"], x["rep"])
        if key in ctrl_val:
            dup_ctrl.append(key)
        ctrl_val[key] = x["median_log2"]

missing_ctrl = []
for x in rows:
    x["anchor_layer"] = anchor[x["layer"]]
    x["offset_anchor"] = x["median_log2"] - anchor[x["layer"]]
    if x["timepoint"] == "CTRL":
        x["offset_ctrl"] = None
    else:
        key = (x["layer"], x["fraction"], x["rep"])
        if key in ctrl_val:
            x["offset_ctrl"] = x["median_log2"] - ctrl_val[key]
        else:
            x["offset_ctrl"] = None
            missing_ctrl.append((x["run"], key))

# 排序：layer, fraction, timepoint(CTRL..90min), rep
order = lambda x: (LAYERS.index(x["layer"]), FRS.index(x["fraction"]) if x["fraction"] in FRS else 99,
                   TPS.index(x["timepoint"]) if x["timepoint"] in TPS else 99, x["rep"])
rows.sort(key=order)

long_df = pd.DataFrame([{
    "layer": x["layer"], "run": x["run"], "timepoint": x["timepoint"],
    "fraction": x["fraction"], "rep": x["rep"],
    "n_precursors": x["n_precursors_str"], "median_log2": x["median_log2_str"],
    "anchor_layer": fmt(x["anchor_layer"]), "offset_anchor": fmt(x["offset_anchor"]),
    "offset_ctrl": fmt(x["offset_ctrl"]),
} for x in rows], columns=LONG_COLS)
long_df.to_csv(OUT_LONG, sep="\t", index=False)

# ======================================================================
# 3. summary（60 cell）
# ======================================================================
cell_rows = defaultdict(list)
for x in rows:
    cell_rows[(x["layer"], x["fraction"], x["timepoint"])].append(x)

summary = []
for layer in LAYERS:
    for fr in FRS:
        for tp in TPS:
            cr = cell_rows.get((layer, fr, tp), [])
            n = len(cr)
            meds = [x["median_log2"] for x in cr]
            npcs = [x["n_precursors"] for x in cr]
            pa, na, za = sign_counts([x["offset_anchor"] for x in cr])
            s = dict(layer=layer, fraction=fr, timepoint=tp, n_reps=n,
                     mean_median_log2=statistics.mean(meds) if n else None,
                     sd_median_log2=statistics.stdev(meds) if n >= 2 else None,   # ddof=1
                     min_median_log2=min(meds) if n else None,
                     max_median_log2=max(meds) if n else None,
                     mean_n_precursors=statistics.mean(npcs) if n else None,
                     n_pos_anchor=pa, n_neg_anchor=na, n_zero_anchor=za,
                     consistent_anchor=consistent(pa, na, za, n))
            if tp == "CTRL":
                s.update(n_pos_ctrl=None, n_neg_ctrl=None, n_zero_ctrl=None, consistent_ctrl=None)
            else:
                pc, nc, zc = sign_counts([x["offset_ctrl"] for x in cr if x["offset_ctrl"] is not None])
                s.update(n_pos_ctrl=pc, n_neg_ctrl=nc, n_zero_ctrl=zc,
                         consistent_ctrl=consistent(pc, nc, zc, n))
            summary.append(s)


def sum_out(s):
    o = {}
    for c in SUM_COLS:
        v = s[c]
        if v is None:
            o[c] = ""
        elif c == "mean_n_precursors":
            o[c] = fmt(v, 2)
        elif isinstance(v, float):
            o[c] = fmt(v)
        else:
            o[c] = str(v)
    return o


sum_df = pd.DataFrame([sum_out(s) for s in summary], columns=SUM_COLS)
sum_df.to_csv(OUT_SUM, sep="\t", index=False)

# ======================================================================
# 4. 一致性计数
# ======================================================================
cons_rows = []
tot = defaultdict(int)
for layer in LAYERS:
    for fr in FRS:
        ss = [s for s in summary if s["layer"] == layer and s["fraction"] == fr]
        ya = [s["timepoint"] for s in ss if s["consistent_anchor"] == "yes"]
        yc = [s["timepoint"] for s in ss if s["consistent_ctrl"] == "yes"]
        za = sum(1 for s in ss if s["n_zero_anchor"] > 0)
        zc = sum(1 for s in ss if s["n_zero_ctrl"] not in (None, 0))
        cons_rows.append((layer, fr, len(ya), ya, len(yc), yc, za, zc))
        for k, v in (("a", len(ya)), ("c", len(yc)), ("za", za), ("zc", zc)):
            tot[(layer, k)] += v
            tot[("all", k)] += v
n_zero_anchor_vals = sum(s["n_zero_anchor"] for s in summary)
n_zero_ctrl_vals = sum(s["n_zero_ctrl"] or 0 for s in summary)

# ======================================================================
# 5. md
# ======================================================================
DEF_TEXT = """\
（照抄 `analysis/06_pilot/00_design.md` 任务 2「解析」「偏移定义」两节）

**解析（照 07 `parse_design`，`:32-38`）**
- `layer` = `run` 在第一个 `/` 之前的部分（phospho | proteome）。
- `timepoint` = 正则 `_(2min|8min|20min|90min|CTRL)_`（忽略大小写）；`fraction` = `_(FR\\d)_`；`rep` = `_(Rep\\d)`。解析失败的行计数并列出（预期 0）。
- 预期 2 layer × 5 tp × 6 FR × 4 rep = 240，每个 cell 恰 4 行。

**偏移定义（两种都算，都只是计数）**
- `anchor_{layer}` = 该 layer 120 个 `median_log2` 的中位数（07 `factors()` `:239-242` 的 anchor）。`offset_anchor` = `median_log2 − anchor_layer`。
- `offset_ctrl` = `median_log2 − median_log2(同 layer、同 fraction、同 rep 的 CTRL run)`；CTRL 行本身留空。
- "方向一致" := 同一 layer×fraction×timepoint 的 4 个重复的 offset 符号全同（全 >0 或全 <0；出现 0 记为不一致，并单独计数）。

**summary 列**：sd 用样本标准差（ddof=1）。CTRL 行的 `*_ctrl` 列留空。

执行细节（不改定义，仅说明实现）：
- 中位数用 07 `median()`（`:100-101`）的同一写法；偶数个取中间两值的平均。
- 正负号在未取整的 float 差值上判断；`n_zero_*` 计 offset 恰等于 0 的行。
- `consistent_*` 取值 `yes` / `no`；CTRL 行 `consistent_ctrl` 为空。
- 输出数值保留 6 位小数（`mean_n_precursors` 2 位）；`median_log2`、`n_precursors` 原样照抄输入字符串。
"""


def md_table(header, body):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for b in body:
        out.append("| " + " | ".join(str(v) for v in b) + " |")
    return "\n".join(out)


md = []
md.append("# 02 归一化趋势：run 级 median_log2\n")
md.append("## ① 输入\n")
md.append(md_table(["文件", "行数(不含表头)", "列"],
                   [[f"`{IN_PATH.relative_to(REPO)}`", n_input, ", ".join(f"`{c}`" for c in df.columns)]]))
md.append("")
md.append(f"- 解析成功 {len(rows)} 行；解析失败 {len(failures)} 行。")
for layer in LAYERS:
    md.append(f"- `{layer}`：{sum(1 for x in rows if x['layer'] == layer)} 行。")
md.append(f"- 输出：`{OUT_LONG.name}`（{len(long_df)} 行）、`{OUT_SUM.name}`（{len(sum_df)} 行）、"
          f"`{OUT_PNG.name}`、`{OUT_PNG_LAYER['phospho'].name}`、`{OUT_PNG_LAYER['proteome'].name}`、`02_run_medians.py`。\n")
md.append("## ② 定义\n")
md.append(DEF_TEXT)
md.append("## ③ anchor 值\n")
md.append(md_table(["layer", "n_runs", "anchor_layer"],
                   [[l, n_anchor_runs[l], fmt(anchor[l])] for l in LAYERS]))
md.append("")
md.append("## ④ summary 表（全文，60 行；与 `02_run_medians_summary.tsv` 相同）\n")
md.append(md_table(SUM_COLS, sum_df.values.tolist()))
md.append("")
md.append("## ⑤ 一致性计数\n")
md.append("每个 layer×fraction：`consistent_anchor=yes` 的时间点数（/5，含 CTRL）、`consistent_ctrl=yes` 的时间点数（/4，仅 EGF 时间点）；"
          "`zero_cells_*` = 该 layer×fraction 中出现 ≥1 个 offset=0 的时间点数（这些 cell 记为不一致）。\n")
body = []
for (layer, fr, na, ya, nc, yc, za, zc) in cons_rows:
    body.append([layer, fr, f"{na}/5", ",".join(ya) or "—", f"{nc}/4", ",".join(yc) or "—", za, zc])
md.append(md_table(["layer", "fraction", "consistent_anchor=yes", "timepoints(anchor)",
                    "consistent_ctrl=yes", "timepoints(ctrl)", "zero_cells_anchor", "zero_cells_ctrl"], body))
md.append("")
md.append(md_table(["范围", "consistent_anchor=yes", "consistent_ctrl=yes", "zero_cells_anchor", "zero_cells_ctrl"],
                   [[l, f"{tot[(l, 'a')]}/30", f"{tot[(l, 'c')]}/24", tot[(l, 'za')], tot[(l, 'zc')]] for l in LAYERS]
                   + [["总计", f"{tot[('all', 'a')]}/60", f"{tot[('all', 'c')]}/48", tot[('all', 'za')], tot[('all', 'zc')]]]))
md.append("")
md.append(f"- offset_anchor 恰为 0 的行数：{n_zero_anchor_vals}（/240）；offset_ctrl 恰为 0 的行数：{n_zero_ctrl_vals}（/192）。\n")
md.append("## ⑥ 解析失败列表 / 跳过 / 缺失项\n")
if failures:
    md.append(md_table(["行号(含表头计)", "run", "原因"], failures))
else:
    md.append("- 解析失败：0 行（列表为空）。")
md.append(f"- 缺 CTRL 对照导致 offset_ctrl 无法计算的 EGF 行：{len(missing_ctrl)}。")
md.append(f"- 重复的 CTRL 键（layer, fraction, rep）：{len(dup_ctrl)}。")
bad_cells = [(s["layer"], s["fraction"], s["timepoint"], s["n_reps"]) for s in summary if s["n_reps"] != 4]
md.append(f"- n_reps ≠ 4 的 cell：{len(bad_cells)}。")
md.append("")
OUT_MD.write_text("\n".join(md), encoding="utf-8")

# ======================================================================
# 6. 图（dataviz 规范：参考调色板 slot 1-4 固定顺序；浅色 surface；
#    细线 + ≥8px 带 surface 环的 marker；发丝网格；文字用 ink token）
# ======================================================================
SURFACE = "#fcfcfb"
INK_PRIMARY, INK_SECONDARY, INK_MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
REP_COLOR = dict(zip(REPS, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]))
REP_MARKER = dict(zip(REPS, ["o", "s", "^", "D"]))          # 次级编码（形状）
LW = 1.5                                                      # pt；细线，round join/cap
MS = 5.5                                                      # pt；marker 直径 ≥8px @ DPI

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK_SECONDARY,
    "xtick.color": AXIS, "ytick.color": AXIS,
    "xtick.labelcolor": INK_MUTED, "ytick.labelcolor": INK_MUTED,
    "font.size": 9,
})

ylim = {}
for layer in LAYERS:
    v = [x["median_log2"] for x in rows if x["layer"] == layer]
    lo, hi = min(v), max(v)
    pad = (hi - lo) * 0.06 or 0.5
    ylim[layer] = (lo - pad, hi + pad)


def draw_panel(ax, layer, fr, show_ylabel):
    for rep in REPS:
        ys = []
        for tp in TPS:
            m = [x["median_log2"] for x in rows
                 if x["layer"] == layer and x["fraction"] == fr and x["timepoint"] == tp and x["rep"] == rep]
            ys.append(m[0] if m else float("nan"))
        ax.plot(range(len(TPS)), ys, color=REP_COLOR[rep], lw=LW, solid_joinstyle="round",
                solid_capstyle="round", marker=REP_MARKER[rep], ms=MS,
                markeredgecolor=SURFACE, markeredgewidth=0.9, zorder=3)
    ax.set_xticks(range(len(TPS)))
    ax.set_xticklabels(TPS)
    ax.set_xlim(-0.35, len(TPS) - 0.65)
    ax.set_ylim(*ylim[layer])
    ax.grid(axis="y", color=GRID, lw=0.6, ls="-", zorder=0)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.6)
    ax.tick_params(length=2.5, width=0.6)
    ax.set_title(f"{layer} · {fr}", color=INK_PRIMARY, fontsize=10, loc="left")
    if show_ylabel:
        ax.set_ylabel("median_log2")


def legend_handles():
    return [Line2D([0], [0], color=REP_COLOR[r], lw=LW, marker=REP_MARKER[r], ms=MS,
                   markeredgecolor=SURFACE, markeredgewidth=0.9, label=r) for r in REPS]


def add_legend(fig):
    fig.legend(handles=legend_handles(), loc="upper right", ncol=4, frameon=False,
               labelcolor=INK_SECONDARY, fontsize=9, handlelength=2.2,
               bbox_to_anchor=(0.995, 0.995))


# 2×6 合图（上 phospho、下 proteome；每行共享 y）
fig, axes = plt.subplots(2, 6, figsize=(18, 6.6), sharey="row")
for ri, layer in enumerate(LAYERS):
    for ci, fr in enumerate(FRS):
        draw_panel(axes[ri, ci], layer, fr, ci == 0)
fig.suptitle("Run-level median_log2 by timepoint (lines = replicates)", x=0.005, ha="left",
             color=INK_PRIMARY, fontsize=11)
add_legend(fig)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(OUT_PNG, dpi=DPI)
plt.close(fig)

# 单 layer 1×6
for layer in LAYERS:
    fig, axes = plt.subplots(1, 6, figsize=(18, 3.6), sharey=True)
    for ci, fr in enumerate(FRS):
        draw_panel(axes[ci], layer, fr, ci == 0)
    fig.suptitle(f"{layer}: run-level median_log2 by timepoint (lines = replicates)", x=0.005,
                 ha="left", color=INK_PRIMARY, fontsize=11)
    add_legend(fig)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(OUT_PNG_LAYER[layer], dpi=DPI)
    plt.close(fig)

print(f"rows parsed={len(rows)} failures={len(failures)} "
      f"anchor={ {l: fmt(anchor[l]) for l in LAYERS} } "
      f"consistent_anchor={tot[('all', 'a')]}/60 consistent_ctrl={tot[('all', 'c')]}/48")
