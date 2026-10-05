#!/usr/bin/env python3
"""03_ctrl_qc.py — 任务 3：CTRL 质检（run 级 n_precursors / median_log2 的稳健 z）。

按 analysis/07_pilot2/00_design.md 开头通用约定与「任务 3：CTRL 质检」执行。
输入（只读）: data/pilot/run_medians.tsv
输出: analysis/07_pilot2/03_ctrl_qc.tsv, 03_ctrl_qc_cells.tsv, 03_ctrl_qc.md

定义：
  组(cell) = 同 layer x fraction 的 20 个 run（5 tp x 4 rep，含非 CTRL）。
  对 n_precursors 与 median_log2 各算：
    med = median(20 值)；MAD = median(|x - med|)；z = (x - med) / (1.4826 * MAD)
  MAD = 0 时 z 留空并记录。只输出 CTRL run 的行。
  flag = |z| > 2 -> yes，否则 no（z 为空时 flag 也留空）。

注：方案写「24 个 CTRL run（phospho 12 + proteome 12）」，但输入中 CTRL run 实为
2 layer x 6 FR x 4 rep = 48 个（phospho 24 + proteome 24）。本脚本输出全部 48 个，
不做任意子集挑选；自检中「24 行」一项按方案原文判定并如实记为 FAIL。

用法: python3 analysis/07_pilot2/03_ctrl_qc.py   （在任意目录下可复跑，确定）
"""
import re
import statistics
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
IN_PATH = REPO / "data" / "pilot" / "run_medians.tsv"
OUT_DIR = REPO / "analysis" / "07_pilot2"
OUT_TSV = OUT_DIR / "03_ctrl_qc.tsv"
OUT_CELLS = OUT_DIR / "03_ctrl_qc_cells.tsv"
OUT_MD = OUT_DIR / "03_ctrl_qc.md"

LAYERS = ["phospho", "proteome"]
EGF_TPS = ["2min", "8min", "20min", "90min"]
TPS = ["CTRL"] + EGF_TPS
FRS = [f"FR{i}" for i in range(1, 7)]
REPS = [f"Rep{i}" for i in range(1, 5)]
MAD_K = 1.4826
Z_THRESH = 2.0
ND = 4  # 保留 4 位小数

QC_COLS = ["layer", "fraction", "rep", "run", "n_precursors", "median_log2",
           "cell_median_n", "cell_mad_n", "z_n_precursors",
           "cell_median_log2", "cell_mad_log2", "z_median_log2",
           "flag_n", "flag_log2"]
CELL_COLS = ["layer", "fraction", "n_runs", "median_n", "mad_n",
             "median_log2", "mad_log2"]


# ---- 照抄 code/Protein contour/Zhihan/Test/07_main_analysis.py parse_design (:32-38) ----
def parse_design(run):
    tp = re.search(r"_(2min|8min|20min|90min|CTRL)_", run, re.I)
    fr = re.search(r"_(FR\d)_", run)
    rep = re.search(r"_(Rep\d)", run)
    if not (tp and fr and rep):
        return None
    return (tp.group(1), fr.group(1), rep.group(1))


def fmt(x, nd=ND):
    return "" if x is None else f"{x:.{nd}f}"


def robust(vals):
    """返回 (med, MAD)。"""
    med = statistics.median(vals)
    mad = statistics.median([abs(v - med) for v in vals])
    return med, mad


def zscore(x, med, mad):
    if mad == 0:
        return None
    return (x - med) / (MAD_K * mad)


def flag(z):
    if z is None:
        return ""
    return "yes" if abs(z) > Z_THRESH else "no"


def md_table(df, cols):
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def main():
    raw = pd.read_csv(IN_PATH, sep="\t", dtype=str, keep_default_na=False)
    n_input = len(raw)

    # ---- 解析 ----
    recs, parse_fail = [], []
    for _, r in raw.iterrows():
        run = r["run"]
        layer = run.split("/", 1)[0]
        pd_ = parse_design(run)
        if pd_ is None or layer not in LAYERS:
            parse_fail.append(run)
            continue
        tp, fr, rep = pd_
        # 07 正则带 re.I；统一成 TPS 中的写法（本数据全为 'CTRL'，此处只为稳妥）
        tp_norm = next((t for t in TPS if t.lower() == tp.lower()), tp)
        missing = [c for c in ("n_precursors", "median_log2") if r[c] == ""]
        recs.append(dict(layer=layer, timepoint=tp_norm, fraction=fr, rep=rep, run=run,
                         n_str=r["n_precursors"], log2_str=r["median_log2"],
                         n=float(r["n_precursors"]) if r["n_precursors"] != "" else None,
                         lg=float(r["median_log2"]) if r["median_log2"] != "" else None,
                         missing=missing))
    df = pd.DataFrame(recs)

    # ---- cell 统计（同 layer x fraction 的全部 run，含非 CTRL）----
    cells, cell_stats, mad0 = [], {}, []
    missing_vals = []
    for layer in LAYERS:
        for fr in FRS:
            sub = df[(df.layer == layer) & (df.fraction == fr)]
            nv = [v for v in sub["n"] if v is not None]
            lv = [v for v in sub["lg"] if v is not None]
            if len(nv) != len(sub) or len(lv) != len(sub):
                missing_vals.append((layer, fr, len(sub) - len(nv), len(sub) - len(lv)))
            med_n, mad_n = robust(nv)
            med_l, mad_l = robust(lv)
            if mad_n == 0:
                mad0.append((layer, fr, "n_precursors"))
            if mad_l == 0:
                mad0.append((layer, fr, "median_log2"))
            cell_stats[(layer, fr)] = (med_n, mad_n, med_l, mad_l)
            cells.append(dict(layer=layer, fraction=fr, n_runs=len(sub),
                              median_n=fmt(med_n), mad_n=fmt(mad_n),
                              median_log2=fmt(med_l), mad_log2=fmt(mad_l)))
    cells_df = pd.DataFrame(cells, columns=CELL_COLS)

    # ---- 24 个 CTRL run ----
    rows, zraw = [], []
    ctrl = df[df.timepoint == "CTRL"].copy()
    ctrl["_li"] = ctrl.layer.map(LAYERS.index)
    ctrl["_fi"] = ctrl.fraction.map(FRS.index)
    ctrl["_ri"] = ctrl.rep.map(lambda x: REPS.index(x) if x in REPS else 99)
    ctrl = ctrl.sort_values(["_li", "_fi", "_ri"])
    for _, r in ctrl.iterrows():
        med_n, mad_n, med_l, mad_l = cell_stats[(r.layer, r.fraction)]
        zn = zscore(r.n, med_n, mad_n) if r.n is not None else None
        zl = zscore(r.lg, med_l, mad_l) if r.lg is not None else None
        zraw.append(dict(layer=r.layer, fraction=r.fraction, rep=r.rep, run=r.run,
                         n=r.n_str, lg=r.log2_str, zn=zn, zl=zl))
        rows.append(dict(layer=r.layer, fraction=r.fraction, rep=r.rep, run=r.run,
                         n_precursors=r.n_str, median_log2=r.log2_str,
                         cell_median_n=fmt(med_n), cell_mad_n=fmt(mad_n),
                         z_n_precursors=fmt(zn),
                         cell_median_log2=fmt(med_l), cell_mad_log2=fmt(mad_l),
                         z_median_log2=fmt(zl),
                         flag_n=flag(zn), flag_log2=flag(zl)))
    qc = pd.DataFrame(rows, columns=QC_COLS)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    qc.to_csv(OUT_TSV, sep="\t", index=False)
    cells_df.to_csv(OUT_CELLS, sep="\t", index=False)

    # ---- |z|>2 清单 ----
    hits = []
    for z in zraw:
        if z["zn"] is not None and abs(z["zn"]) > Z_THRESH:
            hits.append(dict(layer=z["layer"], fraction=z["fraction"], rep=z["rep"], run=z["run"],
                             metric="n_precursors", value=z["n"], z=fmt(z["zn"])))
        if z["zl"] is not None and abs(z["zl"]) > Z_THRESH:
            hits.append(dict(layer=z["layer"], fraction=z["fraction"], rep=z["rep"], run=z["run"],
                             metric="median_log2", value=z["lg"], z=fmt(z["zl"])))
    hits_df = pd.DataFrame(hits, columns=["layer", "fraction", "rep", "run", "metric", "value", "z"])

    # 与阈值距离 < 0.0005（4 位取整可能改变判定）的 z，记录
    near = [(z["run"], m, v) for z in zraw for m, v in (("n_precursors", z["zn"]), ("median_log2", z["zl"]))
            if v is not None and abs(abs(v) - Z_THRESH) < 0.0005]

    # ---- 计数 ----
    cnt = []
    for layer in LAYERS:
        q = qc[qc.layer == layer]
        n_n = int((q.flag_n == "yes").sum())
        n_l = int((q.flag_log2 == "yes").sum())
        n_any = int(((q.flag_n == "yes") | (q.flag_log2 == "yes")).sum())
        n_both = int(((q.flag_n == "yes") & (q.flag_log2 == "yes")).sum())
        n_entries = int((hits_df.layer == layer).sum())
        cnt.append(dict(layer=layer, n_ctrl_runs=len(q), n_flag_n=n_n, n_flag_log2=n_l,
                        n_flag_either=n_any, n_flag_both=n_both, n_entries=n_entries))
    cnt_df = pd.DataFrame(cnt)
    tot = dict(layer="合计", n_ctrl_runs=int(cnt_df.n_ctrl_runs.sum()),
               n_flag_n=int(cnt_df.n_flag_n.sum()), n_flag_log2=int(cnt_df.n_flag_log2.sum()),
               n_flag_either=int(cnt_df.n_flag_either.sum()), n_flag_both=int(cnt_df.n_flag_both.sum()),
               n_entries=int(cnt_df.n_entries.sum()))
    cnt_df = pd.concat([cnt_df, pd.DataFrame([tot])], ignore_index=True)

    # ---- 自检 ----
    checks = []
    checks.append(("输入 240 行", n_input == 240, f"{n_input}"))
    checks.append(("解析失败 0 行", len(parse_fail) == 0, f"{len(parse_fail)}"))
    n_ctrl_data = int((df.timepoint == "CTRL").sum())
    checks.append(("03_ctrl_qc.tsv 24 行（方案原文）", len(qc) == 24, f"{len(qc)}"))
    checks.append(("phospho / proteome CTRL 各 12（方案原文）", list(cnt_df.n_ctrl_runs[:2]) == [12, 12],
                   f"{list(cnt_df.n_ctrl_runs[:2])}"))
    checks.append(("03_ctrl_qc.tsv 行数 = 输入中 CTRL run 数 = 2x6x4", len(qc) == n_ctrl_data == 48,
                   f"{len(qc)} / {n_ctrl_data}"))
    checks.append(("每个 layer x fraction 的 CTRL run = 4（Rep1-4）", all(
        len(qc[(qc.layer == l) & (qc.fraction == f)]) == 4 for l in LAYERS for f in FRS), ""))
    checks.append(("03_ctrl_qc_cells.tsv 12 行", len(cells_df) == 12, f"{len(cells_df)}"))
    checks.append(("每个 cell n_runs = 20", bool((cells_df.n_runs == 20).all()),
                   f"{sorted(set(cells_df.n_runs))}"))
    checks.append(("每个 cell 5 tp x 4 rep 各 1", all(
        len(df[(df.layer == l) & (df.fraction == f) & (df.timepoint == t) & (df.rep == rp)]) == 1
        for l in LAYERS for f in FRS for t in TPS for rp in REPS), ""))
    raw_map = raw.set_index("run")
    same = all(raw_map.loc[r.run, "n_precursors"] == r.n_precursors and
               raw_map.loc[r.run, "median_log2"] == r.median_log2 for r in qc.itertuples())
    checks.append((f"{len(qc)} 个 CTRL run 的 n_precursors / median_log2 与输入原值一致", same, ""))
    checks.append(("|z|>2 条目共 25", len(hits_df) == 25, f"{len(hits_df)}"))
    n_ph = int((hits_df.layer == "phospho").sum())
    n_pr = int((hits_df.layer == "proteome").sum())
    checks.append(("phospho 9 条 / proteome 16 条", (n_ph, n_pr) == (9, 16), f"{n_ph} / {n_pr}"))
    ref = [("phospho", "FR2", "Rep4", -5.153), ("phospho", "FR6", "Rep4", 5.222),
           ("proteome", "FR3", "Rep4", 5.697), ("proteome", "FR6", "Rep4", 5.594)]
    for l, f, rp, zref in ref:
        zz = next(z["zl"] for z in zraw if (z["layer"], z["fraction"], z["rep"]) == (l, f, rp))
        checks.append((f"{l} {f} {rp} z_median_log2 ≈ {zref}",
                       zz is not None and abs(zz - zref) <= 0.001, fmt(zz)))
    checks.append(("MAD=0 的 (cell, 指标)", True, f"{len(mad0)}"))
    for name, ok, info in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {name}  {info}")

    # ---- md ----
    L = []
    L.append("# 03 CTRL 质检（run 级稳健 z）\n")
    L.append("脚本：`analysis/07_pilot2/03_ctrl_qc.py`（确定性，可复跑）。只报计数和数值。\n")
    L.append("## ① 输入\n")
    L.append(f"- `data/pilot/run_medians.tsv`：{n_input} 行（不含表头），列 `run, n_precursors, median_log2`；"
             f"读法 `pd.read_csv(sep='\\t', dtype=str, keep_default_na=False)`，空串 = 缺失（本表空串 0 个）。")
    L.append("- 解析：`layer` = `run` 第一个 `/` 前的前缀（phospho | proteome）；timepoint / fraction / rep 用 "
             "`code/Protein contour/Zhihan/Test/07_main_analysis.py` `parse_design`（`:32-38`）的三个正则 "
             "`_(2min|8min|20min|90min|CTRL)_`（re.I）、`_(FR\\d)_`、`_(Rep\\d)`。")
    L.append(f"- 解析失败：{len(parse_fail)} 行。layer 计数：" +
             "、".join(f"{l} {int((df.layer == l).sum())}" for l in LAYERS) +
             f"；CTRL run：{len(qc)}（" + "、".join(f"{l} {int((qc.layer == l).sum())}" for l in LAYERS) + "）。")
    L.append(f"- **与方案不符**：方案写 24 个 CTRL run（phospho 12 + proteome 12）；输入中 CTRL run 为 "
             f"{len(qc)} 个（2 layer × 6 fraction × 4 rep）。本表输出全部 {len(qc)} 个，未做子集挑选。"
             "方案验收参考值（25 条、phospho 9 / proteome 16、4 个 z）均在这 48 个 run 上复现，见自检。\n")
    L.append("## ② 定义\n")
    L.append("- 组（cell）= 同 layer × fraction 的 20 个 run（5 timepoint × 4 rep，**含非 CTRL**）；共 2 × 6 = 12 个 cell。")
    L.append("- 对 `n_precursors` 和 `median_log2` 分别在每个 cell 的 20 个值上算："
             "`med = median(20 值)`，`MAD = median(|x − med|)`，`z = (x − med) / (1.4826 × MAD)`；"
             "中位数用 `statistics.median`（偶数个取中间两值均值）。")
    L.append("- MAD = 0 时 z 留空（本数据 MAD = 0 的 (cell, 指标) 数见 ⑥）。")
    L.append(f"- 只输出 CTRL run（本数据 {len(qc)} 个：" + "、".join(f"{l} {int((qc.layer == l).sum())}" for l in LAYERS) +
             "）。`flag_n` / `flag_log2`：|z| > 2 → yes，否则 no（判定用未取整的 z；z 为空时 flag 留空）。")
    L.append("- 数值列保留 4 位小数；`n_precursors`、`median_log2` 两列为输入原字符串。\n")
    L.append("### cell 统计（`03_ctrl_qc_cells.tsv`，12 行）\n")
    L.append(md_table(cells_df, CELL_COLS) + "\n")
    L.append(f"## ③ {len(qc)} 个 CTRL run（`03_ctrl_qc.tsv`，{len(qc)} 行）\n")
    short = qc.copy()
    short["run"] = short["run"].map(lambda s: "`" + s + "`")
    L.append(md_table(short, QC_COLS) + "\n")
    L.append(f"## ④ |z| > 2 清单（{len(hits_df)} 条 (run, 指标)）\n")
    hd = hits_df.copy()
    hd.insert(0, "#", range(1, len(hd) + 1))
    hd["run"] = hd["run"].map(lambda s: "`" + s + "`")
    L.append(md_table(hd, ["#", "layer", "fraction", "rep", "run", "metric", "value", "z"]) + "\n")
    L.append("## ⑤ 计数\n")
    L.append("各 layer 的 CTRL run 中 |z| > 2 的 run 数（n = `n_precursors`，log2 = `median_log2`）：\n")
    L.append(md_table(cnt_df.astype(str), ["layer", "n_ctrl_runs", "n_flag_n", "n_flag_log2",
                                           "n_flag_either", "n_flag_both", "n_entries"]))
    L.append("\n列说明：`n_flag_n` = n 上 |z|>2 的 CTRL run 数；`n_flag_log2` = log2 上 |z|>2 的 CTRL run 数；"
             "`n_flag_either` = 任一指标；`n_flag_both` = 两指标同时；`n_entries` = (run, 指标) 条目数。\n")
    L.append("按 layer × fraction 的条目数：\n")
    pv = []
    for l in LAYERS:
        for f in FRS:
            h = hits_df[(hits_df.layer == l) & (hits_df.fraction == f)]
            pv.append(dict(layer=l, fraction=f,
                           n_entries_n=int((h.metric == "n_precursors").sum()),
                           n_entries_log2=int((h.metric == "median_log2").sum())))
    L.append(md_table(pd.DataFrame(pv).astype(str), ["layer", "fraction", "n_entries_n", "n_entries_log2"]) + "\n")
    L.append("## ⑥ 跳过项 / 记录\n")
    L.append(f"- MAD = 0 的 (cell, 指标)：{len(mad0)} 个" +
             ("" if not mad0 else "：" + "；".join(f"{a} {b} {c}" for a, b, c in mad0)) +
             "；因此 z 留空的格：" + str(int((qc.z_n_precursors == "").sum() + (qc.z_median_log2 == "").sum())) + " 个。")
    L.append(f"- 解析失败行：{len(parse_fail)}" + ("" if not parse_fail else "：" + "；".join(parse_fail)) + "。")
    L.append(f"- 输入缺失值（空串）：{len(missing_vals)} 个 cell 有缺失。")
    L.append(f"- |z| 与阈值 2 相差 < 0.0005（4 位取整可能影响判定）的格：{len(near)} 个"
             + ("" if not near else "：" + "；".join(f"{a} {b} {v:.6f}" for a, b, v in near)) + "。")
    L.append(f"- 方案「24 个 CTRL run / 24 行」与输入不符（实为 {len(qc)}），未按 24 截取；无其他跳过项"
             "（本任务只依赖 `run_medians.tsv`，已存在）。\n")
    L.append("## 自检\n")
    L.append("| 检查 | 结果 | 值 |\n|---|---|---|")
    for name, ok, info in checks:
        L.append(f"| {name} | {'PASS' if ok else 'FAIL'} | {info} |")
    L.append("")
    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    print(f"wrote {OUT_TSV}\nwrote {OUT_CELLS}\nwrote {OUT_MD}")
    print(hits_df.to_string(index=False))


if __name__ == "__main__":
    main()
