#!/usr/bin/env python3
"""02_events.py — 08_gate1 任务 2：GRB2 / SHC1 / CBL 三个易位事件（表 E-A..E-D）+ 图 + md。

按 analysis/08_gate1/00_design.md 的 §0、§1（作者 Movement Score 定义）、§2（共用定义）与「任务 2」执行。

输入（只读）:
  analysis/08_gate1/01_shares.tsv          区室份额（primary group）→ E-A、E-B、E-D、图
  analysis/08_gate1/01_groups.tsv          primary group = is_primary == "yes"
  data/pilot/gate1_proteins.tsv            fraction 级 quantity → E-C（Movement Score）
  analysis/07_pilot2/01_run_factors.tsv    layer=proteome 的 a（centered = quantity × 2^a）
  analysis/08_gate1/03_egfr_total.md       ⑦ 节：原样抄入其 ③④ 两节
  作者代码 translocation_plots.R           只用于核对 §1 引用行的原文
输出（analysis/08_gate1/）:
  02_events.tsv   长格式 table, gene, normalization, timepoint, key, value（E-A..E-D 全部数值）
  02_events.png   2 行（raw / centered）× 9 列（GRB2-Cyt … CBL-Nuc），每子图 4 条重复线
  02_events.md    ① 输入 ② 定义 ③ E-A ④ E-B ⑤ E-C ⑥ E-D ⑦ EGFR 总量 ⑧ 跳过项

用法:
  python3 analysis/08_gate1/02_events.py          （任意目录可复跑，输出确定）
  可选参数（仅供代码测试，默认即仓库路径）:
    --proteins PATH  --shares PATH  --groups PATH（默认与 --shares 同目录的 01_groups.tsv）
    --factors PATH   --egfr-md PATH  --outdir DIR
"""
import argparse
import math
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEF_PROT = REPO / "data" / "pilot" / "gate1_proteins.tsv"
DEF_SHARES = REPO / "analysis" / "08_gate1" / "01_shares.tsv"
DEF_FACT = REPO / "analysis" / "07_pilot2" / "01_run_factors.tsv"
DEF_EGFR_MD = REPO / "analysis" / "08_gate1" / "03_egfr_total.md"
DEF_OUT = REPO / "analysis" / "08_gate1"
R_FILE = (REPO / "code" / "Protein contour" / "Phospho" / "SpatialProteoDynamics.github.io-v1.0"
          / "SpatialProteoDynamics-SpatialProteoDynamics.github.io-a6b8aac" / "DataProcessing"
          / "translocation_plots.R")

GENES_ALL = ["GRB2", "SHC1", "CBL", "EGFR"]        # ① 测到 run 数 /120 四个都报
GENES = ["GRB2", "SHC1", "CBL"]                    # 任务 2 只做这三个
TPS = ["CTRL", "2min", "8min", "20min", "90min"]
EGF_TPS = ["2min", "8min", "20min", "90min"]
DIR_TPS = ["2min", "8min"]
FRS = [f"FR{i}" for i in range(1, 7)]
REPS = [f"Rep{i}" for i in range(1, 5)]
COMPS = [("Cyt", ("FR1", "FR2")), ("Mem", ("FR3", "FR4")), ("Nuc", ("FR5", "FR6"))]
CNAMES = [c for c, _ in COMPS]
NORMS = ["raw", "centered"]
N_DESIGN = len(TPS) * len(FRS) * len(REPS)          # 120
EXPECTED_TRIPS = {(t, f, r) for t in TPS for f in FRS for r in REPS}
MS_THRESHOLD = 0.1
TOL = 1e-9

IN_COLS = ["run", "genes", "protein_group", "quantity"]
FACT_COLS = ["layer", "run", "timepoint", "fraction", "rep", "a"]
GROUP_COLS = ["gene", "protein_group", "genes_raw", "n_runs_detected", "is_primary"]
SHARE_COLS = ["gene", "protein_group", "timepoint", "rep", "normalization",
              "n_fractions_present", "Cyt", "Mem", "Nuc"]
EA_COLS = ["gene", "normalization", "timepoint", "compartment", "n", "mean", "sd", "delta_vs_ctrl"]
EB_COLS = ["gene", "normalization", "timepoint", "n_mem_up_cyt_down", "n_mem_up", "n_cyt_down",
           "denominator"]
EC_COLS = ["gene", "normalization", "timepoint", "MS_author", "MS_max1", "MS_max2",
           "MS_perrep_mean", "ge_0.1"]
ED_COLS = ["gene", "normalization", "timepoint", "mem_mean_tp", "mem_mean_ctrl", "t", "p_welch",
           "n_tp", "n_ctrl"]
LONG_COLS = ["table", "gene", "normalization", "timepoint", "key", "value"]

NODATA = "缺数据：`gate1_proteins.tsv` 0 行，跳过"

# §1 引用的作者原文（00_design.md §1），运行时与 R 文件逐行核对
R_QUOTE = {
    19: 'prot_2min<-2^prot_2min',
    20: 'prot_ctrl<-2^prot_ctrl',
    23: 'FR<-c(rep("FR1",4), rep("FR2", 4),  rep("FR3", 4),  rep("FR4", 4),  rep("FR5", 4),  rep("FR6", 4))',
    24: 'prot_2min_mean<-t(apply(prot_2min, 1, function(x) tapply(x,FR,function(x) {mean(x)})))',
    25: 'prot_CTRL_mean<-t(apply(prot_ctrl, 1, function(x) tapply(x,FR,function(x) {mean(x)})))',
    30: 'prot_2min_scaled<-t(apply(prot_2min_mean_raw, 1, function(x) {x/sum(x)}))',
    31: 'prot_Ctrl_scaled<-t(apply(prot_CTRL_mean_raw, 1, function(x) {x/sum(x)}))',
    33: 'MS_2min<-abs(prot_2min_scaled-prot_Ctrl_scaled)',
    34: 'MS_2min_mean<-apply(MS_2min, 1, function(x) mean(x))',
    36: 'MS_max1<-(apply(MS_2min, 1, function(x) {match(max(x),x)}))',
    47: 'MS_max2<-(apply(MS_2min, 1, function(x) {match(maxN(x),x)}))',
    181: 'geom_text_repel(aes(label=Gene), data=temp_table[(temp_table$pval_combi_FDR<0.05) & (temp_table$MS_2min_mean>0.1),])+',
}

# 图：dataviz 参考调色板 categorical slot 1–4（固定顺序；validate_palette.js light/adjacent 通过，
# 对比度 WARN → 每个 rep 另加 marker 形状作二次编码，数值见 E-A 表 / 01_shares.tsv）
SURF, INK, INK2, MUTED, GRIDC = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
REP_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
REP_MARKERS = ["o", "s", "^", "D"]
FIG_DPI = 200


# ---- 照抄 code/Protein contour/Zhihan/Test/07_main_analysis.py parse_design (:32-38) ----
def parse_design(run):
    tp = re.search(r"_(2min|8min|20min|90min|CTRL)_", run, re.I)
    fr = re.search(r"_(FR\d)_", run)
    rep = re.search(r"_(Rep\d)", run)
    if not (tp and fr and rep):
        return None
    return (tp.group(1), fr.group(1), rep.group(1))
# ---- 照抄结束 ----


def parse_run(run):
    """返回 ((tp, fr, rep), None) 或 (None, 失败原因)。tp 统一成 TPS 写法（与任务 1 相同）。"""
    p = parse_design(run)
    if p is None:
        return None, "正则不匹配"
    tp, fr, rep = p
    tp = next(t for t in TPS if t.lower() == tp.lower())
    if fr not in FRS or rep not in REPS:
        return None, f"超出设计范围（{fr}, {rep}）"
    return (tp, fr, rep), None


def to_float(s):
    """空串 / 非数（含 nan、inf）→ None。"""
    s = str(s).strip()
    if s == "":
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def rel(p):
    p = Path(p).resolve()
    try:
        return str(p.relative_to(REPO))
    except ValueError:
        return str(p)


def isnum(x):
    return x is not None and not (isinstance(x, float) and math.isnan(x))


def f12(x):
    """tsv 数值：12 位有效数字；None/nan → 空串；整数原样。"""
    if x is None:
        return ""
    if isinstance(x, (int, np.integer)) and not isinstance(x, bool):
        return str(int(x))
    if isinstance(x, str):
        return x
    x = float(x)
    return "" if math.isnan(x) else format(x, ".12g")


def f4(s):
    return "" if s == "" else f"{float(s):.4f}"


def fsig(s, n=3):
    return "" if s == "" else format(float(s), f".{n}g")


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    for r in rows:
        out.append("| " + " | ".join(str(v) for v in r) + " |")
    return "\n".join(out)


def mean_or_none(v):
    return sum(v) / len(v) if v else None


def sd_or_none(v):
    if len(v) < 2:
        return None
    m = sum(v) / len(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))


def author_ms(m_tp, m_ctrl):
    """§1：s = m / Σ m（6 个 fraction）；MS = mean |s_tp − s_CTRL|；
    MS_max1 = match(max(x), x)；MS_max2 = match(maxN(x), x)，maxN = 第二大的值（R :38-45）。
    m_* 为长度 6 的 list（FR1..FR6），任一为 None 或和为 0 → 返回 None。"""
    if any(v is None for v in m_tp + m_ctrl):
        return None
    st, sc = sum(m_tp), sum(m_ctrl)
    if st == 0 or sc == 0:
        return None
    d = [abs(a / st - b / sc) for a, b in zip(m_tp, m_ctrl)]
    ms = sum(d) / len(d)
    i1 = d.index(max(d))                    # R match()：第一个等于该值的位置
    i2 = d.index(sorted(d)[len(d) - 2])     # R maxN(x, 2) = sort(x)[len-1]
    return ms, FRS[i1], FRS[i2], d


def main():
    ap = argparse.ArgumentParser(description="08_gate1 任务 2：三个易位事件")
    ap.add_argument("--proteins", default=str(DEF_PROT))
    ap.add_argument("--shares", default=str(DEF_SHARES))
    ap.add_argument("--groups", default=None, help="默认：与 --shares 同目录的 01_groups.tsv")
    ap.add_argument("--factors", default=str(DEF_FACT))
    ap.add_argument("--egfr-md", default=str(DEF_EGFR_MD))
    ap.add_argument("--outdir", default=str(DEF_OUT))
    args = ap.parse_args()
    p_prot, p_shares, p_fact = Path(args.proteins), Path(args.shares), Path(args.factors)
    p_groups = Path(args.groups) if args.groups else p_shares.parent / "01_groups.tsv"
    p_egfr = Path(args.egfr_md)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    o_tsv, o_png, o_md = outdir / "02_events.tsv", outdir / "02_events.png", outdir / "02_events.md"

    checks = []

    def chk(name, ok, val=""):
        status = ok if isinstance(ok, str) else ("PASS" if ok else "FAIL")
        checks.append((name, status, val))

    # ================================================================ 读表
    P = read(p_prot)
    miss = [c for c in IN_COLS if c not in P.columns]
    if miss:
        sys.exit(f"输入缺列 {miss}：{p_prot}")
    n_input = len(P)

    F = read(p_fact)
    miss = [c for c in FACT_COLS if c not in F.columns]
    if miss:
        sys.exit(f"run_factors 缺列 {miss}：{p_fact}")
    FP = F[F["layer"] == "proteome"]
    a_map, fact_bad, fact_dup = {}, 0, 0
    for r in FP.itertuples(index=False):
        key = (r.timepoint, r.fraction, r.rep)
        a = to_float(r.a)
        if key not in EXPECTED_TRIPS or a is None:
            fact_bad += 1
            continue
        if key in a_map:
            fact_dup += 1
            continue
        a_map[key] = a

    G = read(p_groups)
    miss = [c for c in GROUP_COLS if c not in G.columns]
    if miss:
        sys.exit(f"01_groups.tsv 缺列 {miss}：{p_groups}")
    S = read(p_shares)
    miss = [c for c in SHARE_COLS if c not in S.columns]
    if miss:
        sys.exit(f"01_shares.tsv 缺列 {miss}：{p_shares}")
    chk("输入 01_shares.tsv / 01_groups.tsv 列 = 任务 1 方案列",
        list(S.columns) == SHARE_COLS and list(G.columns) == GROUP_COLS)
    n_shares = len(S)

    # primary group：01_groups.tsv 中 is_primary == yes
    primary, multi_primary = {}, []
    for g in GENES_ALL:
        ys = sorted(G.loc[(G.gene == g) & (G.is_primary == "yes"), "protein_group"])
        if len(ys) > 1:
            multi_primary.append((g, ys))
        if ys:
            primary[g] = ys[0]
    chk("01_groups.tsv：每个基因至多一个 is_primary=yes", not multi_primary,
        "; ".join(f"{g}: {', '.join(x)}" for g, x in multi_primary))

    # ================================================================ fraction 级数值（规则同任务 1）
    runs = sorted(set(P["run"]))
    run_trip, parse_fail = {}, []
    for run in runs:
        trip, why = parse_run(run)
        if trip is None:
            parse_fail.append((run, why))
        else:
            run_trip[run] = trip
    trip_runs = defaultdict(list)
    for run, trip in run_trip.items():
        trip_runs[trip].append(run)
    n_dup_trips = sum(1 for v in trip_runs.values() if len(v) > 1)

    E = P.reset_index(drop=False).rename(columns={"index": "x_row"})
    E = E.assign(x_tok=E["genes"].str.split(";")).explode("x_tok")
    E["x_tok"] = E["x_tok"].fillna("").astype(str).str.strip()
    M = E[E["x_tok"].isin(GENES_ALL)].drop_duplicates(["x_row", "x_tok"])

    cell_vals = defaultdict(list)
    for r in M.itertuples(index=False):
        trip = run_trip.get(r.run)
        if trip is None:
            continue
        cell_vals[(r.x_tok, r.protein_group) + trip].append(to_float(r.quantity))
    cell_q, n_conflicts = {}, 0
    gene_detected = defaultdict(set)
    for key, lst in cell_vals.items():
        qs = [q for q in lst if q is not None]
        if not qs:
            continue
        gene_detected[key[0]].add(key[2:])
        if len(set(qs)) == 1:
            cell_q[key] = qs[0]
        else:
            n_conflicts += 1

    def qv(g, tp, fr, rep, norm):
        """primary group 的 fraction 级数值（raw 或 centered）；缺 → None。"""
        grp = primary.get(g)
        if grp is None:
            return None
        q = cell_q.get((g, grp, tp, fr, rep))
        if q is None:
            return None
        if norm == "raw":
            return q
        a = a_map.get((tp, fr, rep))
        return None if a is None else q * 2.0 ** a

    # ================================================================ 份额（来自 01_shares.tsv）
    Sg = S[S.gene.isin(GENES)]
    share = {}                        # (g, norm, tp, rep) -> {Cyt, Mem, Nuc: float|None}
    bad_group_rows = 0
    for r in Sg.itertuples(index=False):
        if r.protein_group != primary.get(r.gene):
            bad_group_rows += 1
        share[(r.gene, r.normalization, r.timepoint, r.rep)] = {c: to_float(getattr(r, c)) for c in CNAMES}
    if n_shares == 0:
        chk("01_shares.tsv 的 protein_group = 01_groups.tsv 的 primary", "不适用", "0 行")
    else:
        chk("01_shares.tsv 的 protein_group = 01_groups.tsv 的 primary", bad_group_rows == 0,
            f"不一致 {bad_group_rows} 行")

    # 一致性：用 fraction 级数值按 §2 重算区室份额，应与 01_shares.tsv 相同
    recomputed, diff_max, n_cmp, n_missing_side = {}, 0.0, 0, 0
    for g in GENES:
        for tp in TPS:
            for rep in REPS:
                for norm in NORMS:
                    vals = {fr: qv(g, tp, fr, rep, norm) for fr in FRS}
                    pres = {f: v for f, v in vals.items() if v is not None}
                    if not pres:
                        continue
                    tot = sum(pres.values())
                    recomputed[(g, norm, tp, rep)] = {
                        c: (None if tot == 0 or not any(f in pres for f in frs)
                            else sum(pres[f] for f in frs if f in pres) / tot) for c, frs in COMPS}
    keys_union = set(recomputed) | set(share)
    for k in keys_union:
        a, b = recomputed.get(k), share.get(k)
        if a is None or b is None:
            n_missing_side += 1
            continue
        for c in CNAMES:
            if (a[c] is None) != (b[c] is None):
                n_missing_side += 1
            elif a[c] is not None:
                n_cmp += 1
                diff_max = max(diff_max, abs(a[c] - b[c]))
    if n_input == 0 and n_shares == 0:
        chk("01_shares.tsv 与 fraction 级重算份额一致（GRB2/SHC1/CBL）", "不适用", "0 行")
    else:
        chk("01_shares.tsv 与 fraction 级重算份额一致（GRB2/SHC1/CBL）",
            n_missing_side == 0 and diff_max <= TOL,
            f"比较 {n_cmp} 格，max|差| = {diff_max:.1e}，一侧缺 {n_missing_side}")

    have_shares = n_shares > 0
    have_prot = n_input > 0

    # ================================================================ E-A
    EA = []
    if have_shares:
        for g in GENES:
            for norm in NORMS:
                stat = {}
                for tp in TPS:
                    for c in CNAMES:
                        v = [share[(g, norm, tp, rep)][c] for rep in REPS
                             if (g, norm, tp, rep) in share and share[(g, norm, tp, rep)][c] is not None]
                        stat[(tp, c)] = (len(v), mean_or_none(v), sd_or_none(v))
                for tp in TPS:
                    for c in CNAMES:
                        n, m, sd = stat[(tp, c)]
                        mc = stat[("CTRL", c)][1]
                        d = None if (m is None or mc is None) else m - mc
                        EA.append([g, norm, tp, c, n, m, sd, d])

    # ================================================================ E-B
    EB, eb_short = [], []
    if have_shares:
        for g in GENES:
            for norm in NORMS:
                for tp in DIR_TPS:
                    both = up = down = den = 0
                    for rep in REPS:
                        a, b = share.get((g, norm, tp, rep)), share.get((g, norm, "CTRL", rep))
                        if a is None or b is None or None in (a["Mem"], a["Cyt"], b["Mem"], b["Cyt"]):
                            continue
                        den += 1
                        dm, dc = a["Mem"] - b["Mem"], a["Cyt"] - b["Cyt"]
                        up += dm > 0
                        down += dc < 0
                        both += (dm > 0) and (dc < 0)
                    if den < len(REPS):
                        eb_short.append((g, norm, tp, den))
                    fmt = (lambda n: f"{n}/{den}") if den > 0 else (lambda n: "")
                    EB.append([g, norm, tp, fmt(both), fmt(up), fmt(down), den])

    # ================================================================ E-C
    EC, ec_partial_mean, ec_missing_fr, ec_perrep_short, ec_dvec = [], [], [], [], {}
    if have_prot:
        for g in GENES:
            for norm in NORMS:
                fm = {}                                        # (tp, fr) -> (mean, n_reps)
                for tp in TPS:
                    for fr in FRS:
                        v = [x for x in (qv(g, tp, fr, rep, norm) for rep in REPS) if x is not None]
                        fm[(tp, fr)] = (mean_or_none(v), len(v))
                        if 0 < len(v) < len(REPS):
                            ec_partial_mean.append((g, norm, tp, fr, len(v)))
                for tp in EGF_TPS:
                    m_tp = [fm[(tp, fr)][0] for fr in FRS]
                    m_ct = [fm[("CTRL", fr)][0] for fr in FRS]
                    res = author_ms(m_tp, m_ct)
                    if res is None:
                        ms = mx1 = mx2 = None
                        miss = [f"{t}/{fr}" for t in ("CTRL", tp) for fr in FRS if fm[(t, fr)][0] is None]
                        ec_missing_fr.append((g, norm, tp, ", ".join(miss) if miss else "Σ 均值 = 0"))
                    else:
                        ms, mx1, mx2, dv = res
                        ec_dvec[(g, norm, tp)] = dv
                    per = []
                    for rep in REPS:
                        qt = [qv(g, tp, fr, rep, norm) for fr in FRS]
                        qc = [qv(g, "CTRL", fr, rep, norm) for fr in FRS]
                        r2 = author_ms(qt, qc)
                        if r2 is not None:
                            per.append(r2[0])
                    if len(per) < len(REPS):
                        ec_perrep_short.append((g, norm, tp, len(per)))
                    pm = mean_or_none(per)
                    ge = "" if ms is None else ("yes" if ms >= MS_THRESHOLD else "no")
                    EC.append([g, norm, tp, ms, mx1, mx2, pm, ge])

    # ================================================================ E-D
    ED = []
    if have_shares:
        for g in GENES:
            for norm in NORMS:
                vc = [share[(g, norm, "CTRL", r)]["Mem"] for r in REPS
                      if (g, norm, "CTRL", r) in share and share[(g, norm, "CTRL", r)]["Mem"] is not None]
                for tp in EGF_TPS:
                    vt = [share[(g, norm, tp, r)]["Mem"] for r in REPS
                          if (g, norm, tp, r) in share and share[(g, norm, tp, r)]["Mem"] is not None]
                    t = p = None
                    if len(vt) >= 2 and len(vc) >= 2:
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            res = stats.ttest_ind(vt, vc, equal_var=False)
                        t, p = float(res.statistic), float(res.pvalue)
                        if math.isnan(t) or math.isnan(p):
                            t = p = None
                    ED.append([g, norm, tp, mean_or_none(vt), mean_or_none(vc), t, p, len(vt), len(vc)])

    # ================================================================ 写 tsv（长格式）
    LONG = []
    for r in EA:
        g, norm, tp, c, n, m, sd, d = r
        for k, v in (("n", n), ("mean", m), ("sd", sd), ("delta_vs_ctrl", d)):
            LONG.append(["E-A", g, norm, tp, f"{c}_{k}", f12(v)])
    for r in EB:
        for k, v in zip(EB_COLS[3:], r[3:]):
            LONG.append(["E-B", r[0], r[1], r[2], k, f12(v)])
    for r in EC:
        for k, v in zip(EC_COLS[3:], r[3:]):
            LONG.append(["E-C", r[0], r[1], r[2], k, f12(v)])
    for r in ED:
        for k, v in zip(ED_COLS[3:], r[3:]):
            LONG.append(["E-D", r[0], r[1], r[2], k, f12(v)])
    pd.DataFrame(LONG, columns=LONG_COLS).to_csv(o_tsv, sep="\t", index=False)

    # ================================================================ 图
    fig = plt.figure(figsize=(18.5, 5.6), dpi=FIG_DPI, facecolor=SURF)
    wr = [1, 1, 1, 0.2, 1, 1, 1, 0.2, 1, 1, 1]
    gs = GridSpec(2, len(wr), figure=fig, width_ratios=wr, left=0.055, right=0.995,
                  top=0.80, bottom=0.10, wspace=0.10, hspace=0.30)
    col_pos = [0, 1, 2, 4, 5, 6, 8, 9, 10]
    panels = [(g, c) for g in GENES for c in CNAMES]
    xs = list(range(len(TPS)))
    panel_lines, panel_nodata = {}, []
    for i, norm in enumerate(NORMS):
        for j, (g, c) in enumerate(panels):
            ax = fig.add_subplot(gs[i, col_pos[j]])
            ax.set_facecolor(SURF)
            ax.set_xlim(-0.4, len(TPS) - 0.6)
            ax.set_ylim(0, 1)
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
            ax.set_xticks(xs)
            ax.grid(axis="y", color=GRIDC, linewidth=0.6, linestyle="-")
            ax.set_axisbelow(True)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            for s in ("left", "bottom"):
                ax.spines[s].set_color(MUTED)
                ax.spines[s].set_linewidth(0.6)
            ax.tick_params(colors=MUTED, labelcolor=INK2, labelsize=7, length=2.5, width=0.6)
            ax.set_xticklabels(TPS if i == len(NORMS) - 1 else [], fontsize=6.6)
            if j % 3 != 0:
                ax.set_yticklabels([])
            if j == 0:
                ax.set_ylabel(f"{norm}\nshare", color=INK, fontsize=9)
            if i == 0:
                ax.set_title(f"{g}-{c}", color=INK, fontsize=9.5, pad=5)
            has = any((g, norm, tp, rep) in share and share[(g, norm, tp, rep)][c] is not None
                      for tp in TPS for rep in REPS)
            if not has:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center",
                        color=MUTED, fontsize=9)
                panel_nodata.append((norm, g, c))
                panel_lines[(norm, g, c)] = 0
                continue
            n_lines = 0
            for k, rep in enumerate(REPS):
                ys = [share[(g, norm, tp, rep)][c]
                      if (g, norm, tp, rep) in share and share[(g, norm, tp, rep)][c] is not None
                      else np.nan for tp in TPS]
                ax.plot(xs, ys, color=REP_COLORS[k], linewidth=1.5, marker=REP_MARKERS[k],
                        markersize=4.5, markeredgecolor=SURF, markeredgewidth=0.8,
                        solid_joinstyle="round", solid_capstyle="round", label=rep, zorder=3 + k)
                n_lines += 1
            panel_lines[(norm, g, c)] = n_lines
    handles = [Line2D([0], [0], color=REP_COLORS[k], linewidth=1.5, marker=REP_MARKERS[k],
                      markersize=5, markeredgecolor=SURF, markeredgewidth=0.8, label=rep)
               for k, rep in enumerate(REPS)]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 0.985), ncol=4,
               frameon=False, fontsize=9, labelcolor=INK2, handlelength=2.2, columnspacing=1.4)
    fig.text(0.055, 0.955, "Compartment share per replicate: GRB2, SHC1, CBL (primary protein group)",
             color=INK, fontsize=12, fontweight="semibold", ha="left", va="center")
    fig.text(0.055, 0.905, "Source 01_shares.tsv. Rows: normalization (raw, centered). "
             "Y: share of summed fractions (Cyt = FR1+FR2, Mem = FR3+FR4, Nuc = FR5+FR6), 0 to 1. "
             "X: timepoint, equally spaced. One line per replicate; gaps = missing share.",
             color=INK2, fontsize=8.2, ha="left", va="center")
    n_axes = len(fig.axes)
    fig.savefig(o_png, dpi=FIG_DPI, facecolor=SURF)
    plt.close(fig)

    # ================================================================ ⑦ EGFR 节
    egfr_block, egfr_ok = None, False
    if p_egfr.exists():
        lines = p_egfr.read_text(encoding="utf-8").splitlines()
        heads = [i for i, ln in enumerate(lines) if ln.startswith("## ")]
        blocks = []
        for h in heads:
            if lines[h].startswith("## ③") or lines[h].startswith("## ④"):
                nxt = next((x for x in heads if x > h), len(lines))
                blk = lines[h:nxt]
                while blk and blk[-1].strip() == "":
                    blk = blk[:-1]
                blocks.append(blk)
        if len(blocks) == 2:
            orig = blocks[0] + [""] + blocks[1]
            egfr_block = [("#" + ln) if ln.startswith("#") else ln for ln in orig]
            back = [ln[1:] if ln.startswith("##") else ln for ln in egfr_block]
            egfr_ok = back == orig
    chk("⑦ 节 = 03_egfr_total.md 的 ③④ 原文（仅标题降一级）",
        egfr_ok if egfr_block is not None else False,
        rel(p_egfr) if egfr_block is not None else f"未找到 ③④：{rel(p_egfr)}")

    # ================================================================ 自检
    # R 原文核对
    if R_FILE.exists():
        rl = R_FILE.read_text(encoding="utf-8", errors="replace").splitlines()
        bad = [n for n, s in R_QUOTE.items() if n > len(rl) or rl[n - 1].strip() != s]
        chk("§1 引用的 R 原文与 translocation_plots.R 行号逐行一致", not bad,
            f"{len(R_QUOTE)} 行" + (f"；不一致行 {bad}" if bad else ""))
    else:
        chk("§1 引用的 R 原文与 translocation_plots.R 行号逐行一致", "FAIL", "R 文件不存在")

    R_L = read(o_tsv)
    chk("02_events.tsv 列 = table, gene, normalization, timepoint, key, value", list(R_L.columns) == LONG_COLS)
    exp_rows = {"E-A": 90, "E-B": 12, "E-C": 24, "E-D": 24}
    got_rows = {"E-A": len(EA), "E-B": len(EB), "E-C": len(EC), "E-D": len(ED)}
    if not have_prot and not have_shares:
        chk("0 行输入：E-A..E-D 全为 0 行，02_events.tsv 只有表头",
            all(v == 0 for v in got_rows.values()) and len(R_L) == 0,
            " / ".join(f"{k} {v}" for k, v in got_rows.items()) + f"；tsv {len(R_L)} 行")
    else:
        chk("表行数 E-A 90 / E-B 12 / E-C 24 / E-D 24", got_rows == exp_rows,
            " / ".join(f"{k} {v}" for k, v in got_rows.items()))
    # 长表 ↔ 宽表
    n_long_exp = 4 * len(EA) + 4 * len(EB) + 5 * len(EC) + 6 * len(ED)
    lk = {(r.table, r.gene, r.normalization, r.timepoint, r.key): r.value for r in R_L.itertuples(index=False)}
    rt_bad = sum(1 for row in LONG if lk.get(tuple(row[:5])) != row[5])
    chk("02_events.tsv 覆盖 E-A..E-D 全部数值（行数 = 4·EA + 4·EB + 5·EC + 6·ED，逐格回读一致）",
        len(R_L) == n_long_exp and len(lk) == len(R_L) and rt_bad == 0,
        f"{len(R_L)} / {n_long_exp}，不一致 {rt_bad}")

    # E-A：mean/sd 由回读的 01_shares.tsv 独立复算（pandas groupby）
    if EA:
        S2 = read(p_shares)
        S2 = S2[S2.gene.isin(GENES)]
        L2 = S2.melt(id_vars=["gene", "normalization", "timepoint", "rep"], value_vars=CNAMES,
                     var_name="compartment", value_name="v")
        L2 = L2[L2.v != ""].assign(v=lambda d: d.v.astype(float))
        agg = L2.groupby(["gene", "normalization", "timepoint", "compartment"]).v.agg(["count", "mean", "std"])
        dm = dsd = 0.0
        bad_n = 0
        for g, norm, tp, c, n, m, sd, d in EA:
            if (g, norm, tp, c) in agg.index:
                a = agg.loc[(g, norm, tp, c)]
                bad_n += int(a["count"]) != n
                dm = max(dm, abs(a["mean"] - float(lk[("E-A", g, norm, tp, f"{c}_mean")])))
                if n >= 2:
                    dsd = max(dsd, abs(a["std"] - float(lk[("E-A", g, norm, tp, f"{c}_sd")])))
            else:
                bad_n += n != 0
        chk("E-A：n / mean / sd(ddof=1) 由 01_shares.tsv 复算一致", bad_n == 0 and dm <= TOL and dsd <= TOL,
            f"n 不一致 {bad_n}，max|Δmean| = {dm:.1e}，max|Δsd| = {dsd:.1e}")
        dd = max((abs(float(lk[("E-A", g, norm, tp, f"{c}_delta_vs_ctrl")])
                      - (float(lk[("E-A", g, norm, tp, f"{c}_mean")])
                         - float(lk[("E-A", g, norm, "CTRL", f"{c}_mean")])))
                  for g, norm, tp, c, n, m, sd, d in EA if d is not None), default=0.0)
        chk("E-A：delta_vs_ctrl = mean_tp − mean_CTRL", dd <= TOL, f"max|差| = {dd:.1e}")
    else:
        chk("E-A：n / mean / sd(ddof=1) 由 01_shares.tsv 复算一致", "不适用", "0 行")

    # E-B：分母 = 两边 Mem、Cyt 都有值的 rep 数（pandas merge 独立复算）；计数 ≤ 分母
    if EB:
        S2 = read(p_shares)
        S2 = S2[S2.gene.isin(GENES)]
        ok_rows = 0
        for g, norm, tp, nb, nu, nd, den in EB:
            a = S2[(S2.gene == g) & (S2.normalization == norm) & (S2.timepoint == tp)]
            b = S2[(S2.gene == g) & (S2.normalization == norm) & (S2.timepoint == "CTRL")]
            mm = a.merge(b, on="rep", suffixes=("_t", "_c"))
            mm = mm[(mm.Mem_t != "") & (mm.Cyt_t != "") & (mm.Mem_c != "") & (mm.Cyt_c != "")]
            dmm = mm.Mem_t.astype(float) - mm.Mem_c.astype(float)
            dcc = mm.Cyt_t.astype(float) - mm.Cyt_c.astype(float)
            exp = [f"{int(x)}/{len(mm)}" if len(mm) else "" for x in
                   (((dmm > 0) & (dcc < 0)).sum(), (dmm > 0).sum(), (dcc < 0).sum())]
            ok_rows += (len(mm) == den and exp == [nb, nu, nd])
        chk("E-B：分母 = 01_shares 中两边 Mem、Cyt 都有值的 rep 数；三个计数复算一致", ok_rows == len(EB),
            f"{ok_rows}/{len(EB)} 行一致")
    else:
        chk("E-B：分母 = 01_shares 中两边 Mem、Cyt 都有值的 rep 数", "不适用", "0 行")

    # E-C：MS 由 fraction 级（6 个 fraction）重复均值独立复算（numpy 数组路径），max 为 FR 编号
    if EC:
        bad_ms = bad_max = n_fr_used_bad = 0
        n_ms = 0
        for g, norm, tp, ms, mx1, mx2, pm, ge in EC:
            def arr(t):
                return np.array([[np.nan if qv(g, t, fr, rep, norm) is None else qv(g, t, fr, rep, norm)
                                  for fr in FRS] for rep in REPS], dtype=float)
            At, Ac = arr(tp), arr("CTRL")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mt, mc = np.nanmean(At, axis=0), np.nanmean(Ac, axis=0)
            if ms is None:
                continue
            n_ms += 1
            n_fr_used_bad += (mt.shape[0] != 6 or np.isnan(mt).any() or np.isnan(mc).any())
            dv = np.abs(mt / mt.sum() - mc / mc.sum())
            bad_ms += abs(dv.mean() - ms) > TOL
            order = np.argsort(-dv, kind="stable")
            bad_max += (mx1 != FRS[int(order[0])]) or (mx2 not in FRS) or (mx1 not in FRS)
            if not (np.sort(dv)[-1] == np.sort(dv)[-2]):
                bad_max += mx2 != FRS[int(order[1])]
        chk("E-C：MS_author 由 fraction 级（6 个 fraction）重复均值复算一致", bad_ms == 0 and n_fr_used_bad == 0,
            f"{n_ms} 行有值，不一致 {bad_ms}")
        chk("E-C：MS_max1 / MS_max2 ∈ FR1..FR6，且为 |Δs| 最大 / 次大的 fraction", bad_max == 0,
            f"不一致 {bad_max}")
        # 区室份额版 MS（3 区室）不是本表的值：仅确认表值不等于它（防止误用区室份额）
        same_as_comp = 0
        for g, norm, tp, ms, *_ in EC:
            if ms is None:
                continue
            ct = [mean_or_none([share[(g, norm, tp, r)][c] for r in REPS
                                if (g, norm, tp, r) in share and share[(g, norm, tp, r)][c] is not None])
                  for c in CNAMES]
            cc = [mean_or_none([share[(g, norm, "CTRL", r)][c] for r in REPS
                                if (g, norm, "CTRL", r) in share and share[(g, norm, "CTRL", r)][c] is not None])
                  for c in CNAMES]
            if None not in ct + cc:
                same_as_comp += abs(sum(abs(x - y) for x, y in zip(ct, cc)) / 3 - ms) <= TOL
        chk("E-C：MS_author 不等于用 3 个区室份额算的同式值", same_as_comp == 0, f"相等 {same_as_comp}")
        bad_ge = sum(1 for r in EC if r[3] is not None and (r[7] == "yes") != (r[3] >= MS_THRESHOLD))
        chk("E-C：ge_0.1 = (MS_author ≥ 0.1)", bad_ge == 0, f"不一致 {bad_ge}")
    else:
        chk("E-C：MS_author 由 fraction 级（6 个 fraction）重复均值复算一致", "不适用", "0 行")

    # E-D：Welch 手算（公式）与 scipy 一致
    if ED:
        bad_t = 0
        for g, norm, tp, mt_, mc_, t, p, nt, nc in ED:
            if t is None:
                continue
            vt = [share[(g, norm, tp, r)]["Mem"] for r in REPS
                  if (g, norm, tp, r) in share and share[(g, norm, tp, r)]["Mem"] is not None]
            vc = [share[(g, norm, "CTRL", r)]["Mem"] for r in REPS
                  if (g, norm, "CTRL", r) in share and share[(g, norm, "CTRL", r)]["Mem"] is not None]
            v1, v2 = np.var(vt, ddof=1) / len(vt), np.var(vc, ddof=1) / len(vc)
            tm = (np.mean(vt) - np.mean(vc)) / math.sqrt(v1 + v2)
            df = (v1 + v2) ** 2 / (v1 ** 2 / (len(vt) - 1) + v2 ** 2 / (len(vc) - 1))
            pm_ = 2 * stats.t.sf(abs(tm), df)
            bad_t += abs(tm - t) > 1e-9 * max(1, abs(t)) or abs(pm_ - p) > 1e-9
        chk("E-D：t / p_welch 与 Welch 公式手算一致（双侧）", bad_t == 0, f"不一致 {bad_t}")
    else:
        chk("E-D：t / p_welch 与 Welch 公式手算一致（双侧）", "不适用", "0 行")

    # 图
    try:
        from PIL import Image
        with Image.open(o_png) as im:
            png_dpi = im.info.get("dpi", (0, 0))[0]
            png_size = im.size
    except Exception as e:  # pragma: no cover
        png_dpi, png_size = 0, (0, 0)
        print("PIL 读图失败：", e)
    data_panels = [k for k, v in panel_lines.items() if (k[0], k[1], k[2]) not in panel_nodata]
    chk("02_events.png：子图数 = 18（2 行 × 9 列），dpi ≥ 150", n_axes == 18 and round(png_dpi) >= 150,
        f"子图 {n_axes}，dpi {round(png_dpi)}，{png_size[0]}×{png_size[1]} px")
    if data_panels:
        chk("02_events.png：有数据的子图每图 4 条重复线", all(panel_lines[k] == 4 for k in data_panels),
            f"有数据子图 {len(data_panels)}，'no data' 子图 {len(panel_nodata)}")
    else:
        chk("02_events.png：无数据时 18 个子图均写 'no data'", len(panel_nodata) == 18,
            f"'no data' 子图 {len(panel_nodata)}")

    # ================================================================ md
    L = []
    L.append("# 08_gate1 任务 2：GRB2 / SHC1 / CBL 三个易位事件\n")
    L.append("生成脚本：`analysis/08_gate1/02_events.py`（按 `analysis/08_gate1/00_design.md` §0、§1、§2、任务 2 执行；"
             "`python3`，任意目录可复跑，输出确定）。只报数和方向计数。\n")
    L.append("输出：`02_events.tsv`（长格式，E-A..E-D 全部数值）、`02_events.png`、本文件。\n")

    # ① 输入
    L.append("## ① 输入\n")
    L.append(f"- `{rel(p_prot)}`：数据行 {n_input}（不含表头）；distinct run {len(runs)}；解析成功 "
             f"{len(run_trip)}/{N_DESIGN}；解析失败 {len(parse_fail)}；多 run 同三元组 {n_dup_trips}；"
             f"值冲突格 {n_conflicts}。")
    L.append(f"- `{rel(p_fact)}`：layer=proteome {len(FP)} 行，可用 (timepoint, fraction, rep) → a "
             f"{len(a_map)} 个" + (f"；无效 {fact_bad}、重复 {fact_dup}" if fact_bad or fact_dup else "") + "。")
    L.append(f"- `{rel(p_groups)}`：{len(G)} 行；primary（is_primary=yes）："
             + ("、".join(f"{g} `{primary[g]}`" for g in GENES_ALL if g in primary) or "无") + "。")
    L.append(f"- `{rel(p_shares)}`：{n_shares} 行；其中 "
             + "、".join(f"{g} {int((S.gene == g).sum())}" for g in GENES) + " 行。")
    L.append(f"- `{rel(p_egfr)}`：⑦ 节来源。")
    L.append("- **四个蛋白各在多少 run 测到（/120，任一 protein_group 有有效 quantity 的三元组数，口径同任务 1）**："
             + "，".join(f"{g} {len(gene_detected[g])}/{N_DESIGN}" for g in GENES_ALL) + "。\n")

    # ② 定义
    L.append("## ② 定义\n")
    L.append("### 作者 Movement Score（`00_design.md` §1，原文与行号）\n")
    L.append(f"`{rel(R_FILE)}`（运行时逐行核对，见 ⑧ 自检）：\n")
    L.append("```r")
    for n, s in R_QUOTE.items():
        L.append(f"{n:<3} {'  ' if n == 181 else ''}{s}")
    L.append("```")
    L.append("（:181 行尾的 `+` 是 ggplot 链式续行符，`00_design.md` §1 的引用省略了它；其余 11 行与 §1 引用逐字相同。）\n")
    L.append("即：① log2 值还原为线性（:19-20）；② 每个 fraction 内 4 个重复取**均值**（:23-25），得 6 个 fraction 均值；"
             "③ 份额 = 该 fraction 均值 / 6 个均值之和（:30-31）；④ 每个 fraction 的 |份额_tp − 份额_CTRL|（:33）；"
             "⑤ **Movement Score = 这 6 个绝对差的均值**（:34）；⑥ 差最大和次大的 fraction 记为 MS_max1 / MS_max2（:36, :47）；"
             "⑦ 图上标注阈值 `MS_2min_mean > 0.1` 且 `pval_combi_FDR < 0.05`（:181）。作者只对 2 min vs CTRL 算。\n")
    L.append("### 共用定义（`00_design.md` §2，照抄任务 2 用到的条目）\n")
    L.append("- 份额：对每个 (gene, protein_group, timepoint, rep)：`n_fractions_present` = 6 个 fraction 中有值的个数；"
             "`share_f = q_f / Σ_{present} q_f`；`Cyt = share_FR1 + share_FR2`、`Mem = share_FR3 + share_FR4`、"
             "`Nuc = share_FR5 + share_FR6`；某区室两个 fraction 都缺 → 该区室格留空；不填补。（取自 `01_shares.tsv`，primary group。）")
    L.append("- 区室均值/标准差：按 timepoint 对 4 个重复的份额取 `mean`、`sd(ddof=1)`，n 不足 4 时如实写 n。")
    L.append("- Δ份额 = `mean_tp − mean_CTRL`。")
    L.append("- 方向计数（2min、8min 各一次）：按 rep 编号配对，`Mem_tp,r − Mem_CTRL,r > 0` 且 `Cyt_tp,r − Cyt_CTRL,r < 0` "
             "的 rep 数，写成 `n/4`（分母 = 两边都有份额的 rep 数，若不足 4 写实际分母并注明）；另报\"仅 Mem 上升\"、\"仅 Cyt 下降\"的 n/4。")
    L.append("- Movement Score（作者定义，§1）：对每个 timepoint：fraction 均值 `m_f,tp = mean over reps of q_f`"
             "（线性值；缺失的 rep 不计入均值，与作者在填补后无缺失的前提不同）；`s_f,tp = m_f,tp / Σ_f m_f,tp`；"
             "`MS_tp = mean_f |s_f,tp − s_f,CTRL|`；同时报 `MS_max1`、`MS_max2`（fraction 编号）。阈值 0.1（:181）。"
             "raw、centered 各算。另加对照列 `MS_perrep_mean`：按 rep 配对算 `mean_f |s_f,tp,r − s_f,CTRL,r|` "
             "再对 rep 取均值（非作者定义，仅对照）。")
    L.append("- t 检验：Mem 份额 tp 组（4 个 rep）vs CTRL 组（4 个 rep），`scipy.stats.ttest_ind(equal_var=False)`（Welch），"
             "双侧；n = 4 / 4，只作参考。")
    L.append("- 两版归一化：`raw` = quantity；`centered` = quantity × 2^a（a 取 `01_run_factors.tsv` 中 layer=proteome、"
             "同 (timepoint, fraction, rep) 的 `a`）。fraction 级数值的 run 解析、基因匹配、值冲突规则与任务 1 相同；"
             "primary group 取 `01_groups.tsv` 中 `is_primary = yes`。\n")
    L.append("### 脚本补充规则（方案未写明或有歧义处；未替换方案定义）\n")
    L.append("- **E-B 分母（歧义）**：方案写\"两边都有份额的 rep 数\"。本脚本取 tp 与 CTRL 两边 `Mem`、`Cyt` 四格都非空的 rep；"
             "三个计数都在这同一组 rep 上数，格式 `n/分母`；分母 = 0 时三格留空。差值 = 0 不计入\"上升\"或\"下降\"。")
    L.append("- **\"仅 Mem 上升\" / \"仅 Cyt 下降\"（歧义）**：按表 E-B 列名 `n_mem_up` / `n_cyt_down` 理解为单看一项"
             "（`Mem_tp,r − Mem_CTRL,r > 0`，不看 Cyt；`Cyt_tp,r − Cyt_CTRL,r < 0`，不看 Mem）。"
             "若\"仅\"指排他（Mem 升且 Cyt 不降），其值 = `n_mem_up − n_mem_up_cyt_down`（同一组 rep），"
             "`n_cyt_down` 同理。")
    L.append("- **E-C 缺整个 fraction（方案未定义）**：tp 或 CTRL 中任一 fraction 在 4 个 rep 全缺，"
             "或 6 个均值之和 = 0 → `MS_author`、`MS_max1`、`MS_max2`、`ge_0.1` 留空，在 ⑤ 列出；不在少于 6 个 fraction 上算。")
    L.append("- **MS_max2**：照作者 `maxN`（:38-45）取第二大的值，再 `match` 第一个等于它的位置；"
             "若最大值并列，`MS_max2` 与 `MS_max1` 相同（作者代码行为）。写成 `FRk`，对应作者的整数 k。")
    L.append("- **MS_perrep_mean**：某 rep 只有在 tp 与 CTRL 两边 6 个 fraction 都有值时才计入；计入 rep 数不足 4 的在 ⑤ 列出。")
    L.append("- **阈值（歧义）**：方案列名 `ge_0.1` 与\"≥0.1 的时间点清单\"用 ≥ 0.1；作者 :181 为严格 `> 0.1`。"
             "本脚本按方案用 ≥；两者只在 MS 恰为 0.1 时不同。作者的 `pval_combi_FDR < 0.05` 条件（limma + sumlog + BH）本轮不复现。")
    L.append("- quantity 按线性值使用（§2\"线性值\"），不做 2^x（作者 :19-20 的 2^ 针对其 log2 输入）。")
    L.append("- E-A：CTRL 行的 `delta_vs_ctrl` = 0；n < 2 时 sd 留空；n = 0 时 mean、delta 留空。")
    L.append("- E-D：tp 组在前、CTRL 组在后（t > 0 表示 tp 的 Mem 均值更高）；任一组 n < 2 或结果为 nan 时 t、p 留空；"
             "只用非空 Mem 份额，n_tp、n_ctrl 为实际个数。")
    L.append("- 02_events.tsv：`key` 在 E-A 中为 `<区室>_<统计量>`（如 `Mem_mean`、`Mem_n`），其余表为列名；"
             "数值 12 位有效数字，空值写空串。md 表中份额、均值、sd、MS 保留 4 位小数，t 保留 3 位，p 保留 3 位有效数字。")
    L.append("- 图：份额取自 `01_shares.tsv`；某 rep 某 tp 份额缺 → 该点断开；配色为 dataviz 参考调色板分类色 1–4（固定顺序），"
             "另以 marker 形状区分 rep。\n")

    def skip():
        if not have_prot:
            L.append(NODATA + "。\n")
        else:
            L.append("缺数据：`01_shares.tsv` 0 行，跳过（请先复跑任务 1）。\n")

    # ③ E-A
    L.append("## ③ 表 E-A：区室份额均值 / 标准差 / Δ份额（primary group，两版）\n")
    if not EA:
        skip()
    else:
        L.append(md_table(EA_COLS, [[g, nm, tp, c, n, f4(f12(m)), f4(f12(sd)), f4(f12(d))]
                                     for g, nm, tp, c, n, m, sd, d in EA]))
        short = [(g, nm, tp, c, n) for g, nm, tp, c, n, *_ in EA if n < len(REPS)]
        L.append(f"\nn < 4 的格：{len(short)}" + ("：" + "；".join(f"{g}/{nm}/{tp}/{c} n={n}" for g, nm, tp, c, n in short)
                                                if short else "") + "。\n")

    # ④ E-B
    L.append("## ④ 表 E-B：方向计数（Mem 上升且 Cyt 下降，按 rep 配对 vs CTRL）\n")
    if not EB:
        skip()
    else:
        L.append(md_table(EB_COLS, EB))
        L.append(f"\n分母不足 4 的行：{len(eb_short)}" + ("：" + "；".join(f"{g}/{nm}/{tp} 分母 {d}" for g, nm, tp, d in eb_short)
                                                  if eb_short else "") + "。\n")

    # ⑤ E-C
    L.append("## ⑤ 表 E-C：Movement Score（作者定义，fraction 级重复均值；对照列 MS_perrep_mean）\n")
    if not EC:
        L.append(NODATA + "。\n")
    else:
        L.append(md_table(EC_COLS, [[g, nm, tp, f4(f12(ms)), mx1 or "", mx2 or "", f4(f12(pm)), ge]
                                     for g, nm, tp, ms, mx1, mx2, pm, ge in EC]))
        hits = [(g, nm, tp) for g, nm, tp, ms, *_ in EC if ms is not None and ms >= MS_THRESHOLD]
        L.append(f"\n**MS_author ≥ 0.1 的 (gene, normalization, timepoint)**：{len(hits)} 个"
                 + ("：" + "；".join(f"({g}, {nm}, {tp})" for g, nm, tp in hits) if hits else "") + "。\n")
        L.append(f"- MS_author 留空（tp 或 CTRL 缺整个 fraction）：{len(ec_missing_fr)}"
                 + ("：" + "；".join(f"{g}/{nm}/{tp}（缺 {w}）" for g, nm, tp, w in ec_missing_fr) if ec_missing_fr else "") + "。")
        L.append(f"- fraction 均值所用 rep 数不足 4（缺失 rep 不计入均值）：{len(ec_partial_mean)}"
                 + ("：" + "；".join(f"{g}/{nm}/{tp}/{fr} n={n}" for g, nm, tp, fr, n in ec_partial_mean)
                    if ec_partial_mean else "") + "。")
        L.append(f"- MS_perrep_mean 计入 rep 数不足 4：{len(ec_perrep_short)}"
                 + ("：" + "；".join(f"{g}/{nm}/{tp} n={n}" for g, nm, tp, n in ec_perrep_short)
                    if ec_perrep_short else "") + "。\n")

    # ⑥ E-D
    L.append("## ⑥ 表 E-D：Mem 份额 tp vs CTRL，Welch t 检验（n = 4 / 4，只作参考）\n")
    if not ED:
        skip()
    else:
        L.append(md_table(ED_COLS, [[g, nm, tp, f4(f12(a)), f4(f12(b)),
                                     "" if t is None else f"{t:.3f}", fsig(f12(p)), nt, nc]
                                    for g, nm, tp, a, b, t, p, nt, nc in ED]))
        L.append("\n`scipy.stats.ttest_ind(equal_var=False)`，双侧，未做多重检验校正；每组 n ≤ 4，p 值只作参考。\n")

    # ⑦ EGFR
    L.append("## ⑦ EGFR 总量（抄自 `03_egfr_total.md` 的 ③④ 两节）\n")
    if egfr_block is None:
        L.append(f"缺 `{rel(p_egfr)}` 的 ③④ 节，跳过。\n")
    else:
        L.append(f"以下原样抄自 `{rel(p_egfr)}`（只把标题降一级，正文未改）。\n")
        L.extend(egfr_block)
        L.append("")

    # ⑧ 跳过项
    L.append("## ⑧ 跳过项\n")
    if not have_prot:
        L.append(f"- {NODATA}：③ E-A、④ E-B、⑤ E-C、⑥ E-D 均无数值；`02_events.tsv` 只有表头；"
                 "`02_events.png` 只有 2 × 9 图框，子图内写 \"no data\"。数据到位、任务 1 复跑后直接复跑本脚本即可。")
        L.append(f"- 四个蛋白各在 run 测到：" + "，".join(f"{g} {len(gene_detected[g])}/{N_DESIGN}" for g in GENES_ALL) + "。")
    if have_prot and not have_shares:
        L.append("- `01_shares.tsv` 0 行而 `gate1_proteins.tsv` 非空：E-A、E-B、E-D 与图跳过，请先复跑任务 1。")
    L.append("- 不复现 DAPAR 的 LOESS 归一化与 KNN/MEC/detQuant 填补（`00_design.md` §1），只用 raw 与 centered 两版；不填补缺失。")
    L.append("- 不复现作者 :58-153 的 limma p 值合并与 :173 的 BH 校正（无 `pval_combi_FDR`）；E-C 只按 MS 阈值列出。")
    L.append("- 脚本在 scratchpad 合成夹具上跑通。\n")
    L.append("### 自检（方案验收项）\n")
    L.append(md_table(["检查", "结果", "值"], checks))
    L.append("")
    o_md.write_text("\n".join(L) + "\n", encoding="utf-8")

    # ================================================================ stdout
    print(f"input rows {n_input}; 01_shares rows {n_shares}; runs parsed {len(run_trip)}/{N_DESIGN}")
    print("detected runs: " + ", ".join(f"{g} {len(gene_detected[g])}/{N_DESIGN}" for g in GENES_ALL))
    print("rows: " + ", ".join(f"{k} {v}" for k, v in got_rows.items()) + f"; 02_events.tsv {len(R_L)}")
    for name, st, val in checks:
        print(f"  [{st}] {name} {val}")
    print("wrote:", ", ".join(rel(p) for p in (o_tsv, o_png, o_md)))
    return 0 if all(st != "FAIL" for _, st, _ in checks) else 1


if __name__ == "__main__":
    sys.exit(main())
