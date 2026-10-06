#!/usr/bin/env python3
"""09_gate1b 任务 2（Q2c）：并列表 + 独立复核（Python）。

先运行 02_author_branch.R，再运行本脚本：
    python3 /home/user/phd-script/analysis/09_gate1b/02_author_compare.py

读：analysis/09_gate1b/02_author_branch.tsv、02_author_cells.tsv、02_author_merged_log2.tsv、
    02_author_limma_FR<j>_<tp>vsCTRL.tsv（24 个）；analysis/08_gate1/02_events.tsv（表 E-C）；
    data/pilot/proteome_pg_matrix_120.tsv（只用于核对四蛋白原始非空格数）。
写：analysis/09_gate1b/02_author_compare.tsv；并把 02_author_branch.md 的 ⑥ 节
    （<!-- SECTION6_BEGIN --> 与 <!-- SECTION6_END --> 之间）替换为并列表与 Python 自检。
只用全局 pandas / numpy（不装包）；读表 sep='\t', dtype=str, keep_default_na=False，空串 = 缺失。
"""
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path("/home/user/phd-script")
OUT = REPO / "analysis" / "09_gate1b"
F_BRANCH = OUT / "02_author_branch.tsv"
F_CELLS = OUT / "02_author_cells.tsv"
F_MERGED = OUT / "02_author_merged_log2.tsv"
F_MD = OUT / "02_author_branch.md"
F_COMPARE = OUT / "02_author_compare.tsv"
F_EVENTS = REPO / "analysis" / "08_gate1" / "02_events.tsv"
F_PG = REPO / "data" / "pilot" / "proteome_pg_matrix_120.tsv"

GENES4 = ["GRB2", "SHC1", "CBL", "EGFR"]
TPS = ["2min", "8min", "20min", "90min"]
FRS = [f"FR{i}" for i in range(1, 7)]
REPS = [f"Rep{i}" for i in range(1, 5)]
COLS = ["gene", "timepoint", "MS_author_pipeline", "MS_gate1_raw", "MS_gate1_centered",
        "MS_perrep_mean_gate1_raw", "n_obs_cells_tp", "n_obs_cells_ctrl",
        "n_imputed_cells_tp", "n_imputed_cells_ctrl"]
NO_GATE1 = "无"


def read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def parse_design(run):  # 07 :32-38 原样
    tp = re.search(r"_(2min|8min|20min|90min|CTRL)_", run, re.I)
    fr = re.search(r"_(FR\d)_", run)
    rep = re.search(r"_(Rep\d)", run)
    if not (tp and fr and rep):
        return None
    return (tp.group(1), fr.group(1), rep.group(1))


def bh(p):
    """R p.adjust(method = "BH")。"""
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(-p, kind="mergesort")
    i = np.arange(n, 0, -1)
    q = np.minimum.accumulate(n / i * p[o])
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def sumlog_p(p1, p2):
    """metap::sumlog(c(p1, p2))$p：chisq = -2 Σ ln p，df = 4 → 生存函数闭式 exp(-x/2)(1 + x/2)。"""
    x = -2.0 * (math.log(p1) + math.log(p2))
    return math.exp(-x / 2.0) * (1.0 + x / 2.0)


checks = []


def check(name, ok, value=""):
    checks.append((name, "PASS" if ok else "FAIL", str(value)))


# ------------------------------------------------------------------ 读表
br = read(F_BRANCH)
cells = read(F_CELLS)
ev = read(F_EVENTS)
merged = read(F_MERGED)
assert len(br) == 16 and len(cells) == 480

ec = ev[ev["table"] == "E-C"]
ec_key = {(r.gene, r.normalization, r.timepoint, r.key): r.value for r in ec.itertuples()}
ec_genes = sorted(set(ec["gene"]))

# ------------------------------------------------------------------ 并列表
rows = []
for g in GENES4:
    for tp in TPS:
        b = br[(br["gene"] == g) & (br["timepoint"] == tp)]
        assert len(b) == 1
        c_tp = cells[(cells["gene"] == g) & (cells["timepoint"] == tp)]
        c_ct = cells[(cells["gene"] == g) & (cells["timepoint"] == "CTRL")]
        if g in ec_genes:
            ms_raw = ec_key[(g, "raw", tp, "MS_author")]
            ms_cen = ec_key[(g, "centered", tp, "MS_author")]
            ms_prm = ec_key[(g, "raw", tp, "MS_perrep_mean")]
        else:
            ms_raw = ms_cen = ms_prm = NO_GATE1
        rows.append({
            "gene": g, "timepoint": tp,
            "MS_author_pipeline": b["MS_author_pipeline"].iloc[0],
            "MS_gate1_raw": ms_raw, "MS_gate1_centered": ms_cen, "MS_perrep_mean_gate1_raw": ms_prm,
            "n_obs_cells_tp": int((c_tp["source"] == "observed").sum()),
            "n_obs_cells_ctrl": int((c_ct["source"] == "observed").sum()),
            "n_imputed_cells_tp": int(c_tp["source"].isin(["KNN", "detQuant"]).sum()),
            "n_imputed_cells_ctrl": int(c_ct["source"].isin(["KNN", "detQuant"]).sum()),
            "_fo_tp": int((c_tp["source"] == "filtered_out").sum()),
            "_fo_ctrl": int((c_ct["source"] == "filtered_out").sum()),
        })
cmp_df = pd.DataFrame(rows)
cmp_df[COLS].to_csv(F_COMPARE, sep="\t", index=False)

# ------------------------------------------------------------------ 自检
cmp_back = read(F_COMPARE)
check("02_author_compare.tsv：16 行，列 = 方案列，(gene, tp) 齐全",
      len(cmp_back) == 16 and list(cmp_back.columns) == COLS
      and set(zip(cmp_back["gene"], cmp_back["timepoint"])) == {(g, t) for g in GENES4 for t in TPS},
      f"{len(cmp_back)} 行")
mis = 0
for r in cmp_back.itertuples():
    if r.gene in ec_genes:
        for col, (norm, key) in {"MS_gate1_raw": ("raw", "MS_author"),
                                 "MS_gate1_centered": ("centered", "MS_author"),
                                 "MS_perrep_mean_gate1_raw": ("raw", "MS_perrep_mean")}.items():
            if getattr(r, col) != ec_key[(r.gene, norm, r.timepoint, key)]:
                mis += 1
    else:
        mis += sum(getattr(r, c) != NO_GATE1 for c in ["MS_gate1_raw", "MS_gate1_centered", "MS_perrep_mean_gate1_raw"])
check("gate1 三列与 08_gate1/02_events.tsv（E-C）逐格一致（EGFR 不在 E-C → 无）",
      mis == 0, f"E-C 基因 {','.join(ec_genes)}；不一致 {mis}")
mis = sum(cmp_back["MS_author_pipeline"].values != br.set_index(["gene", "timepoint"]).loc[
    list(zip(cmp_back["gene"], cmp_back["timepoint"])), "MS_author_pipeline"].values)
check("MS_author_pipeline 与 02_author_branch.tsv 一致", mis == 0, f"不一致 {mis}")
ok24 = all(r["n_obs_cells_tp"] + r["n_imputed_cells_tp"] + r["_fo_tp"] == 24
           and r["n_obs_cells_ctrl"] + r["n_imputed_cells_ctrl"] + r["_fo_ctrl"] == 24 for r in rows)
check("n_obs + n_imputed + filtered_out = 24（每 tp / CTRL）", ok24, "")

# cells ↔ merged（独立回读）
run_map = {}
for c in merged.columns[1:]:
    d = parse_design(c)
    assert d is not None
    run_map[d] = c
assert len(run_map) == 120
mrow = {g: i for i, g in enumerate(merged["gene"])}
mis = 0
for r in cells.itertuples():
    col = run_map[(r.timepoint, r.fraction, r.rep)]
    mv = merged.at[mrow[r.gene], col] if r.gene in mrow else ""
    if mv != r.value_log2_final:
        mis += 1
    if (r.source == "filtered_out") != (r.value_log2_final == ""):
        mis += 1
check("02_author_cells 值 = 02_author_merged_log2 对应格（Python 独立回读；filtered_out ⇔ 空）",
      mis == 0, f"480 格，不一致 {mis}")

# 四蛋白 observed 格数 = pg_matrix 原始非空 − 没过 S3 fraction 中的非空
pg = read(F_PG)
runs = [c for c in pg.columns if c not in ("Protein.Group", "Genes")]
des = {c: parse_design(c) for c in runs}
txt = []
ok_obs = True
for g in GENES4:
    sub = pg[pg["Genes"] == g]
    assert len(sub) == 1
    nonempty = {des[c] for c in runs if sub[c].iloc[0] != ""}
    fo = cells[(cells["gene"] == g) & (cells["source"] == "filtered_out")]
    fo_keys = set(zip(fo["timepoint"], fo["fraction"], fo["rep"]))
    obs = cells[(cells["gene"] == g) & (cells["source"] == "observed")]
    obs_keys = set(zip(obs["timepoint"], obs["fraction"], obs["rep"]))
    ok_obs &= obs_keys == (nonempty - fo_keys)
    txt.append(f"{g} 原始非空 {len(nonempty)}，observed {len(obs_keys)}，filtered_out 中原始非空 {len(nonempty & fo_keys)}")
check("四蛋白 observed 格 = pg_matrix 原始非空格 − filtered_out 格（逐格集合相等）", ok_obs, "；".join(txt))

# T1–T10 独立复算（合并矩阵 → MS / max1 / max2 / p / sumlog / BH）
num = merged.set_index("gene").replace("", np.nan).astype(float)
complete = num.dropna(how="any")
lim = {}
for j, fr in enumerate(FRS, start=1):
    for tp in TPS:
        t = read(OUT / f"02_author_limma_FR{j}_{tp}vsCTRL.tsv")
        lim[(fr, tp)] = dict(zip(t["gene"], t["P_Value"].astype(float)))
max_ms = max_sh = max_p = max_fdr = 0.0
idx_mis = 0
n_fdr = set()
br_i = br.set_index(["gene", "timepoint"])
for tp in TPS:
    def frmean(cond):
        cols = [[run_map[(cond, fr, rp)] for rp in REPS] for fr in FRS]
        return np.column_stack([np.power(2.0, complete[c].values).mean(axis=1) for c in cols])
    A, B = frmean(tp), frmean("CTRL")
    As = A / A.sum(axis=1, keepdims=True)
    Bs = B / B.sum(axis=1, keepdims=True)
    MS = np.abs(As - Bs)
    msm = MS.mean(axis=1)
    m1 = np.array([int(np.flatnonzero(x == x.max())[0]) for x in MS])
    m2 = np.array([int(np.flatnonzero(x == np.sort(x)[-2])[0]) for x in MS])
    genes = list(complete.index)
    p1 = np.array([lim[(FRS[k], tp)].get(g, 1.0) for k, g in zip(m1, genes)])
    p2 = np.array([lim[(FRS[k], tp)].get(g, 1.0) for k, g in zip(m2, genes)])
    pc = np.array([sumlog_p(a, b) for a, b in zip(p1, p2)])
    fdr = bh(pc)
    n_fdr.add(len(pc))
    for g in GENES4:
        b = br_i.loc[(g, tp)]
        if g in complete.index:
            i = genes.index(g)
            max_ms = max(max_ms, abs(msm[i] - float(b["MS_author_pipeline"])))
            max_sh = max(max_sh, abs(As[i, m1[i]] - float(b["share_max1_tp"])),
                         abs(Bs[i, m1[i]] - float(b["share_max1_ctrl"])))
            idx_mis += (FRS[m1[i]] != b["MS_max1"]) + (FRS[m2[i]] != b["MS_max2"])
            max_p = max(max_p, abs(p1[i] - float(b["p_max1"])), abs(p2[i] - float(b["p_max2"])),
                        abs(pc[i] - float(b["pval_combi"])) / max(pc[i], 1e-300))
            max_fdr = max(max_fdr, abs(fdr[i] - float(b["pval_combi_FDR"])))
        else:
            idx_mis += int(b["MS_author_pipeline"] != "")
check("T2–T4 用合并矩阵独立复算 MS：16 行一致（不在完整行内的蛋白 tsv 为空）",
      max_ms < 1e-12 and idx_mis == 0, f"max|ΔMS| = {max_ms:.1e}，max|Δ份额| = {max_sh:.1e}，max1/max2 不一致 {idx_mis}")
check("T6–T9 独立复算 p_max1 / p_max2 / pval_combi / pval_combi_FDR 一致",
      max_p < 1e-9 and max_fdr < 1e-9 and n_fdr == {int(x) for x in br["n_proteins_in_FDR"]},
      f"max|Δp|（pval_combi 为相对差）= {max_p:.1e}，max|ΔFDR| = {max_fdr:.1e}，BH 的 n = {sorted(n_fdr)}")

import matplotlib  # noqa: E402  只为记录全局版本
vers = f"pandas {pd.__version__} / numpy {np.__version__} / matplotlib {matplotlib.__version__}"
check("全局 Python 版本 = pandas 2.2.3 / numpy 1.26.4 / matplotlib 3.8.4（本脚本运行时）",
      (pd.__version__, np.__version__, matplotlib.__version__) == ("2.2.3", "1.26.4", "3.8.4"), vers)

# ------------------------------------------------------------------ md ⑥
def md_table(df):
    h = "| " + " | ".join(df.columns) + " |"
    s = "|" + "|".join(["---"] * len(df.columns)) + "|"
    body = ["| " + " | ".join(str(v).replace("|", "\\|") for v in r) + " |" for r in df.itertuples(index=False)]
    return [h, s] + body


show = cmp_df.copy()
for c in ["MS_author_pipeline", "MS_gate1_raw", "MS_gate1_centered", "MS_perrep_mean_gate1_raw"]:
    show[c] = [("" if v == "" else v if v == NO_GATE1 else f"{float(v):.4f}") for v in show[c]]
show["filtered_out_tp/ctrl"] = [f"{a}/{b}" for a, b in zip(show["_fo_tp"], show["_fo_ctrl"])]
show = show[COLS + ["filtered_out_tp/ctrl"]]

sec = [
    "MS_author_pipeline 来自本轮 `02_author_branch.tsv`；MS_gate1_raw / MS_gate1_centered / MS_perrep_mean_gate1_raw 原样抄 "
    "`analysis/08_gate1/02_events.tsv` 表 E-C 的 `MS_author`（raw / centered）与 `MS_perrep_mean`（raw）；"
    "空 = 原表为空串（gate1 未算；MS_author 为空的原因是该 tp 或 CTRL 有整个 fraction 4 个 rep 全缺，MS_perrep_mean 为空的原因是没有 rep 在两边 6 个 fraction 都有值，见 08_gate1/02_events.md ②）；`无` = gate1 E-C 无该基因（EGFR）。"
    "n_obs / n_imputed 为该蛋白在该 tp（或 CTRL）24 格（6 FR × 4 rep）中来源为 observed / KNN + detQuant 的格数；"
    "filtered_out 格两者都不计（最后一列为 filtered_out 格数，只在 md 中列出，tsv 无此列）。MS 4 位小数，全精度见 tsv。",
    "",
]
sec += md_table(show)
sec += ["", "### 自检（Python 部分）", ""]
sec += md_table(pd.DataFrame(checks, columns=["check", "result", "value"]))
sec += ["", f"- 全局 Python（本脚本运行时）：{vers}。", ""]

md = F_MD.read_text(encoding="utf-8").split("\n")
b = md.index("<!-- SECTION6_BEGIN -->")
e = md.index("<!-- SECTION6_END -->")
md = md[:b + 1] + sec + md[e:]
F_MD.write_text("\n".join(md), encoding="utf-8")

for c in checks:
    print("\t".join(c))
print(cmp_df[COLS].to_string(index=False))
sys.exit(0 if all(c[1] == "PASS" for c in checks) else 1)
