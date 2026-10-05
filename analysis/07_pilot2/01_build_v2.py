#!/usr/bin/env python3
"""07_pilot2 任务 1：统一蛋白分母 -> 01_site_traj_v2.tsv 等。

按 analysis/07_pilot2/00_design.md 的 "0. 输入现状" 与 "任务 1" 执行。
规则来源：code/Protein contour/Zhihan/Test/07_main_analysis.py :238-306
          （median :117-118、parse_design :32-38 原样照抄）；
分母来源：code/Protein contour/Zhihan/Test/files1/06_occupancy_variants.py :185-227。

输入（只读）：data/pilot/occupancy_variants_long.tsv, data/pilot/run_medians.tsv,
              data/pilot/site_traj.tsv
输出：analysis/07_pilot2/01_run_factors.tsv, 01_site_traj_v2.tsv,
      01_class_compare.tsv, 01_class_compare.md
复跑：python3 analysis/07_pilot2/01_build_v2.py   （任意目录均可；确定性）
"""
import math
import os
import re
import statistics
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data", "pilot")
P_LONG = os.path.join(DATA, "occupancy_variants_long.tsv")
P_RM = os.path.join(DATA, "run_medians.tsv")
P_TRAJ = os.path.join(DATA, "site_traj.tsv")
O_FACT = os.path.join(HERE, "01_run_factors.tsv")
O_V2 = os.path.join(HERE, "01_site_traj_v2.tsv")
O_CMP = os.path.join(HERE, "01_class_compare.tsv")
O_MD = os.path.join(HERE, "01_class_compare.md")

# 07 :22-25 已有常量，不改
FC_THRESH = 1.0
MIN_REPS = 2
TPS = ["CTRL", "2min", "8min", "20min", "90min"]
EGF_TPS = ["2min", "8min", "20min", "90min"]
FRS = [f"FR{i}" for i in range(1, 7)]
REPS = [f"Rep{i}" for i in range(1, 5)]
CLASSES = ["both", "occupancy_only", "intensity_only", "none"]
RESP = {"both", "occupancy_only"}
OCC_COLS = [f"log2occ_{t}" for t in TPS]
INT_COLS = [f"log2int_{t}" for t in TPS]
COLS = ["site_id", "protein", "fraction", "denom_tier", "class", "peak_tp"] + OCC_COLS + INT_COLS


def median(v):  # 07 :117-118 原样
    return sorted(v)[len(v) // 2] if len(v) % 2 else sum(sorted(v)[len(v) // 2 - 1:len(v) // 2 + 1]) / 2


def parse_design(run):  # 07 :32-38 原样
    tp = re.search(r"_(2min|8min|20min|90min|CTRL)_", run, re.I)
    fr = re.search(r"_(FR\d)_", run)
    rep = re.search(r"_(Rep\d)", run)
    if not (tp and fr and rep):
        return None
    return (tp.group(1), fr.group(1), rep.group(1))


def fmt3(v):
    # 07 :294-295 round(median, 3) 经 csv.writer 写出 = str(float)；None -> ''
    return "" if v is None else str(round(v, 3))


def read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


# ------------------------------------------------------------------ 读表
L = read(P_LONG)
T = read(P_TRAJ)
RM = read(P_RM)
assert len(L) == 306035 and len(T) == 39813 and len(RM) == 240
assert list(T.columns) == COLS
assert not L.duplicated(["site_id", "timepoint", "fraction", "rep"]).any()
assert (L["phos_intensity"] != "").all()
T_KEYS = list(zip(T["site_id"], T["fraction"]))
assert len(set(T_KEYS)) == 39813
assert set(zip(L["site_id"], L["fraction"])) == set(T_KEYS)
TV = {c: T[c].tolist() for c in T.columns}
T_IDX = {k: i for i, k in enumerate(T_KEYS)}

# ------------------------------------------------------------------ 规则 1：a0
RM["layer"] = RM["run"].str.split("/").str[0]
RM["key"] = RM["run"].map(parse_design)
assert RM["key"].notna().all()
RM["m"] = RM["median_log2"].astype(float)
anchor, a0 = {}, {"phospho": {}, "proteome": {}}
run_name, run_m = {}, {}
for lay in ("phospho", "proteome"):
    g = RM[RM["layer"] == lay]
    assert len(g) == 120 and g["key"].nunique() == 120
    anchor[lay] = statistics.median(g["m"].tolist())
    for run, k, m in zip(g["run"], g["key"], g["m"]):
        a0[lay][k] = anchor[lay] - m
        run_name[(lay, k)] = run
        run_m[(lay, k)] = m

# ------------------------------------------------------------------ 规则 2：轨迹/观测
# obs = (tp, rep, log2(phos), log2(denom_residue)|None, log2(denom_protein)|None)
traj = defaultdict(list)
for sid, tp, fr, rep, p, dr, dp in zip(L["site_id"], L["timepoint"], L["fraction"], L["rep"],
                                         L["phos_intensity"], L["denom_residue"], L["denom_protein"]):
    traj[(sid, fr)].append((tp, rep, math.log2(float(p)),
                            math.log2(float(dr)) if dr else None,
                            math.log2(float(dp)) if dp else None))


def old_tier(obs):  # 07 :262-263
    n_res = sum(1 for o in obs if o[3] is not None)
    return "residue" if n_res >= 0.5 * len(obs) else "protein"


def responder(vals):  # 07 :276-285
    if len(vals.get("CTRL", [])) < MIN_REPS:
        return False, ""
    base = median(vals["CTRL"])
    for tp in EGF_TPS:
        v = vals.get(tp, [])
        if len(v) < MIN_REPS:
            continue
        fc = median(v) - base
        if abs(fc) >= FC_THRESH and all((x - base) * fc > 0 for x in v):
            return True, tp
    return False, ""


def classify(fr, obs, ap, aq, mode):
    """mode 'old' = 07 规则 3；mode 'new' = 所有观测 d = den_prot。"""
    tier = old_tier(obs) if mode == "old" else "protein"
    occ, inten = defaultdict(list), defaultdict(list)
    n_occ_obs = 0
    for tp, rep, lp, ldr, ldp in obs:
        k = (tp, fr, rep)
        x = lp + ap[k]                      # log2(num)
        inten[tp].append(x)
        ld = ldr if tier == "residue" else ldp
        if ld is not None:                  # 07 :271 `if d:`
            occ[tp].append(x - (ld + aq[k]))  # log2(num / d)
            n_occ_obs += 1
    o_resp, o_tp = responder(occ)
    i_resp, i_tp = responder(inten)
    cls = ("both" if o_resp and i_resp else
           "occupancy_only" if o_resp else
           "intensity_only" if i_resp else "none")
    om = [median(occ[t]) if len(occ.get(t, [])) >= MIN_REPS else None for t in TPS]
    im = [median(inten[t]) if len(inten.get(t, [])) >= MIN_REPS else None for t in TPS]
    return dict(tier=tier, cls=cls, peak=o_tp or i_tp, om=om, im=im,
                n_obs=len(obs), n_occ_obs=n_occ_obs, occ=occ)


def run_all(ap, aq, mode):
    return [classify(fr, traj[(sid, fr)], ap, aq, mode) for sid, fr in T_KEYS]


def thou(s):
    return int(round(float(s) * 1000))


def compare_old(rows):
    """与 site_traj 比：class/denom_tier/peak_tp 一致行数；log2 列有无与差值（单位 0.001）。"""
    r = Counter()
    mism = []
    for i, x in enumerate(rows):
        r["class"] += x["cls"] == TV["class"][i]
        r["tier"] += x["tier"] == TV["denom_tier"][i]
        r["peak"] += x["peak"] == TV["peak_tp"][i]
        pres_ok = True
        for c, v in zip(OCC_COLS + INT_COLS, x["om"] + x["im"]):
            s = TV[c][i]
            if (s == "") != (v is None):
                pres_ok = False
                r["cells_presence_diff"] += 1
                continue
            if v is None:
                continue
            d = abs(thou(fmt3(v)) - thou(s))
            r["cells_nonempty"] += 1
            r["maxdiff_thou"] = max(r["maxdiff_thou"], d)
            if d == 1:
                r["cells_diff_eq_0.001"] += 1
            elif d > 1:
                r["cells_diff_gt_0.001"] += 1
        r["rows_presence_ok"] += pres_ok
        if (x["cls"] != TV["class"][i] or x["tier"] != TV["denom_tier"][i]
                or x["peak"] != TV["peak_tp"][i]):
            mism.append(i)
    return r, mism


# ------------------------------------------------------------------ 步骤 0a：a0 复现
ap0, aq0 = a0["phospho"], a0["proteome"]
rows_a0 = run_all(ap0, aq0, "old")
cmp_a0, mism_a0 = compare_old(rows_a0)


# ------------------------------------------------------------------ 步骤 0b：解 δ
def solve(eqs, sign):
    """eqs[cell] = [(rep_i, rep_j, diff)]；模型 diff = sign * (δ_i + δ_j) / 2；逐 cell 最小二乘。"""
    out = {}
    for cell in sorted(eqs):
        lst = sorted(eqs[cell])
        A = np.zeros((len(lst), 4))
        b = np.zeros(len(lst))
        for n, (ri, rj, df) in enumerate(lst):
            A[n, REPS.index(ri)] += sign * 0.5
            A[n, REPS.index(rj)] += sign * 0.5
            b[n] = df
        rank = int(np.linalg.matrix_rank(A))
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        res = b - A @ x
        per_run = Counter()
        for ri, rj, _ in lst:
            per_run[ri] += 1
            per_run[rj] += 1
        out[cell] = dict(delta=dict(zip(REPS, x.tolist())), n=len(lst), rank=rank,
                         maxres=float(np.abs(res).max()), per_run=per_run,
                         pairs=Counter((ri, rj) for ri, rj, _ in lst))
    return out


def groups_by_tp(obs):
    by = defaultdict(list)
    for o in obs:
        by[o[0]].append(o)
    return by


# phospho：log2int_tp 非空 且 长表恰 2 个重复
eq_p = defaultdict(list)
for (sid, fr), obs in traj.items():
    i = T_IDX[(sid, fr)]
    for tp, g in groups_by_tp(obs).items():
        if len(g) != 2 or TV[f"log2int_{tp}"][i] == "":
            continue
        (_, r1, l1, _, _), (_, r2, l2, _, _) = sorted(g, key=lambda o: o[1])
        diff = float(TV[f"log2int_{tp}"][i]) - median([l1 + ap0[(tp, fr, r1)], l2 + ap0[(tp, fr, r2)]])
        eq_p[(fr, tp)].append((r1, r2, diff))
sol_p = solve(eq_p, +1.0)
ap = {k: ap0[k] + sol_p[(k[1], k[0])]["delta"][k[2]] if sol_p[(k[1], k[0])]["rank"] == 4 else ap0[k]
      for k in ap0}
rows_p_only = run_all(ap, aq0, "old")
cmp_p_only, _ = compare_old(rows_p_only)

# proteome：旧档位定分母，log2occ_tp 非空 且 恰 2 个重复有分母
eq_q = defaultdict(list)
for (sid, fr), obs in traj.items():
    i = T_IDX[(sid, fr)]
    tier = old_tier(obs)
    by = defaultdict(list)
    for tp, rep, lp, ldr, ldp in obs:
        ld = ldr if tier == "residue" else ldp
        if ld is not None:
            by[tp].append((rep, lp, ld))
    for tp, g in by.items():
        if len(g) != 2 or TV[f"log2occ_{tp}"][i] == "":
            continue
        (r1, l1, d1), (r2, l2, d2) = sorted(g)
        v1 = l1 + ap[(tp, fr, r1)] - d1 - aq0[(tp, fr, r1)]
        v2 = l2 + ap[(tp, fr, r2)] - d2 - aq0[(tp, fr, r2)]
        diff = float(TV[f"log2occ_{tp}"][i]) - median([v1, v2])
        eq_q[(fr, tp)].append((r1, r2, diff))
sol_q = solve(eq_q, -1.0)
aq = {k: aq0[k] + sol_q[(k[1], k[0])]["delta"][k[2]] if sol_q[(k[1], k[0])]["rank"] == 4 else aq0[k]
      for k in aq0}

N_CELLS_LOW_RANK = sum(1 for s in list(sol_p.values()) + list(sol_q.values()) if s["rank"] < 4)
assert len(sol_p) == 30 and len(sol_q) == 30

rows_cal = run_all(ap, aq, "old")
cmp_cal, mism_cal = compare_old(rows_cal)

# 01_run_factors.tsv
fact_rows = []
for lay, sol, a_fin in (("phospho", sol_p, ap), ("proteome", sol_q, aq)):
    for fr in FRS:
        for tp in TPS:
            s = sol[(fr, tp)]
            for rep in REPS:
                k = (tp, fr, rep)
                fact_rows.append([lay, run_name[(lay, k)], tp, fr, rep, f"{run_m[(lay, k)]:.4f}",
                                  f"{a0[lay][k]:.6f}", f"{a_fin[k] - a0[lay][k]:.6f}", f"{a_fin[k]:.6f}",
                                  s["n"], f"{s['maxres']:.6f}"])
FACT = pd.DataFrame(fact_rows, columns=["layer", "run", "timepoint", "fraction", "rep", "median_log2",
                                         "a0", "delta", "a", "n_equations", "max_residual"])
FACT.to_csv(O_FACT, sep="\t", index=False, lineterminator="\n")

# ------------------------------------------------------------------ 步骤 1：蛋白档
rows_new = run_all(ap, aq, "new")
v2 = []
for (sid, fr), x in zip(T_KEYS, rows_new):
    v2.append([sid, sid.rsplit("_", 1)[0], fr, x["tier"], x["cls"], x["peak"]]
              + [fmt3(v) for v in x["om"]] + [fmt3(v) for v in x["im"]])
V2 = pd.DataFrame(v2, columns=COLS)
assert (V2["protein"] == T["protein"]).all()
V2.to_csv(O_V2, sep="\t", index=False, lineterminator="\n")

# 对照：新规则若用未校准 a0（仅用于说明校准对新分类的影响）
rows_new_a0 = run_all(ap0, aq0, "new")
new_a0_cls = [x["cls"] for x in rows_new_a0]

# ------------------------------------------------------------------ 步骤 2：比较
old_cls = TV["class"]
rep_cls = [x["cls"] for x in rows_cal]
new_cls = V2["class"].tolist()
C1 = {c: (old_cls.count(c), rep_cls.count(c), new_cls.count(c)) for c in CLASSES}
C2 = Counter(zip(old_cls, new_cls))
n_changed = sum(1 for a, b in zip(old_cls, new_cls) if a != b)
site_old = {s for s, c in zip(TV["site_id"], old_cls) if c in RESP}
site_new = {s for s, c in zip(TV["site_id"], new_cls) if c in RESP}
C3 = dict(n_old=len(site_old), n_new=len(site_new), both=len(site_old & site_new),
          old_only=len(site_old - site_new), new_only=len(site_new - site_old))
C4 = {}
for fr in FRS + ["all"]:
    sel = [i for i, f in enumerate(TV["fraction"]) if fr == "all" or f == fr]
    C4[fr] = (sum(old_cls[i] in RESP for i in sel), sum(new_cls[i] in RESP for i in sel))

occ_ctrl_old = sum(1 for s in TV["log2occ_CTRL"] if s != "")
occ_ctrl_rep = sum(1 for x in rows_cal if x["om"][0] is not None)
occ_ctrl_new = sum(1 for s in V2["log2occ_CTRL"] if s != "")
n_obs_total = sum(x["n_obs"] for x in rows_new)
n_obs_no_dprot = int((L["denom_protein"] == "").sum())
n_obs_occ_new = sum(x["n_occ_obs"] for x in rows_new)
n_obs_occ_old = sum(x["n_occ_obs"] for x in rows_cal)

# log2int 新(v2) vs 旧(site_traj) 逐格；以及 新 vs 旧复现
int_cells = int_str_diff = int_pres_diff = int_num_diff = int_maxd = 0
int_vs_rep_diff = 0
for i, x in enumerate(rows_new):
    for j, c in enumerate(INT_COLS):
        s_new, s_old, s_rep = V2[c].iat[i], TV[c][i], fmt3(rows_cal[i]["im"][j])
        int_vs_rep_diff += s_new != s_rep
        if s_new == "" and s_old == "":
            continue
        int_cells += 1
        int_str_diff += s_new != s_old
        if (s_new == "") != (s_old == ""):
            int_pres_diff += 1
            continue
        d = abs(thou(s_new) - thou(s_old))
        int_num_diff += d > 0
        int_maxd = max(int_maxd, d)

# 自检
checks = []
checks.append(("v2 行数 39813、16 列", len(V2) == 39813 and V2.shape[1] == 16))
checks.append(("v2 (site_id, fraction) 集合与顺序同 site_traj",
               list(zip(V2["site_id"], V2["fraction"])) == T_KEYS))
checks.append(("v2 denom_tier 全为 protein", (V2["denom_tier"] == "protein").all()))
checks.append(("log2int_* 与 site_traj 逐格相同（字符串）", int_str_diff == 0))
checks.append(("log2int_* 与 步骤0 校准复现 逐格相同（字符串）", int_vs_rep_diff == 0))
checks.append(("校准后 class 39813/39813", cmp_cal["class"] == 39813))
checks.append(("校准后 denom_tier 39813/39813", cmp_cal["tier"] == 39813))
checks.append(("校准后 peak_tp 39813/39813", cmp_cal["peak"] == 39813))
checks.append(("校准后 log2 有无一致（单元格）", cmp_cal["cells_presence_diff"] == 0))
checks.append(("校准后 非空 log2 差 ≤ 0.001", cmp_cal["maxdiff_thou"] <= 1))
checks.append(("所有 cell 残差最大值 ≤ 0.0015",
               max(s["maxres"] for s in list(sol_p.values()) + list(sol_q.values())) <= 0.0015))
checks.append(("C1 三列合计各 39813", all(sum(C1[c][j] for c in CLASSES) == 39813 for j in range(3))))
checks.append(("C2 行边际 = n_old、列边际 = n_new",
               all(sum(C2[(a, b)] for b in CLASSES) == C1[a][0] for a in CLASSES)
               and all(sum(C2[(a, b)] for a in CLASSES) == C1[b][2] for b in CLASSES)))
checks.append(("C3 两者皆+仅旧 = n_old、两者皆+仅新 = n_new",
               C3["both"] + C3["old_only"] == C3["n_old"] and C3["both"] + C3["new_only"] == C3["n_new"]))
checks.append(("C4 六个 fraction 之和 = all",
               sum(C4[f][0] for f in FRS) == C4["all"][0] and sum(C4[f][1] for f in FRS) == C4["all"][1]))
checks.append(("新 log2occ_CTRL 非空行数 ≥ 旧", occ_ctrl_new >= occ_ctrl_old))

# ------------------------------------------------------------------ 01_class_compare.tsv（长格式）
cmp_rows = []
for c in CLASSES:
    for j, col in enumerate(["n_old_site_traj", "n_old_replicated_step0", "n_new"]):
        cmp_rows.append(["C1", c, col, C1[c][j]])
for j, col in enumerate(["n_old_site_traj", "n_old_replicated_step0", "n_new"]):
    cmp_rows.append(["C1", "total", col, sum(C1[c][j] for c in CLASSES)])
for a in CLASSES:
    for b in CLASSES:
        cmp_rows.append(["C2", f"old:{a}", f"new:{b}", C2[(a, b)]])
    cmp_rows.append(["C2", f"old:{a}", "total", sum(C2[(a, b)] for b in CLASSES)])
for b in CLASSES:
    cmp_rows.append(["C2", "total", f"new:{b}", sum(C2[(a, b)] for a in CLASSES)])
cmp_rows.append(["C2", "total", "total", 39813])
cmp_rows.append(["C2", "n_changed_rows", "value", n_changed])
for k in ["n_old", "n_new", "both", "old_only", "new_only"]:
    cmp_rows.append(["C3", "site_occ_response", k, C3[k]])
for fr in FRS + ["all"]:
    cmp_rows.append(["C4", fr, "n_old_both+occupancy_only", C4[fr][0]])
    cmp_rows.append(["C4", fr, "n_new_both+occupancy_only", C4[fr][1]])
S = [
    ("rows_log2occ_CTRL_nonempty", "old_site_traj", occ_ctrl_old),
    ("rows_log2occ_CTRL_nonempty", "old_replicated_step0", occ_ctrl_rep),
    ("rows_log2occ_CTRL_nonempty", "new", occ_ctrl_new),
    ("obs_total", "value", n_obs_total),
    ("obs_occ_missing_due_to_denom_protein_missing", "new", n_obs_no_dprot),
    ("obs_entering_occ", "old_replicated_step0", n_obs_occ_old),
    ("obs_entering_occ", "new", n_obs_occ_new),
    ("log2int_cells_nonempty_either", "new_vs_site_traj", int_cells),
    ("log2int_cells_string_diff", "new_vs_site_traj", int_str_diff),
    ("log2int_cells_presence_diff", "new_vs_site_traj", int_pres_diff),
    ("log2int_cells_numeric_diff", "new_vs_site_traj", int_num_diff),
    ("log2int_max_abs_diff", "new_vs_site_traj", f"{int_maxd / 1000:.3f}"),
    ("log2int_cells_string_diff", "new_vs_old_replicated_step0", int_vs_rep_diff),
] + [(f"class_count_{c}", "new_with_a0_uncalibrated", new_a0_cls.count(c)) for c in CLASSES] + [
    ("rows_class_diff", "new_a0_vs_new_calibrated", sum(1 for a, b in zip(new_a0_cls, new_cls) if a != b)),
]
for r, c, v in S:
    cmp_rows.append(["S", r, c, v])
pd.DataFrame(cmp_rows, columns=["table", "row", "column", "value"]).to_csv(
    O_CMP, sep="\t", index=False, lineterminator="\n")


# ------------------------------------------------------------------ 01_class_compare.md
def md_table(header, rows):
    out = ["| " + " | ".join(map(str, header)) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(map(str, r)) + " |" for r in rows]
    return "\n".join(out)


def cmp_line(name, r):
    return [name, f"{r['class']}/39813", f"{r['tier']}/39813", f"{r['peak']}/39813",
            f"{r['rows_presence_ok']}/39813", r["cells_presence_diff"], r["cells_nonempty"],
            f"{r['maxdiff_thou'] / 1000:.3f}", r["cells_diff_eq_0.001"], r["cells_diff_gt_0.001"]]


def drange(sol):
    v = [d for s in sol.values() for d in s["delta"].values()]
    return min(v), max(v)


dp_min, dp_max = drange(sol_p)
dq_min, dq_max = drange(sol_q)
md = []
md.append("# 01 统一蛋白分母：旧/新分类比较\n")
md.append("## ① 输入\n")
md.append(f"- `data/pilot/occupancy_variants_long.tsv`：{len(L)} 行（不含表头）；取 `site_id, timepoint, fraction, rep, "
          f"phos_intensity, denom_residue, denom_protein`。`phos_intensity` 空 0 行；`denom_residue` 空 "
          f"{int((L['denom_residue'] == '').sum())} 行；`denom_protein` 空 {n_obs_no_dprot} 行；(site_id, fraction) "
          f"{len(traj)} 个，与 site_traj 集合相同。")
md.append(f"- `data/pilot/run_medians.tsv`：{len(RM)} 行（phospho 120、proteome 120）；layer = `run` 路径前缀；"
          f"(timepoint, fraction, rep) 用 07 `parse_design` 解析，240/240 成功。")
md.append(f"- `data/pilot/site_traj.tsv`：{len(T)} 行、16 列；旧分类与步骤 0 校准靶。")
md.append("- 读表：`pd.read_csv(path, sep='\\t', dtype=str, keep_default_na=False)`，空串 = 缺失。")
md.append("- 脚本：`analysis/07_pilot2/01_build_v2.py`；输出：`01_run_factors.tsv`、`01_site_traj_v2.tsv`、"
          "`01_class_compare.tsv`、`01_class_compare.md`（均在 `analysis/07_pilot2/`）。\n")

md.append("## ② 定义\n")
md.append("### 07 occupancy 规则（07 `:238-306`，逐行核对后照做）")
md.append("1. run 因子（log2 加性形式 `a_run`）：`num = phos_intensity × 2^a_phospho_run`，`den = denom × 2^a_proteome_run`。")
md.append("2. 轨迹 = (site_id, fraction)；观测 = 长表该轨迹下全部 (timepoint, rep) 行（`phos_intensity` 均非空）。")
md.append("3. 旧档位：`n_res` = `denom_residue` 非空的观测数；`n_res >= 0.5 * len(obs)` → residue 档（d = den_res），否则 protein 档"
          "（d = den_prot）；d 缺失的观测不进 occ（07 `if d:`）。")
md.append("4. `inten[tp] += log2(num)`；`occ[tp] += log2(num / d)`。")
md.append("5. responder：`len(CTRL) < 2` → 非响应；`base = median(CTRL)`；按 2min→8min→20min→90min，跳过 `len < 2` 的 tp，"
          "取第一个满足 `abs(median − base) >= 1.0` 且所有重复 `(x − base) * fc > 0` 的 tp；`peak_tp = o_tp or i_tp`"
          "（occupancy 达标时取其 tp，否则取 intensity 的 tp）。")
md.append("6. class：o 与 i 皆响应 → both；仅 o → occupancy_only；仅 i → intensity_only；否则 none。`log2occ_tp` / `log2int_tp` = "
          "该 tp 值数 ≥ 2 时 `round(median, 3)`（写出格式同 07 csv.writer，即 `str(round(x, 3))`），否则空。")
md.append("- `median()` 用 07 `:117-118` 原函数（奇数取中间值，偶数取两中间值均值；与 `statistics.median` 同义）。\n")
md.append("### 步骤 0 校准方法")
md.append(f"- 初值 `a0_run = anchor_layer − median_log2_run`，`anchor_layer = statistics.median(该 layer 120 个 median_log2)`："
          f"phospho anchor = {anchor['phospho']:.5f}，proteome anchor = {anchor['proteome']:.5f}。")
md.append("- phospho：取 site_traj `log2int_tp` 非空且长表该 (site, fraction, tp) 恰 2 个重复 {i, j} 的组，"
          "`diff = log2int_tp(07) − median(log2(phos_i) + a0_i, log2(phos_j) + a0_j)`，模型 `diff = (δ_i + δ_j)/2`；"
          "每个 (fraction, tp) cell 4 个未知数（Rep1–Rep4），`numpy.linalg.lstsq` 最小二乘；`a = a0 + δ`。")
md.append("- proteome：按规则 3 定每条观测的分母，取 site_traj `log2occ_tp` 非空且恰 2 个重复有分母的组，"
          "`diff_occ = log2occ_tp(07) − median(log2(phos_i) + a_i^phos − log2(den_i) − a0_i^prot, …)`，模型 "
          "`diff_occ = −(δ_i^prot + δ_j^prot)/2`，同法解（`a^phos` 用已校准值）。")
md.append(f"- 每个 cell 的设计矩阵秩：rank < 4 的 cell 数 = {N_CELLS_LOW_RANK}（60 个 cell 全部 rank = 4），"
          "因此未使用 3 重复组补充，也没有 run 保持 a0。")
md.append("- `01_run_factors.tsv` 的 `n_equations`、`max_residual` 为该 run 所在 cell 的值（同 cell 4 个 run 相同）；"
          "`max_residual = max|diff − 模型值|`。a0/delta/a 写 6 位小数，脚本内部用全精度。")
md.append("- 与 site_traj 比对：log2 值比较 `round(本脚本值, 3)` 与 07 表值，差以 0.001 为单位（整数千分位）计。\n")
md.append("### 步骤 1（新）")
md.append("- 规则 2、4、5、6 不变；规则 3 改为所有观测 `d = denom_protein × 2^a_prot`（a 为校准后值），`denom_protein` 空的观测不进 occ；"
          "`denom_tier` 固定写 `protein`。`protein` 列 = site_id 去掉 `_残基位置` 后缀（与 site_traj 逐行相同）。\n")
md.append("### 步骤 2（比较）")
md.append("- 旧 = `data/pilot/site_traj.tsv`；旧复现 = 步骤 0 校准后重算；新 = `01_site_traj_v2.tsv`。")
md.append("- 位点级 occupancy 响应 = 该 site_id 任一 fraction 行 class ∈ {both, occupancy_only}。\n")

md.append("## ③ 校准结果\n")
md.append("### 与 site_traj 的一致性")
md.append(md_table(["因子", "class 一致", "denom_tier 一致", "peak_tp 一致", "10 个 log2 列有无全一致的行",
                    "有无不一致单元格", "比较的非空单元格", "最大差", "差 = 0.001 的单元格", "差 > 0.001 的单元格"],
                   [cmp_line("a0（run_medians）", cmp_a0),
                    cmp_line("a0 + δ^phos（proteome 仍 a0）", cmp_p_only),
                    cmp_line("a = a0 + δ（phospho + proteome）", cmp_cal)]))
md.append("")
if mism_cal:
    md.append(f"校准后 class/denom_tier/peak_tp 任一不一致的行：{len(mism_cal)} 行。")
    rowsm = []
    for i in mism_cal:
        x = rows_cal[i]
        fr = TV["fraction"][i]
        base = median(x["occ"]["CTRL"]) if len(x["occ"].get("CTRL", [])) >= 2 else None
        fcs = []
        for tp in EGF_TPS:
            v = x["occ"].get(tp, [])
            fcs.append(f"{median(v) - base:.6f}" if base is not None and len(v) >= 2 else "")
        rowsm.append([TV["site_id"][i], fr, TV["class"][i], x["cls"], TV["denom_tier"][i], x["tier"],
                      TV["peak_tp"][i] or "(空)", x["peak"] or "(空)"] + fcs)
    md.append(md_table(["site_id", "fraction", "class 07", "class 本脚本", "tier 07", "tier 本脚本",
                        "peak_tp 07", "peak_tp 本脚本", "本脚本 occ fc 2min", "fc 8min", "fc 20min", "fc 90min"], rowsm))
    md.append("")
    md.append("（fc = median(occ_tp) − median(occ_CTRL)，用校准后 a，未取整；阈值 FC_THRESH = 1.0。）\n")
else:
    md.append("校准后 class/denom_tier/peak_tp 三列全部一致。\n")

md.append("### δ 范围")
md.append(md_table(["layer", "δ min", "δ max", "cell 数", "每 cell 方程数 min", "每 cell 方程数 max", "每 run 方程数 min",
                    "max_residual 最大"],
                   [["phospho", f"{dp_min:.6f}", f"{dp_max:.6f}", len(sol_p), min(s["n"] for s in sol_p.values()),
                     max(s["n"] for s in sol_p.values()),
                     min(s["per_run"][r] for s in sol_p.values() for r in REPS),
                     f"{max(s['maxres'] for s in sol_p.values()):.6f}"],
                    ["proteome", f"{dq_min:.6f}", f"{dq_max:.6f}", len(sol_q), min(s["n"] for s in sol_q.values()),
                     max(s["n"] for s in sol_q.values()),
                     min(s["per_run"][r] for s in sol_q.values() for r in REPS),
                     f"{max(s['maxres'] for s in sol_q.values()):.6f}"]]))
md.append("")
md.append("### 每个 cell 的方程数、秩与最大残差")
cell_rows = []
for fr in FRS:
    for tp in TPS:
        sp, sq = sol_p[(fr, tp)], sol_q[(fr, tp)]
        cell_rows.append([fr, tp, sp["n"], sp["rank"], f"{sp['maxres']:.6f}",
                          " / ".join(f"{sp['delta'][r]:.5f}" for r in REPS),
                          sq["n"], sq["rank"], f"{sq['maxres']:.6f}",
                          " / ".join(f"{sq['delta'][r]:.5f}" for r in REPS)])
md.append(md_table(["fraction", "timepoint", "phospho n_eq", "rank", "max_res", "δ Rep1/2/3/4",
                    "proteome n_eq", "rank", "max_res", "δ Rep1/2/3/4"], cell_rows))
md.append("")

md.append("## ④ 比较表\n")
md.append("### C1 四分类计数")
md.append(md_table(["class", "n_old (site_traj)", "n_old_replicated (步骤0, 校准后 a)", "n_new"],
                   [[c] + list(C1[c]) for c in CLASSES]
                   + [["合计"] + [sum(C1[c][j] for c in CLASSES) for j in range(3)]]))
md.append("")
md.append("### C2 交叉表（行 = 旧 site_traj class，列 = 新 class）")
md.append(md_table(["旧 \\ 新"] + CLASSES + ["行合计"],
                   [[a] + [C2[(a, b)] for b in CLASSES] + [sum(C2[(a, b)] for b in CLASSES)] for a in CLASSES]
                   + [["列合计"] + [sum(C2[(a, b)] for a in CLASSES) for b in CLASSES] + [39813]]))
md.append("")
md.append(f"`n_changed_rows`（旧 class ≠ 新 class）= {n_changed}\n")
md.append("### C3 位点级 occupancy 响应（任一 fraction 为 both/occupancy_only）")
md.append(md_table(["n_old", "n_new", "两者皆", "仅旧", "仅新"],
                   [[C3["n_old"], C3["n_new"], C3["both"], C3["old_only"], C3["new_only"]]]))
md.append("")
md.append("### C4 按 fraction 的 both + occupancy_only 行数")
md.append(md_table(["fraction", "旧 (site_traj)", "新"], [[fr, C4[fr][0], C4[fr][1]] for fr in FRS + ["all"]]))
md.append("")
md.append("### 补充计数")
md.append(md_table(["项目", "口径", "值"], [[r, c, v] for r, c, v in S]))
md.append("")
md.append("### 自检")
md.append(md_table(["检查", "结果"], [[n, "通过" if ok else "未通过"] for n, ok in checks]))
md.append("")

md.append("## ⑤ 跳过/缺失项\n")
md.append("- 3 重复组补充：60 个 cell 均 rank = 4，未触发，未使用。")
md.append("- 保持 a0 的 run：0 个。")
md.append("- `data/pilot/sites.tsv`、`summary*.txt`：缺，本任务不依赖，跳过。")
with open(O_MD, "w") as f:
    f.write("\n".join(md) + "\n")

# ------------------------------------------------------------------ 控制台摘要
print("a0  :", dict(cmp_a0))
print("p   :", dict(cmp_p_only))
print("cal :", dict(cmp_cal))
print("delta phospho", dp_min, dp_max, "proteome", dq_min, dq_max)
print("C1", C1)
print("C3", C3, "changed", n_changed)
print("log2int new vs site_traj: cells", int_cells, "strdiff", int_str_diff, "presdiff", int_pres_diff,
      "numdiff", int_num_diff, "maxd", int_maxd / 1000, "| new vs rep strdiff", int_vs_rep_diff)
for n, ok in checks:
    print(("PASS " if ok else "FAIL ") + n)
