#!/usr/bin/env python3
"""08_q1q3_cross.py — Q1(occupancy 解离) x Q3(轮廓动态) 位点级交叉。

输入 (07 产物, 不重扫 parquet):
  analysis07/site_traj.tsv, analysis07/site_profiles.tsv
类名重映射 (07 -> 08, 纯描述筛选结果):
  both -> occ_and_int | occupancy_only -> occ_not_int | intensity_only -> int_not_occ | none -> neither
聚合规则:
  位点级 occ_responder = 该位点任一组分轨迹的 occupancy 达筛选标准; int 同理
  profile_mover = 该位点任一 EGF 时间点的轮廓位移 >= 全体位移的第 PCT 百分位 (默认 90)
  divergence_mover = |位点vs母蛋白偏离度变化| >= 其全体第 PCT 百分位
输出: crosstab_summary.txt (仅计数/分类/阈值数值), top_candidates.tsv (数值清单)
用法: python3 08_q1q3_cross.py [analysis07目录] [PCT]
"""
import sys, os, csv
from collections import defaultdict

REMAP = {"both": "occ_and_int", "occupancy_only": "occ_not_int",
         "intensity_only": "int_not_occ", "none": "neither"}
TPS = ["CTRL", "2min", "8min", "20min", "90min"]
EGF_TPS = ["2min", "8min", "20min", "90min"]
FORBIDDEN = ["证明了假设", "验证了假设", "confirms the hypothesis", "发现了隐藏",
             "hidden signal", "breakthrough", "突破", "隐藏的信号"]

def pctl(sorted_vals, p):
    if not sorted_vals:
        return None
    k = (len(sorted_vals) - 1) * p / 100
    f = int(k)
    return sorted_vals[f] if f == k else sorted_vals[f] + (sorted_vals[f + 1] - sorted_vals[f]) * (k - f)

def rank_pct(sorted_vals, x):
    lo, hi = 0, len(sorted_vals)
    while lo < hi:
        mid = (lo + hi) // 2
        if sorted_vals[mid] <= x:
            lo = mid + 1
        else:
            hi = mid
    return lo / len(sorted_vals) * 100

def main():
    a7 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.environ["VSC_DATA"], "rerun_2026-08/occupancy/analysis07")
    PCT = float(sys.argv[2]) if len(sys.argv) > 2 else 90.0
    traj_p = os.path.join(a7, "site_traj.tsv")
    prof_p = os.path.join(a7, "site_profiles.tsv")
    for p in (traj_p, prof_p):
        if not os.path.isfile(p):
            print(f"FAIL: 缺输入 {p}"); sys.exit(1)

    # ---- Q1 位点级聚合 ----
    q1 = {}   # site -> dict
    with open(traj_p) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sid = r["site_id"]
            d = q1.setdefault(sid, {"protein": r["protein"], "occ": False, "int": False,
                                    "res_evid": False, "max_occ_fc": None, "occ_fc_tp": "",
                                    "occ_fc_fr": "", "max_int_fc": None})
            cls = REMAP.get(r["class"], r["class"])
            o = cls in ("occ_and_int", "occ_not_int")
            it = cls in ("occ_and_int", "int_not_occ")
            d["occ"] |= o
            d["int"] |= it
            if o and r["denom_tier"] == "residue":
                d["res_evid"] = True
            try:
                base = float(r["log2occ_CTRL"])
                for tp in EGF_TPS:
                    v = r.get(f"log2occ_{tp}", "")
                    if v != "":
                        fc = abs(float(v) - base)
                        if d["max_occ_fc"] is None or fc > d["max_occ_fc"]:
                            d["max_occ_fc"], d["occ_fc_tp"], d["occ_fc_fr"] = fc, tp, r["fraction"]
            except ValueError:
                pass
            try:
                base = float(r["log2int_CTRL"])
                for tp in EGF_TPS:
                    v = r.get(f"log2int_{tp}", "")
                    if v != "":
                        fc = abs(float(v) - base)
                        if d["max_int_fc"] is None or fc > d["max_int_fc"]:
                            d["max_int_fc"] = fc
            except ValueError:
                pass

    # ---- Q3 位点级聚合 ----
    q3 = {}
    all_shift, all_ddelta = [], []
    with open(prof_p) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sid = r["site_id"]
            d = q3.setdefault(sid, {"t3": r["proteome_tier"] == "T3", "max_shift": None,
                                    "shift_tp": "", "max_ddelta": None, "ddelta_tp": ""})
            if r["phos_shift_vs_ctrl"] != "":
                v = float(r["phos_shift_vs_ctrl"])
                all_shift.append(v)
                if d["max_shift"] is None or v > d["max_shift"]:
                    d["max_shift"], d["shift_tp"] = v, r["timepoint"]
            if r["divergence_delta_vs_ctrl"] != "":
                v = abs(float(r["divergence_delta_vs_ctrl"]))
                all_ddelta.append(v)
                if d["max_ddelta"] is None or v > d["max_ddelta"]:
                    d["max_ddelta"], d["ddelta_tp"] = v, r["timepoint"]

    all_shift.sort(); all_ddelta.sort()
    thr_shift = pctl(all_shift, PCT)
    thr_ddelta = pctl(all_ddelta, PCT)
    occ_fcs = sorted(d["max_occ_fc"] for d in q1.values() if d["max_occ_fc"] is not None)

    # ---- 交叉 ----
    def site_class(d):
        return ("occ_and_int" if d["occ"] and d["int"] else
                "occ_not_int" if d["occ"] else
                "int_not_occ" if d["int"] else "neither")

    cross = defaultdict(int)
    cand = []
    n_occ = n_occ_mov = n_occ_mov_res = n_occ_mov_dd = 0
    n_t3 = n_t3_mov = 0
    for sid, d in q1.items():
        p = q3.get(sid, {})
        mover = p.get("max_shift") is not None and thr_shift is not None and p["max_shift"] >= thr_shift
        dmover = p.get("max_ddelta") is not None and thr_ddelta is not None and p["max_ddelta"] >= thr_ddelta
        cls = site_class(d)
        cross[(cls, "mover" if mover else "non_mover")] += 1
        if p.get("t3"):
            n_t3 += 1
            if mover:
                n_t3_mov += 1
        if d["occ"]:
            n_occ += 1
            if mover:
                n_occ_mov += 1
                if d["res_evid"]:
                    n_occ_mov_res += 1
                if dmover:
                    n_occ_mov_dd += 1
                score = (rank_pct(occ_fcs, d["max_occ_fc"]) if d["max_occ_fc"] is not None else 0) \
                        + rank_pct(all_shift, p["max_shift"])
                cand.append([sid, d["protein"], site_class(d), int(d["res_evid"]),
                             round(d["max_occ_fc"], 3) if d["max_occ_fc"] is not None else "",
                             d["occ_fc_tp"], d["occ_fc_fr"],
                             round(p["max_shift"], 4), p["shift_tp"],
                             round(p["max_ddelta"], 4) if p.get("max_ddelta") is not None else "",
                             int(dmover), round(score, 1)])

    cand.sort(key=lambda r: -r[-1])
    with open(os.path.join(a7, "top_candidates.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["site_id", "protein", "site_class", "residue_evidence",
                    "max_abs_log2occ_fc", "occ_fc_tp", "occ_fc_fraction",
                    "max_profile_shift", "shift_tp", "max_abs_divergence_delta",
                    "divergence_ge_p{:g}".format(PCT), "rank_score"])
        w.writerows(cand)

    lines = [
        "类名映射: both->occ_and_int | occupancy_only->occ_not_int | intensity_only->int_not_occ | none->neither",
        f"阈值: profile_shift P{PCT:g} = {thr_shift:.4f} | |divergence_delta| P{PCT:g} = {thr_ddelta:.4f}",
        f"位点总数(有轨迹): {len(q1)} | 有轮廓数据: {len(q3)}",
        "",
        "== 交叉列联 (位点级) ==",
    ]
    for cls in ("occ_and_int", "occ_not_int", "int_not_occ", "neither"):
        m, nm = cross[(cls, "mover")], cross[(cls, "non_mover")]
        lines.append(f"  {cls}: mover {m} | non_mover {nm} | 合计 {m + nm}")
    lines += [
        "",
        "== 核心交叉子集 ==",
        f"occ 达标位点(位点级): {n_occ}",
        f"  其中 profile_shift >= P{PCT:g} 的: {n_occ_mov} (占 occ 达标位点的 {n_occ_mov / n_occ * 100:.1f}%)" if n_occ else "  (无)",
        f"  其中 residue 档位证据的: {n_occ_mov_res}",
        f"  其中 |divergence_delta| 同时 >= P{PCT:g} 的: {n_occ_mov_dd}",
        f"T3 档位点(有轨迹且 proteome 无该蛋白): {n_t3} | 其中 profile_shift >= P{PCT:g} 的: {n_t3_mov}",
        "",
        f"候选清单: top_candidates.tsv ({len(cand)} 行, 按 rank_score 降序)",
    ]
    txt = "\n".join(lines)
    for wbad in FORBIDDEN:
        assert wbad not in txt, f"禁用词出现: {wbad}"
    with open(os.path.join(a7, "crosstab_summary.txt"), "w") as f:
        f.write(txt + "\n")
    print(txt)

if __name__ == "__main__":
    main()
