#!/usr/bin/env python3
"""07_pilot2 任务 4(b)：ERK 打分（kinase-library，ser_thr）—— 只用独立 venv 的 python 跑。

按 analysis/07_pilot2/00_design.md "任务 4 (b)" 执行：
  - 位点 = data/pilot/site_traj.tsv 的 18268 个 site_id（ACC_S123 → acc, 残基, 位置）；
  - FASTA 读法照 code/Protein contour/Zhihan/Test/05_build_occupancy.py:40-57（load_fasta）；
  - 15-mer 取法照 code/Protein contour/Phospho/martinez_network_check/kl_core.py:55-65（cut_window, half=7, '_' 补齐）；
  - 打分：优先批量接口 kl.PhosphoProteomics(...).predict(kin_type='ser_thr')（库默认：全部 ser_thr 激酶），
    取 ERK1 / ERK2 的 score 与 percentile（库自带、库默认取整位数）；批量接口失败才逐条 kl.Substrate。
输入 : data/pilot/site_traj.tsv, data/pilot/human_2026-08-19_UP2026_02.fasta
输出 : analysis/07_pilot2/04_erk_sites.tsv（第一阶段：13 列；resp_*/early_* 四列由 04_time_rules_v2.py 补上并重写）
       <scratchpad>/t4_erk_score_log.json（接口、版本、耗时、成功数等，供 04_time_rules_v2.py 写进 md）
复跑 : <scratchpad>/venv_kl/bin/python analysis/07_pilot2/04_erk_score.py
       然后  python3 analysis/07_pilot2/04_time_rules_v2.py
"""
import csv
import json
import math
import os
import re
import sys
import time
import traceback

SP = "/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad"
VENV = os.path.join(SP, "venv_kl")
# ---- 环境硬约束：必须在独立 venv 里跑，不碰全局 Python
assert sys.prefix != sys.base_prefix, f"不在 venv 中（sys.prefix={sys.prefix}），拒绝运行"
assert os.path.realpath(sys.prefix) == os.path.realpath(VENV), f"venv 不是 {VENV}（sys.prefix={sys.prefix}）"

import pandas as pd  # noqa: E402  （venv 内的 pandas）

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data", "pilot")
P_TRAJ = os.path.join(DATA, "site_traj.tsv")
P_FASTA = os.path.join(DATA, "human_2026-08-19_UP2026_02.fasta")
OUT_SITES = os.path.join(HERE, "04_erk_sites.tsv")
OUT_LOG = os.path.join(SP, "t4_erk_score_log.json")

HALF = 7
KIN_TYPE = "ser_thr"
ERKS = ["ERK1", "ERK2"]          # kinase_name_mapping.csv：MAPK3→ERK1，MAPK1→ERK2
P95 = 95.0
TOP_FRAC = 0.05
STAGE1_COLS = ["site_id", "protein", "residue", "position", "window15", "center_ok", "plus1_P",
               "erk1_score", "erk1_percentile", "erk2_score", "erk2_percentile",
               "erk_candidate_p95", "erk_top5pct_set"]
SITE_RE = re.compile(r"^(.+)_([A-Z])(\d+)$")


def load_fasta(path):
    """照 05_build_occupancy.py:40-57。"""
    seqs = {}
    acc = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                m = re.match(r">\w+\|([^|]+)\|", line)
                acc = m.group(1) if m else line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def cut_window(seq, pos, half=HALF):
    """照 kl_core.py:55-65：以 1-based pos 为中心的 (2*half+1)-mer，两端 '_' 补齐。"""
    if seq is None:
        return None, "NO_SEQ"
    if pos < 1 or pos > len(seq):
        return None, "POS_OOR"
    center = seq[pos - 1]
    left = seq[max(0, pos - 1 - half):pos - 1]
    right = seq[pos:pos + half]
    left = "_" * (half - len(left)) + left
    right = right + "_" * (half - len(right))
    return left + center + right, center


def tf(b):
    return "True" if b else "False"


def fmt(x, nd):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    s = f"{x:.{nd}f}"
    return s[1:] if s.startswith("-") and float(s) == 0 else s


t_all = time.time()
log = {"python": sys.executable, "sys_prefix": sys.prefix, "base_prefix": sys.base_prefix,
       "pandas_in_venv": pd.__version__}

# ================================================================ 位点与窗口
site_prot = {}
with open(P_TRAJ, newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        site_prot.setdefault(r["site_id"], r["protein"])
fasta = load_fasta(P_FASTA)
log["n_fasta_entries"] = len(fasta)

recs = []
n_parse_fail = 0
for sid, prot in site_prot.items():
    m = SITE_RE.match(sid)
    if not m:
        n_parse_fail += 1
        recs.append({"site_id": sid, "protein": prot, "residue": "", "position": "", "window15": "",
                     "_status": "PARSE_FAIL"})
        continue
    acc, res, pos = m.group(1), m.group(2), int(m.group(3))
    win, st = cut_window(fasta.get(acc), pos)
    rec = {"site_id": sid, "protein": prot, "residue": res, "position": str(pos),
           "window15": win or "", "_acc": acc,
           "_status": "OK" if win else st}
    if win:
        rec["_center_ok"] = (win[HALF] == res)
        rec["_plus1_P"] = (win[HALF + 1] == "P")
    recs.append(rec)

n_sites = len(recs)
n_win = sum(1 for r in recs if r["window15"])
n_center_ok = sum(1 for r in recs if r.get("_center_ok") is True)
n_center_bad = sum(1 for r in recs if r.get("_center_ok") is False)
log.update({
    "n_sites": n_sites, "n_parse_fail": n_parse_fail, "n_window": n_win,
    "n_no_window_by_status": {k: sum(1 for r in recs if r["_status"] == k)
                              for k in sorted({r["_status"] for r in recs}) if k != "OK"},
    "n_acc_not_in_fasta": sum(1 for r in recs if r["_status"] == "NO_SEQ"),
    "n_center_ok": n_center_ok, "n_center_mismatch": n_center_bad,
    "center_mismatch_examples": [(r["site_id"], r["window15"]) for r in recs if r.get("_center_ok") is False][:10],
    "n_plus1_P": sum(1 for r in recs if r.get("_plus1_P") is True),
    "residue_counts": {k: sum(1 for r in recs if r["residue"] == k) for k in sorted({r["residue"] for r in recs})},
})

# ================================================================ 打分
to_score = [r for r in recs if r.get("_center_ok") is True]   # 中心残基不一致的排除
log["n_submitted"] = len(to_score)
log["n_submitted_by_residue"] = {k: sum(1 for r in to_score if r["residue"] == k)
                                 for k in sorted({r["residue"] for r in to_score})}
scores = {}   # site_id -> dict(erk1_score, erk1_percentile, erk2_score, erk2_percentile)
log["attempts"] = []

try:
    import kinase_library as kl
    log["kinase_library_version"] = getattr(kl, "__version__", "")
    log["kinase_library_file"] = kl.__file__
    st_list = list(kl.get_kinase_list(kin_type=KIN_TYPE))
    log["n_ser_thr_kinases"] = len(st_list)
    log["has_ERK1"], log["has_ERK2"] = "ERK1" in st_list, "ERK2" in st_list
    import_ok = True
except Exception as e:  # noqa: BLE001
    import_ok = False
    log["import_error"] = f"{type(e).__name__}: {e}"
    log["import_traceback"] = traceback.format_exc()

if import_ok:
    # ---- 首选：批量接口
    t0 = time.time()
    att = {"interface": "kl.PhosphoProteomics(df, seq_col='window15').predict(kin_type='ser_thr')",
           "kwargs_note": "predict 全部取库默认值（kinases=None → 全部 ser_thr 激酶；st_fav=True；"
                          "score_round_digits=3；percentile_round_digits=2）；PhosphoProteomics 默认 pp=False、"
                          "drop_invalid_subs=True"}
    try:
        df_in = pd.DataFrame({"site_id": [r["site_id"] for r in to_score],
                              "window15": [r["window15"] for r in to_score]})
        t_init = time.time()
        pps = kl.PhosphoProteomics(df_in, seq_col="window15")
        att["init_s"] = round(time.time() - t_init, 2)
        att["n_input"] = len(df_in)
        att["n_omitted_invalid"] = int(len(pps.omited_entries))
        att["n_ser_thr_data"] = int(len(pps.ser_thr_data))
        att["n_tyrosine_data"] = int(len(pps.tyrosine_data))
        att["omitted_examples"] = (pps.omited_entries.head(10).astype(str).values.tolist()
                                   if len(pps.omited_entries) else [])
        t_pred = time.time()
        pred = pps.predict(kin_type=KIN_TYPE)
        att["predict_s"] = round(time.time() - t_pred, 2)
        att["pred_shape"] = list(pred.shape)
        att["n_kinases_in_pred"] = sum(1 for c in pred.columns if c.endswith("_percentile_rank"))
        for kin in ERKS:
            for met in ("score", "percentile"):
                assert f"{kin}_{met}" in pred.columns, f"预测表缺列 {kin}_{met}"
        n_nan = 0
        for row in pred[["site_id"] + [f"{k}_{m}" for k in ERKS for m in ("score", "percentile")]].itertuples(index=False):
            vals = [float(v) for v in row[1:]]
            if any(math.isnan(v) for v in vals):
                n_nan += 1
                continue
            scores[row[0]] = {"erk1_score": vals[0], "erk1_percentile": vals[1],
                              "erk2_score": vals[2], "erk2_percentile": vals[3]}
        att["n_pred_rows_with_nan"] = n_nan
        att["n_scored"] = len(scores)
        att["status"] = "成功"
        att["elapsed_s"] = round(time.time() - t0, 2)
        log["interface_used"] = att["interface"]
        log["score_round_digits"], log["percentile_round_digits"] = 3, 2
        # ---- 抽查：同一批前 20 个已打分窗口，逐条 kl.Substrate(seq, kin_type='ser_thr').predict() 的 ERK 值差
        try:
            chk_ids = [r["site_id"] for r in to_score if r["site_id"] in scores][:20]
            wmap = {r["site_id"]: r["window15"] for r in to_score}
            dmax = {"score": 0.0, "percentile": 0.0}
            t_s = time.time()
            for sid in chk_ids:
                sp = kl.Substrate(wmap[sid], kin_type=KIN_TYPE).predict()
                for kin, key in (("ERK1", "erk1"), ("ERK2", "erk2")):
                    dmax["score"] = max(dmax["score"], abs(float(sp.loc[kin, "Score"]) - scores[sid][f"{key}_score"]))
                    dmax["percentile"] = max(dmax["percentile"],
                                             abs(float(sp.loc[kin, "Percentile"]) - scores[sid][f"{key}_percentile"]))
            log["substrate_crosscheck"] = {"n": len(chk_ids), "max_abs_diff_score": round(dmax["score"], 6),
                                           "max_abs_diff_percentile": round(dmax["percentile"], 6),
                                           "elapsed_s": round(time.time() - t_s, 2),
                                           "call": "kl.Substrate(window15, kin_type='ser_thr').predict()（库默认 "
                                                   "score_round_digits=4、percentile_round_digits=2）"}
        except Exception as e:  # noqa: BLE001
            log["substrate_crosscheck"] = {"error": f"{type(e).__name__}: {e}"}
    except Exception as e:  # noqa: BLE001
        att["status"] = "失败"
        att["error"] = f"{type(e).__name__}: {e}"
        att["traceback"] = traceback.format_exc()
        att["elapsed_s"] = round(time.time() - t0, 2)
        scores = {}
    log["attempts"].append(att)

    # ---- 回退：逐条 Substrate（仅在批量接口失败时）
    if not scores:
        t0 = time.time()
        att = {"interface": "kl.Substrate(window15, kin_type='ser_thr').predict()（逐条）"}
        errs = {}
        for r in to_score:
            try:
                sp = kl.Substrate(r["window15"], kin_type=KIN_TYPE).predict()
                scores[r["site_id"]] = {"erk1_score": float(sp.loc["ERK1", "Score"]),
                                        "erk1_percentile": float(sp.loc["ERK1", "Percentile"]),
                                        "erk2_score": float(sp.loc["ERK2", "Score"]),
                                        "erk2_percentile": float(sp.loc["ERK2", "Percentile"])}
            except Exception as e:  # noqa: BLE001
                k = f"{type(e).__name__}: {str(e)[:150]}"
                errs[k] = errs.get(k, 0) + 1
        att["n_scored"] = len(scores)
        att["errors"] = errs
        att["elapsed_s"] = round(time.time() - t0, 2)
        att["status"] = "成功" if scores else "失败"
        log["attempts"].append(att)
        if scores:
            log["interface_used"] = att["interface"]
            log["score_round_digits"], log["percentile_round_digits"] = 4, 2

log["n_scored"] = len(scores)
log["n_scored_by_residue"] = {k: sum(1 for r in recs if r["site_id"] in scores and r["residue"] == k)
                              for k in sorted({r["residue"] for r in recs})}
log["n_submitted_not_scored_by_residue"] = {
    k: sum(1 for r in to_score if r["site_id"] not in scores and r["residue"] == k)
    for k in sorted({r["residue"] for r in to_score})}

# ================================================================ 定义 1 与 前 5% 对照
cand = {}
for sid, v in scores.items():
    cand[sid] = max(v["erk1_percentile"], v["erk2_percentile"]) >= P95
top5 = {}
if scores:
    mx = {sid: max(v["erk1_score"], v["erk2_score"]) for sid, v in scores.items()}
    vals = sorted(mx.values(), reverse=True)
    k = int(math.floor(TOP_FRAC * len(vals)))
    thr = vals[k - 1] if k >= 1 else float("inf")
    for sid, x in mx.items():
        top5[sid] = x >= thr
    log["top5"] = {"n_scored": len(vals), "k_floor_5pct": k, "threshold_max_score": thr,
                   "n_in_set": sum(top5.values()),
                   "n_ties_at_threshold": sum(1 for x in vals if x == thr)}
log["n_erk1_p95"] = sum(1 for v in scores.values() if v["erk1_percentile"] >= P95)
log["n_erk2_p95"] = sum(1 for v in scores.values() if v["erk2_percentile"] >= P95)
log["n_both_p95"] = sum(1 for v in scores.values() if v["erk1_percentile"] >= P95 and v["erk2_percentile"] >= P95)
log["n_candidate_p95_union"] = sum(cand.values())

# ================================================================ 写第一阶段 tsv（13 列）
sd = log.get("score_round_digits", 4)
pdg = log.get("percentile_round_digits", 2)
out = []
for r in recs:
    sid = r["site_id"]
    v = scores.get(sid)
    out.append({
        "site_id": sid, "protein": r["protein"], "residue": r["residue"], "position": r["position"],
        "window15": r["window15"],
        "center_ok": tf(r["_center_ok"]) if "_center_ok" in r else "",
        "plus1_P": tf(r["_plus1_P"]) if "_plus1_P" in r else "",
        "erk1_score": fmt(v["erk1_score"], sd) if v else "",
        "erk1_percentile": fmt(v["erk1_percentile"], pdg) if v else "",
        "erk2_score": fmt(v["erk2_score"], sd) if v else "",
        "erk2_percentile": fmt(v["erk2_percentile"], pdg) if v else "",
        "erk_candidate_p95": tf(cand[sid]) if sid in cand else "",
        "erk_top5pct_set": tf(top5[sid]) if sid in top5 else "",
    })
with open(OUT_SITES, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=STAGE1_COLS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(out)

log["elapsed_total_s"] = round(time.time() - t_all, 2)
log["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
with open(OUT_LOG, "w", encoding="utf-8") as f:
    json.dump(log, f, ensure_ascii=False, indent=1, default=str)

print(json.dumps({k: v for k, v in log.items() if k not in ("attempts",)}, ensure_ascii=False, default=str)[:3000])
for a in log["attempts"]:
    print({k: v for k, v in a.items() if k != "traceback"})
