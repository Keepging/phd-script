#!/usr/bin/env python3
"""07_pilot2 任务 4：时间规律重跑（旧/新分类并列）+ ERK 阳性对照比例表。用系统 python3 跑。

按 analysis/07_pilot2/00_design.md "任务 4" 与 analysis/06_pilot/00_design.md "任务 3"（3a/3b 定义沿用）执行；
3a/3b 逐行逻辑复用 analysis/06_pilot/03_time_pilots.py（同一算法，旧列须与上一轮完全相同）。
输入 : data/pilot/site_traj.tsv（旧）, analysis/07_pilot2/01_site_traj_v2.tsv（新）,
       analysis/07_pilot2/04_erk_sites.tsv（04_erk_score.py 产出的第一阶段 13 列）,
       analysis/06_pilot/03_time_pilots_rows.tsv（仅用于自检：旧逐行表须与之逐格相同）,
       data/pilot/kinase_name_mapping.csv（只列 ERK 映射）,
       <scratchpad>/t4_ver_before.txt、t4_venv_install.log、t4_erk_score_log.json（环境与打分记录）
输出 : analysis/07_pilot2/04_time_rules_v2_rows.tsv（新 R 逐行表）,
       analysis/07_pilot2/04_erk_sites.tsv（补 resp_old/resp_new/early_old/early_new 四列后重写为 17 列）,
       analysis/07_pilot2/04_time_rules_v2.md
复跑 : <scratchpad>/venv_kl/bin/python analysis/07_pilot2/04_erk_score.py   （先）
       python3 analysis/07_pilot2/04_time_rules_v2.py                       （后；可单独重复跑，结果不变）
"""
import json
import os
import subprocess
import sys
import time
from collections import Counter

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data", "pilot")
SP = "/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad"
VENV = os.path.join(SP, "venv_kl")
P_OLD = os.path.join(DATA, "site_traj.tsv")
P_NEW = os.path.join(HERE, "01_site_traj_v2.tsv")
P_FASTA = os.path.join(DATA, "human_2026-08-19_UP2026_02.fasta")
P_MAP = os.path.join(DATA, "kinase_name_mapping.csv")
P_PREV_ROWS = os.path.join(ROOT, "analysis", "06_pilot", "03_time_pilots_rows.tsv")
P_ERK = os.path.join(HERE, "04_erk_sites.tsv")
P_VER_BEFORE = os.path.join(SP, "t4_ver_before.txt")
P_INSTALL_LOG = os.path.join(SP, "t4_venv_install.log")
P_SCORE_LOG = os.path.join(SP, "t4_erk_score_log.json")
OUT_ROWS = os.path.join(HERE, "04_time_rules_v2_rows.tsv")
OUT_MD = os.path.join(HERE, "04_time_rules_v2.md")
VER_CMD = "import pandas,numpy,matplotlib;print(pandas.__version__,numpy.__version__,matplotlib.__version__)"

FC_THRESH = 1.0
MIN_REPS = 2
EGF_TPS = ["2min", "8min", "20min", "90min"]
FRS = [f"FR{i}" for i in range(1, 7)]
RESP_CLASSES = {"both", "occupancy_only"}
COMP = {"FR1": "Cyt", "FR2": "Cyt", "FR3": "Mem", "FR4": "Mem", "FR5": "Nuc", "FR6": "Nuc"}
T90 = ["transient", "sustained", "90min缺失"]
SITE_TYPES = ["transient", "sustained", "mixed", "untyped"]
EARLY = {"2min", "8min"}
ROW_COLS = ["site_id", "protein", "fraction", "class", "denom_tier", "peak_tp_07", "peak_tp_max",
            "fc_2min", "fc_8min", "fc_20min", "fc_90min", "fc_at_peak", "type_90"]
STAGE1_COLS = ["site_id", "protein", "residue", "position", "window15", "center_ok", "plus1_P",
               "erk1_score", "erk1_percentile", "erk2_score", "erk2_percentile",
               "erk_candidate_p95", "erk_top5pct_set"]
ERK_COLS = STAGE1_COLS + ["resp_old", "resp_new", "early_old", "early_new"]

# 上一轮（analysis/06_pilot/03_time_pilots.md）的数：旧列必须完全相同
PREV = {
    "A1": {"2min": [704, 255, 316], "8min": [768, 428, 261], "20min": [513, 298, 224], "90min": [0, 1882, 0]},
    "site_types": [1114, 1683, 538, 452],
    "nfr": [2547, 814, 426], "nfr_3to6": [283, 99, 35, 9],
    "ncr": [2874, 738, 175],
    "R": (5649, 3787, 1837),
}
NEW_EXPECT = (5419, 3711)   # 任务 1 的 C3 n_new 与 C4 all（新 R 行数 / 位点数）
ADVISOR_B = {"n_window": 18268, "n_center_bad": 0, "n_plus1P": 7148, "n_plus1P_resp_old": 2300}


def rel(p):
    return os.path.relpath(p, ROOT)


def read_tsv(p):
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)


def md_table(header, rows):
    out = ["| " + " | ".join(str(h) for h in header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(out)


def fmt3(x):
    s = f"{x:.3f}"
    return "0.000" if s == "-0.000" else s


def ratio(a, b):
    return f"{a / b:.4f}" if b else ""


def tf(b):
    return "True" if b else "False"


checks = []  # (项, 实际, 期望, 通过)


def check(desc, actual, expected):
    checks.append((desc, actual, expected, actual == expected))


# ================================================================ 3a / 3b（同 06_pilot 03_time_pilots.py 的算法）
def analyze(traj):
    R = traj[traj["class"].isin(RESP_CLASSES)].copy()
    assert (R["log2occ_CTRL"] != "").all(), "R 中存在 log2occ_CTRL 为空的行"
    assert set(R["peak_tp"]) <= set(EGF_TPS)
    out = []
    for rd in R.to_dict("records"):
        ctrl = rd["log2occ_CTRL"]
        fc = {tp: float(rd[f"log2occ_{tp}"]) - float(ctrl) for tp in EGF_TPS if rd[f"log2occ_{tp}"] != ""}
        assert fc, f"{rd['site_id']} {rd['fraction']} 无任何 EGF fc"
        peak = None
        for tp in EGF_TPS:   # |fc| 最大，并列取最早
            if tp in fc and (peak is None or abs(fc[tp]) > abs(fc[peak])):
                peak = tp
        if rd["log2occ_90min"] == "":
            t90 = "90min缺失"
        elif abs(fc["90min"]) < FC_THRESH:
            t90 = "transient"
        else:
            t90 = "sustained"
        out.append({"site_id": rd["site_id"], "protein": rd["protein"], "fraction": rd["fraction"],
                    "class": rd["class"], "denom_tier": rd["denom_tier"],
                    "peak_tp_07": rd["peak_tp"], "peak_tp_max": peak,
                    **{f"fc_{tp}": (fmt3(fc[tp]) if tp in fc else "") for tp in EGF_TPS},
                    "fc_at_peak": fmt3(fc[peak]), "type_90": t90})
    rows = pd.DataFrame(out, columns=ROW_COLS)
    A1 = Counter(zip(rows["peak_tp_max"], rows["type_90"]))
    site_type = {}
    for sid, g in rows.groupby("site_id", sort=False):
        ts = set(g["type_90"]) - {"90min缺失"}
        site_type[sid] = "untyped" if not ts else "mixed" if len(ts) == 2 else ts.pop()
    fr_resp = rows.groupby("site_id", sort=False)["fraction"].apply(set).to_dict()
    nfr = Counter(len(v) for v in fr_resp.values())
    ncr = Counter(len({COMP[f] for f in v}) for v in fr_resp.values())
    early_site = rows.assign(_e=rows["peak_tp_max"].isin(EARLY)).groupby("site_id", sort=False)["_e"].any().to_dict()
    return {"R": R, "rows": rows, "A1": A1, "site_type": site_type, "site_type_cnt": Counter(site_type.values()),
            "nfr": nfr, "ncr": ncr, "early_site": early_site,
            "n_rows": len(R), "n_sites": R["site_id"].nunique(), "n_prot": R["protein"].nunique(),
            "n_90_empty": int((R["log2occ_90min"] == "").sum()),
            "class_cnt": R["class"].value_counts().to_dict()}


old_t = read_tsv(P_OLD)
new_t = read_tsv(P_NEW)
for d in (old_t, new_t):
    assert d.isna().sum().sum() == 0
    assert not d.duplicated(["site_id", "fraction"]).any()
assert list(old_t.columns) == list(new_t.columns)
same_keys = (old_t[["site_id", "fraction"]].values.tolist() == new_t[["site_id", "fraction"]].values.tolist())
OLD = analyze(old_t)
NEW = analyze(new_t)

# 旧逐行表须与上一轮 03_time_pilots_rows.tsv 逐格相同
prev_rows = read_tsv(P_PREV_ROWS)
old_rows_eq_prev = (list(prev_rows.columns) == ROW_COLS and prev_rows.shape == OLD["rows"].shape
                    and (prev_rows.values == OLD["rows"][ROW_COLS].values).all())

NEW["rows"][ROW_COLS].to_csv(OUT_ROWS, sep="\t", index=False)
rows_back = read_tsv(OUT_ROWS)


def a1_list(X, tp):
    return [X["A1"].get((tp, t), 0) for t in T90]


def nfr_list(X):
    return [X["nfr"][1], X["nfr"][2], sum(X["nfr"][k] for k in range(3, 7))]


# ================================================================ (b) ERK：补四列并重写 04_erk_sites.tsv
erk_in = read_tsv(P_ERK)
missing_stage1 = [c for c in STAGE1_COLS if c not in erk_in.columns]
assert not missing_stage1, f"04_erk_sites.tsv 缺列 {missing_stage1}（先跑 04_erk_score.py）"
erk = erk_in[STAGE1_COLS].copy()
assert erk["site_id"].is_unique
resp_old_set, resp_new_set = set(OLD["early_site"]), set(NEW["early_site"])
erk["resp_old"] = [tf(s in resp_old_set) for s in erk["site_id"]]
erk["resp_new"] = [tf(s in resp_new_set) for s in erk["site_id"]]
erk["early_old"] = [tf(OLD["early_site"][s]) if s in resp_old_set else "" for s in erk["site_id"]]
erk["early_new"] = [tf(NEW["early_site"][s]) if s in resp_new_set else "" for s in erk["site_id"]]
erk[ERK_COLS].to_csv(P_ERK, sep="\t", index=False)
erk_back = read_tsv(P_ERK)

# ---- (b) 统计（全部从写出的 tsv 上数）
eb = erk_back
is_T = lambda col: eb[col] == "True"  # noqa: E731
scored = eb["erk1_score"] != ""
B = {
    "n_sites": len(eb),
    "n_window": int((eb["window15"] != "").sum()),
    "n_window_len15": int((eb["window15"].str.len() == 15).sum()),
    "n_center_ok": int(is_T("center_ok").sum()),
    "n_center_bad": int((eb["center_ok"] == "False").sum()),
    "n_pad": int(eb["window15"].str.contains("_").sum()),
    "n_plus1P": int(is_T("plus1_P").sum()),
    "n_plus1_pad": int((eb["window15"].str[8] == "_").sum()),
    "n_scored": int(scored.sum()),
    "n_unscored": int((~scored).sum()),
    "res_all": eb["residue"].value_counts().to_dict(),
    "res_scored": eb.loc[scored, "residue"].value_counts().to_dict(),
    "res_unscored": eb.loc[~scored, "residue"].value_counts().to_dict(),
    "n_erk1_p95": int((scored & (pd.to_numeric(eb["erk1_percentile"], errors="coerce") >= 95)).sum()),
    "n_erk2_p95": int((scored & (pd.to_numeric(eb["erk2_percentile"], errors="coerce") >= 95)).sum()),
    "n_cand": int(is_T("erk_candidate_p95").sum()),
    "n_top5": int(is_T("erk_top5pct_set").sum()),
    "n_cand_and_top5": int((is_T("erk_candidate_p95") & is_T("erk_top5pct_set")).sum()),
    "n_cand_and_P": int((is_T("erk_candidate_p95") & is_T("plus1_P")).sum()),
    "n_top5_and_P": int((is_T("erk_top5pct_set") & is_T("plus1_P")).sum()),
    "n_plus1P_scored": int((is_T("plus1_P") & scored).sum()),
}
p1 = pd.to_numeric(eb["erk1_percentile"], errors="coerce")
p2 = pd.to_numeric(eb["erk2_percentile"], errors="coerce")
B["n_both_p95"] = int(((p1 >= 95) & (p2 >= 95)).sum())
B["n_only_erk1"] = int(((p1 >= 95) & ~(p2 >= 95)).sum())
B["n_only_erk2"] = int((~(p1 >= 95) & (p2 >= 95)).sum())
for lab in ("old", "new"):
    rs = is_T(f"resp_{lab}")
    B[f"resp_{lab}"] = int(rs.sum())
    B[f"resp_{lab}_unscored"] = int((rs & ~scored).sum())
    B[f"resp_{lab}_unscored_res"] = eb.loc[rs & ~scored, "residue"].value_counts().to_dict()
    B[f"resp_{lab}_plus1P"] = int((rs & is_T("plus1_P")).sum())
    B[f"resp_{lab}_cand"] = int((rs & is_T("erk_candidate_p95")).sum())
    B[f"resp_{lab}_erk1"] = int((rs & (p1 >= 95)).sum())
    B[f"resp_{lab}_erk2"] = int((rs & (p2 >= 95)).sum())
    B[f"resp_{lab}_top5"] = int((rs & is_T("erk_top5pct_set")).sum())


# ---- 表 E1
E1_ROWS = [("全部响应位点", None), ("ERK 候选(定义1)", "erk_candidate_p95"),
           ("+1P(定义2)", "plus1_P"), ("集合前5%(对照)", "erk_top5pct_set")]


def e1_site(df):
    """位点口径：从 04_erk_sites.tsv 过滤。"""
    res = {}
    for lab in ("old", "new"):
        base = df[df[f"resp_{lab}"] == "True"]
        for name, col in E1_ROWS:
            sub = base if col is None else base[base[col] == "True"]
            res[(name, lab)] = (len(sub), int((sub[f"early_{lab}"] == "True").sum()))
    return res


def e1_row(rows_by_lab, df):
    """行口径：响应行（逐行表）按 site_id 连上 04_erk_sites.tsv 的位点属性；早峰 = 该行 peak_tp_max ∈ {2min, 8min}。"""
    attr = df.set_index("site_id")
    res = {}
    for lab, rws in rows_by_lab.items():
        j = rws[["site_id", "peak_tp_max"]].join(attr[["erk_candidate_p95", "plus1_P", "erk_top5pct_set"]], on="site_id")
        assert j[["erk_candidate_p95", "plus1_P", "erk_top5pct_set"]].notna().all().all()
        for name, col in E1_ROWS:
            sub = j if col is None else j[j[col] == "True"]
            res[(name, lab)] = (len(sub), int(sub["peak_tp_max"].isin(EARLY).sum()))
    return res


E1S = e1_site(erk)
E1R = e1_row({"old": OLD["rows"], "new": NEW["rows"]}, erk)
# 复算：位点口径从写出的 04_erk_sites.tsv；行口径用写出的逐行表（新：04_time_rules_v2_rows.tsv；旧：06 的 03_time_pilots_rows.tsv）
E1S_re = e1_site(erk_back)
E1R_re = e1_row({"old": prev_rows, "new": rows_back}, erk_back)


def e1_table(E):
    out = []
    for name, _ in E1_ROWS:
        on, oe = E[(name, "old")]
        nn, ne = E[(name, "new")]
        out.append([name, on, oe, ratio(oe, on), nn, ne, ratio(ne, nn)])
    return out


E1_HDR = ["", "旧 n", "旧 早峰 n", "旧 比例", "新 n", "新 早峰 n", "新 比例"]

# ================================================================ 环境记录
ver_before_txt = open(P_VER_BEFORE, encoding="utf-8").read().strip().splitlines() if os.path.exists(P_VER_BEFORE) else []
install_lines = open(P_INSTALL_LOG, encoding="utf-8").read().splitlines() if os.path.exists(P_INSTALL_LOG) else []
score_log = json.load(open(P_SCORE_LOG, encoding="utf-8")) if os.path.exists(P_SCORE_LOG) else None
pyvenv_cfg = open(os.path.join(VENV, "pyvenv.cfg"), encoding="utf-8").read().splitlines() if os.path.exists(os.path.join(VENV, "pyvenv.cfg")) else []
kmap = pd.read_csv(P_MAP, sep=",", dtype=str, keep_default_na=False)
erk_map = kmap[kmap["our_name"].isin(["MAPK1", "MAPK3"])]


def grab(prefix):
    return next((ln for ln in install_lines if ln.startswith(prefix)), "")


inst_installed = next((ln for ln in install_lines if ln.startswith("Successfully installed")), "")
inst_kl = [w for w in inst_installed.split() if w.startswith(("kinase-library", "kinase_library", "pandas-", "numpy-", "matplotlib-"))]

# ---- 结束后版本检查（全部计算与 tsv 写出之后；系统 python3）
cp = subprocess.run(["python3", "-c", VER_CMD], capture_output=True, text=True)
ver_after = cp.stdout.strip()
ver_after_time = time.strftime("%Y-%m-%dT%H:%M:%S+00:00", time.gmtime())
ver_before = ver_before_txt[1] if len(ver_before_txt) > 1 else ""

# ================================================================ 自检
check("旧 R 行/位点/蛋白", (OLD["n_rows"], OLD["n_sites"], OLD["n_prot"]), PREV["R"])
for tp in EGF_TPS:
    check(f"旧 A1 {tp} 行（transient/sustained/90min缺失）= 上一轮", a1_list(OLD, tp), PREV["A1"][tp])
check("旧 位点类型 transient/sustained/mixed/untyped = 上一轮", [OLD["site_type_cnt"][t] for t in SITE_TYPES], PREV["site_types"])
check("旧 n_frac_resp 1/2/3+ = 上一轮", nfr_list(OLD), PREV["nfr"])
check("旧 n_frac_resp 3/4/5/6 = 上一轮", [OLD["nfr"][k] for k in range(3, 7)], PREV["nfr_3to6"])
check("旧 n_comp_resp 1/2/3 = 上一轮", [OLD["ncr"][k] for k in (1, 2, 3)], PREV["ncr"])
check("旧逐行表与 06_pilot/03_time_pilots_rows.tsv 逐格相同", bool(old_rows_eq_prev), True)
check("新旧表 (site_id, fraction) 顺序相同", same_keys, True)
check("新 R 行数 / 位点数 = 任务 1（C4 all / C3 n_new）", (NEW["n_rows"], NEW["n_sites"]), NEW_EXPECT)
check("新 A1 合计 = 新 R 行数", sum(NEW["A1"].values()), NEW["n_rows"])
check("新 A1 90min缺失 列合计 = 新 R 中 log2occ_90min 为空行数",
      sum(NEW["A1"].get((tp, "90min缺失"), 0) for tp in EGF_TPS), NEW["n_90_empty"])
check("新 位点类型四类合计 = 新 R 位点数", sum(NEW["site_type_cnt"][t] for t in SITE_TYPES), NEW["n_sites"])
check("新 n_frac_resp 1/2/3+ 合计 = 新 R 位点数", sum(nfr_list(NEW)), NEW["n_sites"])
check("新 n_comp_resp 1/2/3 合计 = 新 R 位点数", sum(NEW["ncr"][k] for k in (1, 2, 3)), NEW["n_sites"])
check("新 peak_tp_max=90min 且 transient 行数", NEW["A1"].get(("90min", "transient"), 0), 0)
check("04_time_rules_v2_rows.tsv 行数 / 列名", (len(rows_back), list(rows_back.columns) == ROW_COLS), (NEW["n_rows"], True))
check("04_time_rules_v2_rows.tsv 中 nan/NaN/None/NA 文本个数", int(rows_back.isin(["nan", "NaN", "None", "NA"]).sum().sum()), 0)
check("04_erk_sites.tsv 行数 / 列名严格", (len(erk_back), list(erk_back.columns) == ERK_COLS), (18268, True))
check("04_erk_sites.tsv site_id 集合 = site_traj 位点集合", set(erk_back["site_id"]) == set(old_t["site_id"]), True)
check("18268 位点全有 15-mer 窗口", (B["n_window"], B["n_window_len15"]), (18268, 18268))
check("中心残基不一致数（advisor 参考 0）", B["n_center_bad"], ADVISOR_B["n_center_bad"])
check("+1 位是 P 位点数（advisor 参考 7148）", B["n_plus1P"], ADVISOR_B["n_plus1P"])
check("旧响应位点中 +1P 数（advisor 参考 2300）", B["resp_old_plus1P"], ADVISOR_B["n_plus1P_resp_old"])
check("resp_old / resp_new True 数 = 旧/新 R 位点数", (B["resp_old"], B["resp_new"]), (OLD["n_sites"], NEW["n_sites"]))
check("04_erk_sites.tsv 布尔列只含 True/False/空串",
      all(set(erk_back[c]) <= {"True", "False", ""} for c in
          ["center_ok", "plus1_P", "erk_candidate_p95", "erk_top5pct_set", "resp_old", "resp_new", "early_old", "early_new"]), True)
check("04_erk_sites.tsv 中 nan/NaN/None/NA 文本个数", int(erk_back.isin(["nan", "NaN", "None", "NA"]).sum().sum()), 0)
check("定义1 候选数 = max(pct_ERK1, pct_ERK2) ≥ 95 复算", B["n_cand"], int(((p1 >= 95) | (p2 >= 95)).sum()))
if score_log:
    check("打分日志与 tsv 一致：成功数 / ERK1≥95 / ERK2≥95 / 并集 / 前5%",
          (score_log.get("n_scored"), score_log.get("n_erk1_p95"), score_log.get("n_erk2_p95"),
           score_log.get("n_candidate_p95_union"), (score_log.get("top5") or {}).get("n_in_set")),
          (B["n_scored"], B["n_erk1_p95"], B["n_erk2_p95"], B["n_cand"], B["n_top5"]))
check("E1 位点口径 n 列由 04_erk_sites.tsv 过滤复算一致", E1S_re == E1S, True)
check("E1 行口径 n 列由逐行表 + 04_erk_sites.tsv 复算一致", E1R_re == E1R, True)
check("E1 位点口径'全部响应位点' n = 旧/新 R 位点数",
      (E1S[("全部响应位点", "old")][0], E1S[("全部响应位点", "new")][0]), (OLD["n_sites"], NEW["n_sites"]))
check("E1 行口径'全部响应位点' n = 旧/新 R 行数",
      (E1R[("全部响应位点", "old")][0], E1R[("全部响应位点", "new")][0]), (OLD["n_rows"], NEW["n_rows"]))
check("E1 旧口径早峰数 = 上一轮基线（行 2732 / 位点 2179）",
      (E1R[("全部响应位点", "old")][1], E1S[("全部响应位点", "old")][1]), (2732, 2179))
check("全局 pandas/numpy/matplotlib 版本 开始前 = 结束后", ver_before == ver_after and ver_after != "", True)

# ================================================================ md
L = []
A = L.append
A("# 04 时间规律重跑（旧/新分类并列）+ ERK 阳性对照（任务 4）\n")

A("## ① 输入与环境\n")
A("### 输入\n")
A(md_table(["文件", "行数（不含表头）", "说明"], [
    [f"`{rel(P_OLD)}`", len(old_t), f"旧分类；{old_t['site_id'].nunique()} 位点、{old_t['protein'].nunique()} 蛋白"],
    [f"`{rel(P_NEW)}`", len(new_t), f"新分类（任务 1，蛋白档分母）；(site_id, fraction) 顺序与旧表相同：{same_keys}"],
    [f"`{rel(P_FASTA)}`", f"{score_log['n_fasta_entries'] if score_log else ''} 条序列", "05 `load_fasta`（`05_build_occupancy.py:40-57`）读法"],
    [f"`{rel(P_MAP)}`", len(kmap), "ERK 映射：" + "；".join(f"{r.our_name}→{r.atlas_name}（{r.kin_type}，{r.status}）" for r in erk_map.itertuples())],
    [f"`{rel(P_PREV_ROWS)}`", len(prev_rows), "仅用于自检（旧逐行表须逐格相同）及 E1 行口径旧列复算"],
]))
A("")
A("- 读表：`pd.read_csv(path, sep='\\t', dtype=str, keep_default_na=False)`，空串 = 缺失。阈值 FC_THRESH = 1.0，MIN_REPS = 2（未改）。")
A(f"- 旧 R（class ∈ {{both, occupancy_only}}）：{OLD['n_rows']} 行、{OLD['n_sites']} 位点、{OLD['n_prot']} 蛋白；"
  f"class：" + "，".join(f"{k} {v}" for k, v in sorted(OLD['class_cnt'].items())) + f"；`log2occ_90min` 为空 {OLD['n_90_empty']} 行。")
A(f"- 新 R：{NEW['n_rows']} 行、{NEW['n_sites']} 位点、{NEW['n_prot']} 蛋白；"
  f"class：" + "，".join(f"{k} {v}" for k, v in sorted(NEW['class_cnt'].items())) + f"；`log2occ_90min` 为空 {NEW['n_90_empty']} 行。\n")
A("### 环境\n")
A(md_table(["项", "记录"], [
    ["全局版本检查命令", f"`python3 -c \"{VER_CMD}\"`（`python3` = `/usr/local/bin/python3`）"],
    ["开始前（任务开始、建 venv 之前）", f"`{ver_before}`（{ver_before_txt[0] if ver_before_txt else '缺记录'}）"],
    ["结束后（本脚本全部计算与 tsv 写出之后，子进程执行）", f"`{ver_after}`（{ver_after_time}）"],
    ["前后是否相同", "相同" if ver_before == ver_after and ver_after else "**不同**"],
    ["venv 路径", f"`{VENV}`（`python3 -m venv` 创建；pyvenv.cfg：" + "；".join(f"`{x.strip()}`" for x in pyvenv_cfg if x.strip()) + "）"],
    ["venv 安装命令", f"`{VENV}/bin/pip install {SP}/kl/kinase_library-1.8.0-py3-none-any.whl`（联网装依赖，限时 900 s）"],
    ["venv 安装结果", f"{grab('pip_rc=')}；{grab('pip_end')}；venv 创建 {grab('venv_create_start')}"],
    ["venv 内相关包", "、".join(f"`{x}`" for x in inst_kl) if inst_kl else "（见安装日志）"],
    ["kinase-library 版本 / 导入位置", (f"{score_log.get('kinase_library_version')} / `{score_log.get('kinase_library_file')}`" if score_log else "缺打分日志")],
    ["打分进程", (f"`{score_log.get('python')}`（sys.prefix = venv，base_prefix = `{score_log.get('base_prefix')}`；venv 内 pandas {score_log.get('pandas_in_venv')}）" if score_log else "")],
    ["本脚本进程", f"`{sys.executable}`（系统 python3，pandas {pd.__version__}）"],
    ["全局环境备注", "全局 site-packages 中原已有 kinase-library 1.8.0（上一轮 `analysis/06_pilot/03_time_pilots.py` 用全局 pip 安装，见 06 的 md 第 186-190 行）；本轮未调用、未改动它"],
    ["安装/打分日志", f"`{P_INSTALL_LOG}`、`{P_SCORE_LOG}`"],
]))
A("")

A("## ② 定义\n")
A("- **R**：class ∈ {both, occupancy_only} 的行；旧 = `data/pilot/site_traj.tsv`，新 = `analysis/07_pilot2/01_site_traj_v2.tsv`。")
A("- **3a**（沿用 06_pilot 任务 3）：`fc_tp = float(log2occ_tp) − float(log2occ_CTRL)`（tp ∈ 2min/8min/20min/90min 且非空）；"
  "`peak_tp_max` = |fc| 最大的 tp（按 2→8→20→90 遍历、严格大于才替换，即并列取最早）；"
  "`type_90`：`log2occ_90min` 空 → `90min缺失`；|fc_90| < 1.0 → `transient`；否则 `sustained`。"
  "表 A1 = peak_tp_max × type_90 行数。位点级类型：忽略 `90min缺失` 行后，全 transient → transient，全 sustained → sustained，两者都有 → mixed，无可判行 → untyped。")
A("- **3b**：`n_frac_resp` = 位点 class ∈ {both, occupancy_only} 的 fraction 数；区室 FR1/FR2→Cyt、FR3/FR4→Mem、FR5/FR6→Nuc，`n_comp_resp` = 响应 fraction 覆盖的区室数。")
A("- **窗口**：site_id 按 `^(.+)_([A-Z])(\\d+)$` 解析为 acc、残基、位置；FASTA 按 05 `load_fasta`；"
  "15-mer 按 `kl_core.cut_window(acc, pos, half=7)`（中心 ±7，两端 `_` 补齐）；`center_ok` = 窗口第 8 位 == site_id 残基（不一致的不送打分）；"
  "`plus1_P` = 窗口第 9 位（+1 位）== `P`（C 端补齐的 `_` 记 False）。")
A("- **打分**：kinase-library ser_thr，实际接口与参数见 ④；`erk1_*`/`erk2_*` = 库输出的 `ERK1_score`、`ERK1_percentile`、`ERK2_score`、`ERK2_percentile`（库默认取整：score 3 位、percentile 2 位，tsv 原样写出）。")
A("- **定义 1（ERK 候选）**：`max(erk1_percentile, erk2_percentile) ≥ 95`；未打分位点该列为空串。另报 ERK1、ERK2 各自 ≥ 95 的数。")
A("- **定义 2（+1P）**：`plus1_P == True`。")
A("- **对照（集合前 5%）**：在全部打分成功的位点中按 `max(erk1_score, erk2_score)` 降序，k = floor(0.05 × 打分成功数)，"
  "取 `max_score ≥ 第 k 名的值`（边界并列全纳入）；未打分位点该列为空串。")
A("- **早峰**：响应行 `peak_tp_max ∈ {2min, 8min}`。位点口径：位点任一响应行早峰 → `early_* = True`；非响应位点 `early_*` 为空串。"
  "`resp_old`/`resp_new` = 位点在旧/新 R 中（True/False）。")
A("- **表 E1**：位点口径 = 在 `resp_* == True` 的位点上按行定义过滤（全部 / 定义1 True / plus1_P True / 前5% True），早峰 n = 其中 `early_* == True` 数；"
  "行口径 = 旧/新 R 的响应行按 site_id 带上位点属性后同样过滤，早峰 n = 该行 peak_tp_max ∈ {2min, 8min} 的行数。比例 = 早峰 n / n（4 位小数）。未打分位点只出现在“全部响应位点”行。")
A("- 只报计数与数值，不做检验。\n")

A("## ③ (a) 3a / 3b 旧新并列\n")
A("### 表 A1：peak_tp_max × type_90（行数）\n")
hdr = ["peak_tp_max"] + [f"旧 {t}" for t in T90] + ["旧 合计"] + [f"新 {t}" for t in T90] + ["新 合计"]
rows_a1 = []
for tp in EGF_TPS:
    o, n = a1_list(OLD, tp), a1_list(NEW, tp)
    rows_a1.append([tp] + o + [sum(o)] + n + [sum(n)])
ot = [sum(a1_list(OLD, tp)[i] for tp in EGF_TPS) for i in range(3)]
nt = [sum(a1_list(NEW, tp)[i] for tp in EGF_TPS) for i in range(3)]
rows_a1.append(["合计"] + ot + [sum(ot)] + nt + [sum(nt)])
A(md_table(hdr, rows_a1))
A("")
A("### 位点级类型（位点数）\n")
A(md_table(["site_type_90", "旧", "新"], [[t, OLD["site_type_cnt"][t], NEW["site_type_cnt"][t]] for t in SITE_TYPES]
           + [["合计", sum(OLD["site_type_cnt"].values()), sum(NEW["site_type_cnt"].values())]]))
A("")
A("### 3b：n_frac_resp（位点数）\n")
A(md_table(["n_frac_resp", "旧", "新"],
           [[1, OLD["nfr"][1], NEW["nfr"][1]], [2, OLD["nfr"][2], NEW["nfr"][2]], ["3+", nfr_list(OLD)[2], nfr_list(NEW)[2]]]
           + [[f"  其中 {k}", OLD["nfr"][k], NEW["nfr"][k]] for k in range(3, 7)]
           + [["合计", sum(nfr_list(OLD)), sum(nfr_list(NEW))]]))
A("")
A("### 3b：n_comp_resp（位点数）\n")
A(md_table(["n_comp_resp", "旧", "新"], [[k, OLD["ncr"][k], NEW["ncr"][k]] for k in (1, 2, 3)]
           + [["合计", sum(OLD["ncr"][k] for k in (1, 2, 3)), sum(NEW["ncr"][k] for k in (1, 2, 3))]]))
A("")
A(f"- 新逐行表 `{rel(OUT_ROWS)}`：{len(rows_back)} 行，列同 06_pilot `03_time_pilots_rows.tsv`（`denom_tier` 全为 "
  + "/".join(sorted(set(rows_back['denom_tier']))) + "）。")
A(f"- 旧列由同一算法对 site_traj 重算，逐行表与 `{rel(P_PREV_ROWS)}` 逐格相同：{old_rows_eq_prev}。\n")

A("## ④ (b) 窗口与打分统计\n")
att0 = (score_log or {}).get("attempts", [{}])
att_ok = next((a for a in att0 if a.get("status") == "成功"), att0[-1] if att0 else {})
A(md_table(["项", "值"], [
    ["site_traj 位点数", B["n_sites"]],
    ["取到 15-mer 窗口", f"{B['n_window']}（长度 15：{B['n_window_len15']}；含 `_` 补齐 {B['n_pad']}）"],
    ["acc 不在 FASTA / 位置越界", f"{(score_log or {}).get('n_acc_not_in_fasta', '')} / " + str(((score_log or {}).get('n_no_window_by_status') or {}).get('POS_OOR', 0))],
    ["中心残基一致 / 不一致", f"{B['n_center_ok']} / {B['n_center_bad']}"],
    ["残基构成（全部）", "，".join(f"{k} {v}" for k, v in sorted(B['res_all'].items()))],
    ["+1 位是 P（定义 2）", f"{B['n_plus1P']}（+1 位为 `_` 补齐：{B['n_plus1_pad']}）"],
    ["送打分（center_ok）", (score_log or {}).get("n_submitted", "")],
    ["实际接口", f"`{(score_log or {}).get('interface_used', '无（打分失败）')}`"],
    ["接口参数", att_ok.get("kwargs_note", "")],
    ["库内过滤", f"无效序列剔除 {att_ok.get('n_omitted_invalid', '')}；进入 ser_thr 集合 {att_ok.get('n_ser_thr_data', '')}；"
                f"中心为 Y、归入 tyrosine 集合（ser_thr 不打分）{att_ok.get('n_tyrosine_data', '')}"],
    ["ser_thr 激酶数（预测表）", f"{att_ok.get('n_kinases_in_pred', '')}（`get_kinase_list('ser_thr')` {(score_log or {}).get('n_ser_thr_kinases', '')}；含 ERK1 {(score_log or {}).get('has_ERK1')}、ERK2 {(score_log or {}).get('has_ERK2')}）"],
    ["耗时", f"PhosphoProteomics 初始化 {att_ok.get('init_s', '')} s；predict {att_ok.get('predict_s', '')} s；批量接口合计 {att_ok.get('elapsed_s', '')} s；04_erk_score.py 全程 {(score_log or {}).get('elapsed_total_s', '')} s"],
    ["打分成功（4 个 ERK 值均非空）", f"{B['n_scored']}（" + "，".join(f"{k} {v}" for k, v in sorted(B['res_scored'].items())) + "）"],
    ["未打分", f"{B['n_unscored']}（" + "，".join(f"{k} {v}" for k, v in sorted(B['res_unscored'].items())) + "）"],
    ["逐条接口抽查", (lambda c: f"前 {c.get('n')} 个打分位点用 {c.get('call')} 重算：abs(Δscore) 最大 {c.get('max_abs_diff_score')}，abs(Δpercentile) 最大 {c.get('max_abs_diff_percentile')}" if c and 'error' not in c else str(c))((score_log or {}).get("substrate_crosscheck"))],
    ["ERK1 percentile ≥ 95", B["n_erk1_p95"]],
    ["ERK2 percentile ≥ 95", B["n_erk2_p95"]],
    ["两者皆 ≥ 95 / 仅 ERK1 / 仅 ERK2", f"{B['n_both_p95']} / {B['n_only_erk1']} / {B['n_only_erk2']}"],
    ["ERK 候选（定义 1，并集）", B["n_cand"]],
    ["集合前 5%（对照）", f"{B['n_top5']}（k = {((score_log or {}).get('top5') or {}).get('k_floor_5pct', '')}，阈值 max_score ≥ {((score_log or {}).get('top5') or {}).get('threshold_max_score', '')}，阈值处并列 {((score_log or {}).get('top5') or {}).get('n_ties_at_threshold', '')}）"],
    ["交叠：定义1∩前5% / 定义1∩+1P / 前5%∩+1P", f"{B['n_cand_and_top5']} / {B['n_cand_and_P']} / {B['n_top5_and_P']}"],
    ["+1P 且打分成功", B["n_plus1P_scored"]],
]))
A("")
A("按响应集合（位点数）：\n")
A(md_table(["", "旧 R 位点", "新 R 位点"], [
    ["响应位点", B["resp_old"], B["resp_new"]],
    ["其中未打分（Y 中心）", f"{B['resp_old_unscored']}（" + "，".join(f"{k} {v}" for k, v in sorted(B['resp_old_unscored_res'].items())) + "）",
     f"{B['resp_new_unscored']}（" + "，".join(f"{k} {v}" for k, v in sorted(B['resp_new_unscored_res'].items())) + "）"],
    ["ERK1 ≥ 95", B["resp_old_erk1"], B["resp_new_erk1"]],
    ["ERK2 ≥ 95", B["resp_old_erk2"], B["resp_new_erk2"]],
    ["ERK 候选（定义 1）", B["resp_old_cand"], B["resp_new_cand"]],
    ["+1P（定义 2）", B["resp_old_plus1P"], B["resp_new_plus1P"]],
    ["集合前 5%（对照）", B["resp_old_top5"], B["resp_new_top5"]],
]))
A("")

A("## ⑤ 表 E1\n")
A("### E1-位点口径（位点数；早峰 = 位点任一响应行 peak_tp_max ∈ {2min, 8min}）\n")
A(md_table(E1_HDR, e1_table(E1S)))
A("")
A("### E1-行口径（响应行数；早峰 = 该行 peak_tp_max ∈ {2min, 8min}）\n")
A(md_table(E1_HDR, e1_table(E1R)))
A("")
A(f"- 复算：位点口径 n 列由写出的 `{rel(P_ERK)}` 过滤复算一致：{E1S_re == E1S}；"
  f"行口径由 `{rel(OUT_ROWS)}`（新）/ `{rel(P_PREV_ROWS)}`（旧）按 site_id 连 `{rel(P_ERK)}` 复算一致：{E1R_re == E1R}。\n")

A("## ⑥ 跳过 / 缺失项\n")
A(f"- 中心为 Y 的 {B['res_all'].get('Y', 0)} 个位点：kinase-library `PhosphoProteomics` 把它们归入 tyrosine 集合，`predict(kin_type='ser_thr')` 不对其打分；"
  "这些位点 `erk1_*`/`erk2_*`/`erk_candidate_p95`/`erk_top5pct_set` 为空串，只计入 E1 的“全部响应位点”行。")
A(f"- 中心残基不一致：{B['n_center_bad']} 个（无排除）。取不到窗口：{B['n_sites'] - B['n_window']} 个。")
A("- 方案写的逐条调用 `kl.Substrate(seq).predict(kin_type='ser_thr')` 与 1.8.0 的签名不符（`kin_type` 在 `Substrate(...)` 构造参数里，`predict()` 无此参数）；"
  "本轮用批量接口，未走逐条回退；抽查用 `kl.Substrate(seq, kin_type='ser_thr').predict()`。")
A("- 不做 backbone 对照、A2/A3、3b 交叉表（任务 4 未要求）；不做统计检验。")
A("- `sites.tsv`、`summary*.txt`：本任务不依赖，未读。\n")

A("### 附：验收自检（脚本自动核对）\n")
A(md_table(["项", "实际", "期望", "通过"], [[d, a, e, "是" if ok else "否"] for d, a, e, ok in checks]))
A("")
with open(OUT_MD, "w", encoding="utf-8") as f:
    f.write("\n".join(L))

# ================================================================ 控制台
print("OLD", OLD["n_rows"], OLD["n_sites"], OLD["n_prot"], "NEW", NEW["n_rows"], NEW["n_sites"], NEW["n_prot"])
for tp in EGF_TPS:
    print("A1", tp, a1_list(OLD, tp), a1_list(NEW, tp))
print("types old", [OLD["site_type_cnt"][t] for t in SITE_TYPES], "new", [NEW["site_type_cnt"][t] for t in SITE_TYPES])
print("nfr old", nfr_list(OLD), [OLD["nfr"][k] for k in range(3, 7)], "new", nfr_list(NEW), [NEW["nfr"][k] for k in range(3, 7)])
print("ncr old", [OLD["ncr"][k] for k in (1, 2, 3)], "new", [NEW["ncr"][k] for k in (1, 2, 3)])
print("B", {k: v for k, v in B.items()})
print("E1 site", e1_table(E1S))
print("E1 row", e1_table(E1R))
print("ver", ver_before, "->", ver_after)
nf = 0
for d, a, e, ok in checks:
    nf += not ok
    print("PASS" if ok else "FAIL", d, a, "expected", e)
print("FAILED", nf)
