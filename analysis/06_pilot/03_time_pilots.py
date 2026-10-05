#!/usr/bin/env python3
"""06_pilot 任务 3：时间规律小实验（3a 形状 / 3b 空间特异性 / 3c ERK 阳性对照）。

按 analysis/06_pilot/00_design.md 的 "0. 通用约定" 与 "任务 3" 执行；
另按协调方更正：site_features.tsv 的缺失值是字面 "NA"，当缺失处理（不参与 backbone 中位数，sites.tsv 留空串）。
输入 : data/pilot/site_traj.tsv, data/pilot/site_features.tsv, data/pilot/site_profiles.tsv（只列列名）,
       data/pilot/kinase_name_mapping.csv, kinase-library wheel（见 WHEEL）
输出 : analysis/06_pilot/03_time_pilots_rows.tsv, 03_time_pilots_sites.tsv, 03_time_pilots.md
复跑 : python3 analysis/06_pilot/03_time_pilots.py            （会执行 pip install <wheel>，已装则为 already satisfied）
       python3 analysis/06_pilot/03_time_pilots.py --skip-pip （仅调试用：不执行 pip install）
"""
import json
import os
import re
import statistics
import subprocess
import sys
import time
from collections import Counter
from decimal import Decimal

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data", "pilot")
P_TRAJ = os.path.join(DATA, "site_traj.tsv")
P_FEAT = os.path.join(DATA, "site_features.tsv")
P_PROF = os.path.join(DATA, "site_profiles.tsv")
P_MAP = os.path.join(DATA, "kinase_name_mapping.csv")
P_DESIGN = os.path.join(HERE, "00_design.md")
WHEEL = ("/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/"
         "scratchpad/kl/kinase_library-1.8.0-py3-none-any.whl")
KL_CORE_REF = "code/Protein contour/Phospho/martinez_network_check/kl_core.py:7-8,35-42"
OUT_ROWS = os.path.join(HERE, "03_time_pilots_rows.tsv")
OUT_SITES = os.path.join(HERE, "03_time_pilots_sites.tsv")
OUT_MD = os.path.join(HERE, "03_time_pilots.md")

# 脚本已有常量（07 :22-25；08 :45），不改。PCT 在任务 3 中无用到之处。
FC_THRESH = 1.0
MIN_REPS = 2
PCT = 90
EGF_TPS = ["2min", "8min", "20min", "90min"]
FRS = [f"FR{i}" for i in range(1, 7)]
RESP_CLASSES = {"both", "occupancy_only"}
COMP = {"FR1": "Cyt", "FR2": "Cyt", "FR3": "Mem", "FR4": "Mem", "FR5": "Nuc", "FR6": "Nuc"}
T90 = ["transient", "sustained", "90min缺失"]
SITE_TYPES = ["transient", "sustained", "mixed", "untyped"]
EARLY = {"2min", "8min"}
SKIP_SEQ = "缺序列，跳过"

SKIP_PIP = "--skip-pip" in sys.argv[1:]


def rel(p):
    return os.path.relpath(p, ROOT)


def read_tsv(p):
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)


def md_table(header, rows):
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(out)


def crosstab(pairs, row_keys, col_keys, row_name, col_name):
    """pairs: list of (row_value, col_value). 返回 (header, rows, cnt, row_tot, col_tot, total)。"""
    cnt = Counter(pairs)
    header = [f"{row_name} \\ {col_name}"] + list(col_keys) + ["合计"]
    rows, row_tot, col_tot = [], {}, {c: 0 for c in col_keys}
    for r in row_keys:
        vals = [cnt.get((r, c), 0) for c in col_keys]
        row_tot[r] = sum(vals)
        for c, v in zip(col_keys, vals):
            col_tot[c] += v
        rows.append([r] + vals + [row_tot[r]])
    total = sum(row_tot.values())
    rows.append(["合计"] + [col_tot[c] for c in col_keys] + [total])
    unlisted = sum(v for (r, c), v in cnt.items() if r not in row_keys or c not in col_keys)
    assert unlisted == 0, f"交叉表有未列出的键: {[k for k in cnt if k[0] not in row_keys or k[1] not in col_keys]}"
    return header, rows, cnt, row_tot, col_tot, total


def fmt3(x):
    s = f"{x:.3f}"
    return "0.000" if s == "-0.000" else s


checks = []  # (编号, 描述, 实际, 期望, 通过)


def check(no, desc, actual, expected):
    checks.append((no, desc, actual, expected, actual == expected))


# ================================================================ 读表
traj = read_tsv(P_TRAJ)
feat = read_tsv(P_FEAT)
prof = read_tsv(P_PROF)
kmap = pd.read_csv(P_MAP, sep=",", dtype=str, keep_default_na=False)
for d in (traj, feat, prof, kmap):
    assert d.isna().sum().sum() == 0
assert not traj.duplicated(["site_id", "fraction"]).any(), "site_id×fraction 不唯一"
assert (traj.groupby("site_id")["protein"].nunique() == 1).all(), "site_id 对应多个 protein"
assert set(traj["fraction"]) <= set(FRS)

base = {
    "traj_rows": len(traj), "traj_sites": traj["site_id"].nunique(), "traj_prot": traj["protein"].nunique(),
    "class": traj["class"].value_counts().to_dict(), "tier": traj["denom_tier"].value_counts().to_dict(),
}

R = traj[traj["class"].isin(RESP_CLASSES)].copy()
nR_rows, nR_sites, nR_prot = len(R), R["site_id"].nunique(), R["protein"].nunique()
n_ctrl_empty = int((R["log2occ_CTRL"] == "").sum())
n_90_empty = int((R["log2occ_90min"] == "").sum())
assert n_ctrl_empty == 0, "R 中存在 log2occ_CTRL 为空的行"
assert set(R["peak_tp"]) <= set(EGF_TPS), "R 中 07 peak_tp 不在 EGF_TPS 内"

# ================================================================ 3a 行级
row_out = []
prec_peak_diff = prec_type_diff = n_fc90_eq1_exact = n_tie_float = n_tie_dec = 0
for rd in R.to_dict("records"):
    ctrl = rd["log2occ_CTRL"]
    fc = {tp: float(rd[f"log2occ_{tp}"]) - float(ctrl) for tp in EGF_TPS if rd[f"log2occ_{tp}"] != ""}
    assert fc, f"{rd['site_id']} {rd['fraction']} 无任何 EGF fc"
    # |fc| 最大的 tp，并列取最早（EGF_TPS 顺序遍历，严格大于才替换）
    peak = None
    for tp in EGF_TPS:
        if tp in fc and (peak is None or abs(fc[tp]) > abs(fc[peak])):
            peak = tp
    if rd["log2occ_90min"] == "":
        t90 = "90min缺失"
    elif abs(fc["90min"]) < FC_THRESH:
        t90 = "transient"
    else:
        t90 = "sustained"
    # 精度核对（不改变结果）：同一定义用十进制精确减法重算，比较 peak_tp_max / type_90 是否一致
    fcd = {tp: Decimal(rd[f"log2occ_{tp}"]) - Decimal(ctrl) for tp in fc}
    peak_d = None
    for tp in EGF_TPS:
        if tp in fcd and (peak_d is None or abs(fcd[tp]) > abs(fcd[peak_d])):
            peak_d = tp
    t90_d = t90 if t90 == "90min缺失" else ("transient" if abs(fcd["90min"]) < Decimal(1) else "sustained")
    prec_peak_diff += peak_d != peak
    prec_type_diff += t90_d != t90
    n_fc90_eq1_exact += ("90min" in fcd and abs(fcd["90min"]) == Decimal(1))
    m = max(abs(v) for v in fc.values())
    n_tie_float += sum(abs(v) == m for v in fc.values()) > 1
    md_ = max(abs(v) for v in fcd.values())
    n_tie_dec += sum(abs(v) == md_ for v in fcd.values()) > 1
    row_out.append({
        "site_id": rd["site_id"], "protein": rd["protein"], "fraction": rd["fraction"],
        "class": rd["class"], "denom_tier": rd["denom_tier"],
        "peak_tp_07": rd["peak_tp"], "peak_tp_max": peak,
        **{f"fc_{tp}": (fmt3(fc[tp]) if tp in fc else "") for tp in EGF_TPS},
        "fc_at_peak": fmt3(fc[peak]), "type_90": t90,
        "_absmax": m, "_fc90": fc.get("90min"),
    })
rows = pd.DataFrame(row_out)
ROW_COLS = ["site_id", "protein", "fraction", "class", "denom_tier", "peak_tp_07", "peak_tp_max",
            "fc_2min", "fc_8min", "fc_20min", "fc_90min", "fc_at_peak", "type_90"]

# 表 A1 / A2 / A3
A1 = crosstab(list(zip(rows["peak_tp_max"], rows["type_90"])), EGF_TPS, T90, "peak_tp_max", "type_90")
A2 = crosstab(list(zip(rows["peak_tp_07"], rows["type_90"])), EGF_TPS, T90, "peak_tp_07", "type_90")
A3 = crosstab(list(zip(rows["peak_tp_max"], rows["peak_tp_07"])), EGF_TPS, EGF_TPS, "peak_tp_max", "peak_tp_07")
n_rows_max_lt1 = int((rows["_absmax"] < FC_THRESH).sum())

# ================================================================ 3a 位点类型
site_order = list(dict.fromkeys(R["site_id"]))
site_prot = R.drop_duplicates("site_id").set_index("site_id")["protein"]
site_type = {}
for sid, g in rows.groupby("site_id", sort=False):
    ts = set(g["type_90"]) - {"90min缺失"}
    site_type[sid] = ("untyped" if not ts else "mixed" if len(ts) == 2 else ts.pop())
site_type_cnt = Counter(site_type.values())

# ================================================================ 3a backbone 对照
feat_i = feat.set_index("site_id")
assert feat_i.index.is_unique
missing_feat_rows = [s for s in site_order if s not in feat_i.index]
assert not missing_feat_rows, f"R 位点在 site_features 中缺行: {missing_feat_rows[:5]}"


def to_num(s):
    """'NA' / '' / 非数值 → None（缺失）。"""
    if s in ("", "NA"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


bb_res = {}
for col in ("single_backbone", "win5_backbone"):
    val = {s: to_num(feat_i.at[s, col]) for s in site_order}
    miss_R = [s for s in site_order if val[s] is None]
    raw_na_R = sum(feat_i.at[s, col] == "NA" for s in site_order)
    # 只用有 backbone 值的 transient / sustained 位点
    by_prot = {}
    by_prot_all = {}   # 不剔除缺失时的资格（仅用于报告缺失对蛋白资格的影响）
    for s in site_order:
        t = site_type[s]
        if t not in ("transient", "sustained"):
            continue
        p = site_prot[s]
        by_prot_all.setdefault(p, {"transient": 0, "sustained": 0})[t] += 1
        if val[s] is None:
            continue
        by_prot.setdefault(p, {"transient": [], "sustained": []})[t].append(val[s])
    elig = {p: v for p, v in by_prot.items() if v["transient"] and v["sustained"]}
    elig_all = {p for p, v in by_prot_all.items() if v["transient"] and v["sustained"]}
    ds = {p: statistics.median(v["transient"]) - statistics.median(v["sustained"]) for p, v in elig.items()}
    dvals = list(ds.values())
    miss_in_elig_all = [s for s in miss_R if site_prot[s] in elig_all and site_type[s] in ("transient", "sustained")]
    bb_res[col] = {
        "n_prot": len(ds),
        "d_median": statistics.median(dvals) if dvals else None,
        "pos": sum(d > 0 for d in dvals), "neg": sum(d < 0 for d in dvals), "zero": sum(d == 0 for d in dvals),
        "zero_tol": sum(abs(d) < 1e-9 for d in dvals),
        "n_tr_sites": sum(len(v["transient"]) for v in elig.values()),
        "n_su_sites": sum(len(v["sustained"]) for v in elig.values()),
        "miss_R": miss_R, "raw_na_R": raw_na_R,
        "miss_types": Counter(site_type[s] for s in miss_R),
        "n_prot_if_no_drop": len(elig_all),
        "miss_in_elig_all": miss_in_elig_all,
        "lost_prot": sorted(elig_all - set(ds)),
        "val": val,
    }

# ================================================================ 3b 空间特异性（在全表上数 evaluable）
egf_occ = [f"log2occ_{tp}" for tp in EGF_TPS]
traj["_eval"] = (traj["log2occ_CTRL"] != "") & (traj[egf_occ] != "").any(axis=1)
traj["_resp"] = traj["class"].isin(RESP_CLASSES)
Rset = set(site_order)
sub = traj[traj["site_id"].isin(Rset)]
n_resp_not_eval = int((sub["_resp"] & ~sub["_eval"]).sum())
fr_resp = sub[sub["_resp"]].groupby("site_id")["fraction"].apply(set).to_dict()
fr_eval = sub[sub["_eval"]].groupby("site_id")["fraction"].apply(set).to_dict()
site_rec = []
for s in site_order:
    fr_r, fr_e = fr_resp.get(s, set()), fr_eval.get(s, set())
    site_rec.append({"site_id": s, "protein": site_prot[s],
                     "n_frac_resp": len(fr_r), "n_frac_evaluable": len(fr_e),
                     "n_comp_resp": len({COMP[f] for f in fr_r}),
                     "n_comp_evaluable": len({COMP[f] for f in fr_e})})
early_any = rows.assign(_e=rows["peak_tp_max"].isin(EARLY)).groupby("site_id")["_e"].any().to_dict()
for rec in site_rec:
    s = rec["site_id"]
    rec["site_type_90"] = site_type[s]
    for col in ("single_backbone", "win5_backbone"):
        raw = feat_i.at[s, col]
        rec[col] = "" if bb_res[col]["val"][s] is None else raw
    rec["peak_early_any"] = "True" if early_any[s] else "False"
sites = pd.DataFrame(site_rec)
SITE_COLS = ["site_id", "protein", "n_frac_resp", "n_frac_evaluable", "n_comp_resp", "n_comp_evaluable",
             "site_type_90", "single_backbone", "win5_backbone", "peak_early_any"]

nfr = Counter(sites["n_frac_resp"])
ncr = Counter(sites["n_comp_resp"])
XF = crosstab(list(zip(sites["n_frac_resp"], sites["n_frac_evaluable"])), range(1, 7), range(1, 7),
              "n_frac_resp", "n_frac_evaluable")
XC = crosstab(list(zip(sites["n_comp_resp"], sites["n_comp_evaluable"])), range(1, 4), range(1, 4),
              "n_comp_resp", "n_comp_evaluable")

# ================================================================ 3c 序列来源检查
SEQ_PAT = re.compile(r"seq|window|flank|motif|peptide|fasta|aa_", re.I)
col_lists = {rel(P_FEAT): list(feat.columns),
             rel(P_TRAJ): [c for c in traj.columns if not c.startswith("_")],
             rel(P_PROF): list(prof.columns)}
seq_cols = {k: [c for c in v if SEQ_PAT.search(c)] for k, v in col_lists.items()}
fasta_hits, biophys_hits = [], []
for dp, dns, fns in os.walk(ROOT):
    dns[:] = [d for d in dns if d != ".git"]
    for n in dns + fns:
        if n == "biophys_json" or "biophys_json" in n:
            biophys_hits.append(rel(os.path.join(dp, n)))
    for n in fns:
        if n.lower().endswith(".fasta"):
            fasta_hits.append(rel(os.path.join(dp, n)))
has_seq = bool(fasta_hits or biophys_hits or any(seq_cols.values()))

# ---- pip install（限时 10 分钟，失败只重试一次）
VER_PKGS = ["pandas", "numpy", "matplotlib"]


def env_versions():
    """当前环境中已安装的版本（子进程查询，反映 pip 之后的状态）。"""
    code = ("import json, importlib.metadata as m\n"
            f"pk = {VER_PKGS!r}\n"
            "out = {}\n"
            "for p in pk:\n"
            "    try: out[p] = m.version(p)\n"
            "    except Exception: out[p] = '未安装'\n"
            "print(json.dumps(out))\n")
    cp = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    return json.loads(cp.stdout.strip().splitlines()[-1])


ver_before = env_versions()
ver_inproc = {"pandas": pd.__version__}
pip_attempts = []
if SKIP_PIP:
    pip_attempts.append({"attempt": 0, "status": "未执行（--skip-pip 调试运行）"})
else:
    for attempt in (1, 2):
        t0 = time.time()
        try:
            cp = subprocess.run([sys.executable, "-m", "pip", "install", WHEEL],
                                capture_output=True, text=True, timeout=600)
            el = time.time() - t0
            tail = [ln for ln in (cp.stdout + "\n" + cp.stderr).splitlines() if ln.strip()][-3:]
            pip_attempts.append({"attempt": attempt, "returncode": cp.returncode,
                                 "elapsed_s": round(el, 1), "tail": tail,
                                 "status": "成功" if cp.returncode == 0 else "失败"})
            if cp.returncode == 0:
                break
        except subprocess.TimeoutExpired:
            pip_attempts.append({"attempt": attempt, "returncode": None, "elapsed_s": 600.0,
                                 "tail": ["TimeoutExpired (600 s)"], "status": "超时"})
pip_ok = any(a.get("returncode") == 0 for a in pip_attempts)
ver_after = env_versions()

# ---- import kinase_library + get_kinase_list（子进程执行，调用方式同 kl_core.py:7-8）
KL_SNIPPET = (
    "import json, sys\n"
    "out = {'import_ok': False}\n"
    "try:\n"
    "    import kinase_library as kl\n"
    "    out['import_ok'] = True\n"
    "    out['version'] = getattr(kl, '__version__', '')\n"
    "    st = list(kl.get_kinase_list(kin_type='ser_thr'))\n"
    "    out['n_ser_thr'] = len(st)\n"
    "    out['has_ERK1'] = 'ERK1' in st\n"
    "    out['has_ERK2'] = 'ERK2' in st\n"
    "except Exception as e:\n"
    "    out['error'] = f'{type(e).__name__}: {e}'\n"
    "print('KLRESULT ' + json.dumps(out))\n"
)
t0 = time.time()
try:
    cp = subprocess.run([sys.executable, "-c", KL_SNIPPET], capture_output=True, text=True, timeout=600)
    line = [ln for ln in cp.stdout.splitlines() if ln.startswith("KLRESULT ")]
    kl_res = json.loads(line[-1][len("KLRESULT "):]) if line else {
        "import_ok": False, "error": (cp.stderr.strip().splitlines() or ["无输出"])[-1]}
except subprocess.TimeoutExpired:
    kl_res = {"import_ok": False, "error": "TimeoutExpired (600 s)"}
kl_res["elapsed_s"] = round(time.time() - t0, 1)

erk_map = kmap[kmap["our_name"].isin(["MAPK1", "MAPK3"])][["our_name", "atlas_name", "kin_type", "status"]]

# ---- 基线比例
n_row_early = int(rows["peak_tp_max"].isin(EARLY).sum())
n_site_early = int((sites["peak_early_any"] == "True").sum())

# ================================================================ 写 tsv
rows[ROW_COLS].to_csv(OUT_ROWS, sep="\t", index=False)
sites[SITE_COLS].to_csv(OUT_SITES, sep="\t", index=False)
rows_back = read_tsv(OUT_ROWS)
sites_back = read_tsv(OUT_SITES)


def has_nan_text(df):
    return int(df.isin(["nan", "NaN", "None", "NA"]).sum().sum())


# ================================================================ 验收自检
check("1", "R 行数", nR_rows, 5649)
check("1", "R 位点数", nR_sites, 3787)
check("1", "R 蛋白数", nR_prot, 1837)
check("1", "A1 合计", A1[5], 5649)
check("1", "A1 中 90min缺失 列合计", A1[4]["90min缺失"], 801)
check("1", "R 中 log2occ_90min 为空行数", n_90_empty, 801)
check("1", "A1 合计 = A2 合计", A1[5] == A2[5], True)
check("1", "A1 与 A2 的 type_90 列合计逐列相同", A1[4] == A2[4], True)
check("1", "A3 合计", A3[5], 5649)
check("2", "位点类型四类之和", sum(site_type_cnt[t] for t in SITE_TYPES), 3787)
check("2", "n_frac_resp 1/2/3+ 之和", nfr[1] + nfr[2] + sum(nfr[k] for k in range(3, 7)), 3787)
check("2", "n_comp_resp 1/2/3 之和", ncr[1] + ncr[2] + ncr[3], 3787)
check("2", "n_frac 交叉表行边际 = n_frac_resp 分布", [XF[3][k] for k in range(1, 7)], [nfr[k] for k in range(1, 7)])
check("2", "n_frac 交叉表列边际 = n_frac_evaluable 分布",
      [XF[4][k] for k in range(1, 7)], [Counter(sites["n_frac_evaluable"])[k] for k in range(1, 7)])
check("2", "n_comp 交叉表行边际 = n_comp_resp 分布", [XC[3][k] for k in range(1, 4)], [ncr[k] for k in range(1, 4)])
check("2", "n_comp 交叉表列边际 = n_comp_evaluable 分布",
      [XC[4][k] for k in range(1, 4)], [Counter(sites["n_comp_evaluable"])[k] for k in range(1, 4)])
check("2", "交叉表合计（frac / comp）", (XF[5], XC[5]), (3787, 3787))
check("3", "peak_tp_max=90min 且 transient 的行数", A1[2].get(("90min", "transient"), 0), 0)
for col in ("single_backbone", "win5_backbone"):
    b = bb_res[col]
    check("4", f"{col}: d>0 + d<0 + d=0 = 蛋白数", b["pos"] + b["neg"] + b["zero"], b["n_prot"])
    check("4", f"{col}: 缺 backbone 位点数（方案验收 4：5；初稿写 0）", len(b["miss_R"]), 5)
check("5", "3c 写明“缺序列，跳过”且有 pip/import 记录", (not has_seq) and bool(pip_attempts) and "import_ok" in kl_res, True)
check("6", "rows.tsv 数据行数", len(rows_back), 5649)
check("6", "sites.tsv 数据行数", len(sites_back), 3787)
check("6", "rows.tsv / sites.tsv 中 nan/NaN/None/NA 文本个数", (has_nan_text(rows_back), has_nan_text(sites_back)), (0, 0))
check("6", "rows.tsv / sites.tsv 列名与方案一致", (list(rows_back.columns) == ROW_COLS, list(sites_back.columns) == SITE_COLS), (True, True))

# ================================================================ 写 md
design = open(P_DESIGN, encoding="utf-8").read().splitlines()
i0 = next(i for i, ln in enumerate(design) if ln.startswith("### 基集"))
i1 = next(i for i, ln in enumerate(design) if i > i0 and ln.startswith("### 输出"))
design_def = "\n".join(design[i0:i1]).strip()

L = []
A = L.append
A("# 03 时间规律小实验（任务 3）\n")
A("## ① 输入\n")
A(f"- `{rel(P_TRAJ)}`：{base['traj_rows']} 行（不含表头），{base['traj_sites']} 个 site_id，{base['traj_prot']} 个 protein。"
  f"class：" + "，".join(f"{k} {v}" for k, v in sorted(base['class'].items(), key=lambda x: -x[1]))
  + "；denom_tier：" + "，".join(f"{k} {v}" for k, v in sorted(base['tier'].items(), key=lambda x: -x[1])) + "。")
A(f"- `{rel(P_FEAT)}`：{len(feat)} 行，{feat['site_id'].nunique()} 个 site_id，{len(feat.columns)} 列；"
  f"site_traj 位点在其中的覆盖：{len(set(traj['site_id']) & set(feat['site_id']))}/{base['traj_sites']}。"
  f"`single_backbone` 为字面 \"NA\" 的行 {int((feat['single_backbone'] == 'NA').sum())}，"
  f"`win5_backbone` 为 \"NA\" 的行 {int((feat['win5_backbone'] == 'NA').sum())}（全表）。")
A(f"- `{rel(P_PROF)}`：{len(prof)} 行（仅在 3c 列出列名）。")
A(f"- `{rel(P_MAP)}`：{len(kmap)} 行（逗号分隔，`pd.read_csv(..., sep=',', dtype=str, keep_default_na=False)`）。")
A(f"- kinase-library wheel：`{WHEEL}`（{'存在' if os.path.exists(WHEEL) else '不存在'}）。")
A("- 读表：`pd.read_csv(path, sep='\\t', dtype=str, keep_default_na=False)`，空串 = 缺失，数值处显式 `float()`。")
A(f"- 基集 R（class ∈ {{both, occupancy_only}}）：{nR_rows} 行、{nR_sites} 位点、{nR_prot} 蛋白；"
  f"其中 `log2occ_90min` 为空 {n_90_empty} 行，`log2occ_CTRL` 为空 {n_ctrl_empty} 行。")
A(f"- 阈值（脚本已有，未改）：FC_THRESH = {FC_THRESH}，MIN_REPS = {MIN_REPS}，PCT = {PCT}（任务 3 无用到 PCT 之处）。")
A(f"- 脚本：`{rel(__file__)}`；输出：`{rel(OUT_ROWS)}`、`{rel(OUT_SITES)}`、`{rel(OUT_MD)}`。\n")

A("## ② 定义（照抄 00_design.md 任务 3 的 基集/3a/3b/3c，原文如下）\n")
A(design_def)
A("")
A("**执行补充（协调方更正，现已并入方案）**：site_features.tsv 的特征列缺失值为字面字符串 `\"NA\"`（非空串）。"
  "`\"NA\"` 视为缺失：该位点不参与 backbone 中位数计算，sites.tsv 中该位点的 `single_backbone`/`win5_backbone` 留空串；"
  "“缺 backbone 的位点数”预期为 5（方案初稿写 0）。")
A("")
A("**实现细节**：`peak_tp_max` 按 EGF_TPS 顺序（2min→8min→20min→90min）遍历，仅当 |fc| 严格大于当前最大值时替换，即并列取最早。"
  "rows.tsv 中 fc 列以 3 位小数写出（输入 `log2occ_*` 均 ≤3 位小数）；缺该时间点则为空串。`peak_tp_07` = site_traj 的 `peak_tp` 原值。"
  "`peak_early_any` 写作 `True`/`False`。区室：FR1/FR2→Cyt，FR3/FR4→Mem，FR5/FR6→Nuc。"
  "backbone 对照中，蛋白的入选（≥1 transient 且 ≥1 sustained）按剔除 backbone 缺失位点后的位点判断；中位数用 `statistics.median`。\n")

A("## ③ 3a 形状\n")
A("### 表 A1：peak_tp_max × type_90（行数）\n")
A(md_table(A1[0], A1[1]))
A("")
A("### 表 A2：07 peak_tp × type_90（行数）\n")
A(md_table(A2[0], A2[1]))
A("")
A("### 表 A3：peak_tp_max × 07 peak_tp（行数）\n")
A(md_table(A3[0], A3[1]))
A("")
A(f"- 对角线（peak_tp_max = peak_tp_07）行数：{sum(A3[2].get((t, t), 0) for t in EGF_TPS)}；非对角线：{A3[5] - sum(A3[2].get((t, t), 0) for t in EGF_TPS)}。")
A(f"- 每行 max|fc| < FC_THRESH 的响应行数：{n_rows_max_lt1}。")
A(f"- 精度核对（同一定义改用十进制精确减法重算）：peak_tp_max 不一致 {prec_peak_diff} 行，type_90 不一致 {prec_type_diff} 行；"
  f"十进制下 |fc_90| 恰为 1.000 的行 {n_fc90_eq1_exact}（按 `<1.0` 规则归 sustained）；"
  f"max|fc| 出现并列的行：浮点 {n_tie_float}，十进制 {n_tie_dec}（两种算法下 peak_tp_max 结果均一致）。\n")
A("### 位点级类型（忽略 90min缺失 行）\n")
A(md_table(["site_type_90", "位点数"], [[t, site_type_cnt[t]] for t in SITE_TYPES] + [["合计", sum(site_type_cnt.values())]]))
A("")
for col in ("single_backbone", "win5_backbone"):
    b = bb_res[col]
    A(f"### backbone 对照：`{col}`\n")
    A(md_table(["项", "值"], [
        ["同时有 ≥1 transient 与 ≥1 sustained 位点（均有 backbone 值）的蛋白数", b["n_prot"]],
        ["参与计算的 transient 位点数 / sustained 位点数", f"{b['n_tr_sites']} / {b['n_su_sites']}"],
        ["d 的中位数", f"{b['d_median']:.6g}" if b["d_median"] is not None else ""],
        ["d > 0 蛋白数", b["pos"]],
        ["d < 0 蛋白数", b["neg"]],
        ["d = 0 蛋白数（精确相等）", b["zero"]],
        ["abs(d) < 1e-9 蛋白数（容差核对）", b["zero_tol"]],
        ["缺 backbone 的位点数（R 的 3787 位点中；字面 \"NA\"）", f"{len(b['miss_R'])}（\"NA\" {b['raw_na_R']}）"],
        ["缺 backbone 位点的 site_type_90 分布", "，".join(f"{t} {b['miss_types'][t]}" for t in SITE_TYPES)],
        ["缺 backbone 位点 site_id", "，".join(b["miss_R"])],
        ["若不剔除缺失位点、按位点类型判断入选的蛋白数", b["n_prot_if_no_drop"]],
        ["因剔除缺失位点而落选的蛋白", "，".join(b["lost_prot"]) if b["lost_prot"] else "无"],
    ]))
    A("")

A("## ④ 3b 空间特异性\n")
A(f"- 响应行中不满足 evaluable 条件的行数：{n_resp_not_eval}。\n")
A("### 表 B1：n_frac_resp 分布（位点数）\n")
A(md_table(["n_frac_resp", "位点数"],
           [[1, nfr[1]], [2, nfr[2]], ["3+", sum(nfr[k] for k in range(3, 7))]]
           + [[f"  其中 {k}", nfr[k]] for k in range(3, 7)]
           + [["合计", nfr[1] + nfr[2] + sum(nfr[k] for k in range(3, 7))]]))
A("")
A("### 表 B2：n_frac_resp × n_frac_evaluable（位点数）\n")
A(md_table(XF[0], XF[1]))
A("")
A("### 表 B3：n_comp_resp 分布（位点数；FR1-2 Cyt，FR3-4 Mem，FR5-6 Nuc）\n")
A(md_table(["n_comp_resp", "位点数"], [[k, ncr[k]] for k in (1, 2, 3)] + [["合计", sum(ncr.values())]]))
A("")
A("### 表 B4：n_comp_resp × n_comp_evaluable（位点数）\n")
A(md_table(XC[0], XC[1]))
A("")

A("## ⑤ 3c ERK 阳性对照\n")
A("### 序列来源检查\n")
for k, v in col_lists.items():
    A(f"- `{k}` 列名（{len(v)}）：" + ", ".join(f"`{c}`" for c in v))
A("- 列名匹配 `seq|window|flank|motif|peptide|fasta|aa_`（不分大小写）：" +
  "；".join(f"`{k}` {len(v)} 个" + (f"（{', '.join(v)}）" if v else "") for k, v in seq_cols.items()))
A(f"- 仓库内（排除 .git）`find -name '*.fasta'`：{len(fasta_hits)} 个" + (f"（{'; '.join(fasta_hits)}）" if fasta_hits else ""))
A(f"- 仓库内（排除 .git）名称含 `biophys_json` 的文件/目录：{len(biophys_hits)} 个" + (f"（{'; '.join(biophys_hits)}）" if biophys_hits else ""))
A(f"- 结论：{'有序列来源（需复核）' if has_seq else '无序列 / 窗口来源'}。site_id（如 `P08195_S33`）只含残基与位置，不含侧翼序列。")
A(f"- kinase-library 打分：**{SKIP_SEQ}**。")
A(f"- “+1 位是 P” 代理：**{SKIP_SEQ}**。\n")
A("### pip / import / kinase 列表检查\n")
for a in pip_attempts:
    if a.get("attempt", 0) == 0:
        A(f"- `pip install`：{a['status']}")
    else:
        A(f"- `{os.path.basename(sys.executable)} -m pip install <wheel>` 第 {a['attempt']} 次：{a['status']}，"
          f"returncode={a['returncode']}，耗时 {a['elapsed_s']} s，输出末行：" + " / ".join(f"`{t[:200]}`" for t in a["tail"]))
A(f"- `import kinase_library`：{'成功' if kl_res.get('import_ok') else '失败'}"
  + (f"（version {kl_res.get('version')}）" if kl_res.get("version") else "")
  + (f"；错误：`{kl_res.get('error')}`" if kl_res.get("error") else "") + f"；耗时 {kl_res['elapsed_s']} s。")
if "n_ser_thr" in kl_res:
    A(f"- `kl.get_kinase_list(kin_type='ser_thr')`：{kl_res['n_ser_thr']} 个激酶；含 `ERK1`：{kl_res['has_ERK1']}；含 `ERK2`：{kl_res['has_ERK2']}。")
else:
    A("- `kl.get_kinase_list(kin_type='ser_thr')`：未得到结果（见上方 import 结果）。")
A(f"- 调用方式参考：`{KL_CORE_REF}`。")
A("- pip 前后环境版本（`importlib.metadata`）：" + "；".join(
    f"{k} {ver_before[k]} → {ver_after[k]}" for k in VER_PKGS)
  + f"。本脚本进程内的 pandas 为 pip 前已加载的 {ver_inproc['pandas']}。")
A("- kinase_name_mapping.csv 中 ERK 映射：" + "；".join(
    f"{r.our_name}→{r.atlas_name}（{r.kin_type}，{r.status}）" for r in erk_map.itertuples()) + "\n")
A("### 基线比例\n")
A(md_table(["口径", "peak_tp_max ∈ {2min, 8min}", "分母", "比例", "ERK 候选比例"], [
    ["行（全部响应行）", n_row_early, nR_rows, f"{n_row_early / nR_rows:.4f}", SKIP_SEQ],
    ["位点（任一响应行）", n_site_early, nR_sites, f"{n_site_early / nR_sites:.4f}", SKIP_SEQ],
]))
A("")

A("## ⑥ 跳过项\n")
A(f"- 3c kinase-library 打分（ERK1/ERK2）：缺序列，跳过。")
A(f"- 3c “+1 位是 P” 代理：缺序列，跳过。")
A(f"- 3c ERK 候选比例（行口径、位点口径）：缺序列，跳过。")
A("- 未做任何统计检验（按方案）。")
A(f"- backbone 对照中缺 backbone（字面 \"NA\"）的位点：single {len(bb_res['single_backbone']['miss_R'])} 个、"
  f"win5 {len(bb_res['win5_backbone']['miss_R'])} 个，不参与中位数计算。\n")

A("### 附：验收自检（脚本自动核对）\n")
A(md_table(["#", "项", "实际", "期望", "通过"],
           [[n, d, a if not isinstance(a, list) or len(str(a)) < 120 else "（列表，见脚本）",
             e if not isinstance(e, list) or len(str(e)) < 120 else "（列表，见脚本）", "是" if ok else "否"]
            for n, d, a, e, ok in checks]))
A("")

with open(OUT_MD, "w", encoding="utf-8") as f:
    f.write("\n".join(L))

# ================================================================ 控制台
print(f"R: {nR_rows} rows / {nR_sites} sites / {nR_prot} proteins; 90min empty {n_90_empty}")
print("A1:", A1[1])
print("site types:", dict(site_type_cnt))
for col in ("single_backbone", "win5_backbone"):
    b = bb_res[col]
    print(col, {k: b[k] for k in ("n_prot", "d_median", "pos", "neg", "zero", "zero_tol", "n_prot_if_no_drop", "lost_prot")},
          "missing", len(b["miss_R"]))
print("nfr", dict(nfr), "ncr", dict(ncr))
print("pip", pip_attempts)
print("kl", kl_res)
print("early rows", n_row_early, "sites", n_site_early)
for n, d, a, e, ok in checks:
    print("PASS" if ok else "FAIL", n, d, a if len(str(a)) < 100 else "...", "expected", e if len(str(e)) < 100 else "...")
