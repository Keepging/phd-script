#!/usr/bin/env python3
"""07_pilot2 任务 2：漏斗重跑 + 蛋白丰度。

按 analysis/07_pilot2/00_design.md（通用约定 + "任务 2"）与
analysis/06_pilot/00_design.md（"任务 1：漏斗表" 的 c2/c3/c4/c5/c6/S3/响应位点定义与按 fraction 报法）执行。
输入（只读）：
  analysis/07_pilot2/01_site_traj_v2.tsv      新分类（denom_tier 全为 protein）
  data/pilot/occupancy_variants_long.tsv      蛋白丰度（CTRL 的 denom_protein）
  data/pilot/site_traj.tsv                    site_id→protein 映射；旧 log2int 对照
  analysis/06_pilot/01_funnel.tsv             旧漏斗数（直接抄进并列表）
输出：
  analysis/07_pilot2/02_protein_abundance.tsv
  analysis/07_pilot2/02_funnel_v2.tsv
  analysis/07_pilot2/02_funnel_v2.md
复跑：python3 analysis/07_pilot2/02_funnel_v2.py   （任意目录）
"""
import math
import os
import statistics

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
P_V2 = os.path.join(HERE, "01_site_traj_v2.tsv")
P_LONG = os.path.join(ROOT, "data", "pilot", "occupancy_variants_long.tsv")
P_ST = os.path.join(ROOT, "data", "pilot", "site_traj.tsv")
P_OLD = os.path.join(ROOT, "analysis", "06_pilot", "01_funnel.tsv")
REL = {
    P_V2: "analysis/07_pilot2/01_site_traj_v2.tsv",
    P_LONG: "data/pilot/occupancy_variants_long.tsv",
    P_ST: "data/pilot/site_traj.tsv",
    P_OLD: "analysis/06_pilot/01_funnel.tsv",
}
OUT_ABD = os.path.join(HERE, "02_protein_abundance.tsv")
OUT_TSV = os.path.join(HERE, "02_funnel_v2.tsv")
OUT_MD = os.path.join(HERE, "02_funnel_v2.md")

# 脚本已有常量，不改
FC_THRESH = 1.0
MIN_REPS = 2
EGF_TPS = ["2min", "8min", "20min", "90min"]
TPS = ["CTRL"] + EGF_TPS
FRS = [f"FR{i}" for i in range(1, 7)]
OCC_COLS = [f"log2occ_{tp}" for tp in TPS]
INT_COLS = [f"log2int_{tp}" for tp in TPS]
NUM_COLS = OCC_COLS + INT_COLS
RESP_CLASSES = {"both", "occupancy_only"}
NONE_MARK = "无"

COLUMNS = ["version", "step", "condition", "scope", "n_rows", "n_sites", "n_proteins",
           "n_proteins_with_abundance", "protein_abundance_median", "site_log2int_CTRL_median"]


def read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


missing = [REL[p] for p in (P_V2, P_LONG, P_ST, P_OLD) if not os.path.exists(p)]
if missing:
    raise SystemExit("缺 " + "、".join(missing) + "，跳过")

v2 = read(P_V2)
lg = read(P_LONG)
st = read(P_ST)
old = read(P_OLD)
for name, d in (("v2", v2), ("long", lg), ("site_traj", st), ("old_funnel", old)):
    assert d.isna().sum().sum() == 0, name

# ---------------------------------------------------------------- 结构核对
assert list(v2.columns) == list(st.columns)
assert len(v2) == len(st) == 39813
assert not v2.duplicated(["site_id", "fraction"]).any()
assert (v2[["site_id", "fraction"]].values == st[["site_id", "fraction"]].values).all()
assert (st.groupby("site_id")["protein"].nunique() == 1).all()
site2prot = st.drop_duplicates("site_id").set_index("site_id")["protein"]
# 蛋白归属一律用 site_traj 的映射；核对 v2 自带 protein 列与之相同
v2_prot_mapped = v2["site_id"].map(site2prot)
n_v2_prot_diff = int((v2_prot_mapped != v2["protein"]).sum())
v2["protein"] = v2_prot_mapped
tier_counts = v2["denom_tier"].value_counts().to_dict()
class_counts = v2["class"].value_counts().to_dict()

# log2int 新旧核对（分母不影响强度：有无必须一致；数值记录差）
int_cmp = []
for c in INT_COLS:
    a, b = v2[c] != "", st[c] != ""
    both = a & b
    diff = (v2.loc[both, c].astype(float) - st.loc[both, c].astype(float)).abs()
    int_cmp.append({"列": c, "v2 非空": int(a.sum()), "site_traj 非空": int(b.sum()),
                    "有无不一致格": int((a != b).sum()),
                    "数值字符串不同格": int((v2.loc[both, c] != st.loc[both, c]).sum()),
                    "最大绝对差": f"{diff.max():.3f}" if len(diff) else NONE_MARK})

# ---------------------------------------------------------------- 蛋白丰度
n_long = len(lg)
assert lg["site_id"].isin(site2prot.index).all(), "长表有 site_id 不在 site_traj"
lg["protein"] = lg["site_id"].map(site2prot)
ctrl = lg[(lg["timepoint"] == "CTRL") & (lg["denom_protein"] != "")]
n_ctrl_rows_all = int((lg["timepoint"] == "CTRL").sum())
n_ctrl_rows_denom = len(ctrl)
# 同一 protein×run 的 denom_protein 必须唯一（advisor 已核全部 timepoint；此处再核）
nn = lg[lg["denom_protein"] != ""]
g_all = nn.groupby(["protein", "timepoint", "fraction", "rep"])["denom_protein"].nunique()
n_prun_groups = len(g_all)
n_prun_multi = int((g_all > 1).sum())
g_ctrl = ctrl.groupby(["protein", "fraction", "rep"])["denom_protein"].nunique()
assert (g_ctrl == 1).all(), "CTRL 内同一 (protein, fraction, rep) 的 denom_protein 不唯一"
ded = ctrl.drop_duplicates(["protein", "fraction", "rep"]).copy()
ded["log2_denom"] = [math.log2(float(x)) for x in ded["denom_protein"]]
assert all(float(x) > 0 for x in ded["denom_protein"])
abd_rows = []
for prot, sub in ded.groupby("protein", sort=True):
    abd_rows.append({"protein": prot,
                     "n_ctrl_runs_with_denom": int(len(sub)),
                     "abundance_full": statistics.median(sub["log2_denom"].tolist())})
abd = pd.DataFrame(abd_rows)
assert abd["n_ctrl_runs_with_denom"].max() <= 24
abd_out = abd[["protein", "n_ctrl_runs_with_denom"]].copy()
abd_out["abundance_median_log2"] = [f"{x:.4f}" for x in abd["abundance_full"]]
abd_out.to_csv(OUT_ABD, sep="\t", index=False, lineterminator="\n")
prot2abd = dict(zip(abd["protein"], abd["abundance_full"]))
prot2abd_r4 = dict(zip(abd_out["protein"], abd_out["abundance_median_log2"].astype(float)))

# ---------------------------------------------------------------- 行级条件（v2）
df = v2


def nonempty(col):
    return df[col] != ""


c2_tp = pd.Series(False, index=df.index)
c3_tp = pd.Series(False, index=df.index)
for tp in EGF_TPS:
    c2_tp |= nonempty(f"log2int_{tp}")
    c3_tp |= nonempty(f"log2occ_{tp}")
c2 = nonempty("log2int_CTRL") & c2_tp
c3 = nonempty("log2occ_CTRL") & c3_tp
c23 = c2 & c3
any_occ = pd.Series(False, index=df.index)
for col in OCC_COLS:
    any_occ |= nonempty(col)
all_num_empty = pd.Series(True, index=df.index)
for col in NUM_COLS:
    all_num_empty &= ~nonempty(col)

# ---------------------------------------------------------------- 位点级 / 蛋白级条件
S3 = set(df.loc[c23, "site_id"])
n_fr_c23 = df.loc[c23].groupby("site_id")["fraction"].nunique()
c4_sites = set(n_fr_c23[n_fr_c23 >= 2].index)

prot_n_s3 = site2prot.loc[sorted(S3)].value_counts()
c5_prots = set(prot_n_s3[prot_n_s3 >= 2].index)
c5_sites = {s for s in S3 if site2prot[s] in c5_prots}

resp_sites_all = set(df.loc[df["class"].isin(RESP_CLASSES), "site_id"])
resp_in_s3 = resp_sites_all & S3
nonresp_sites = S3 - resp_sites_all
prots_with_resp = {site2prot[s] for s in resp_in_s3}
prots_with_nonresp = {site2prot[s] for s in nonresp_sites}
c6_prots = prots_with_resp & prots_with_nonresp
c6_sites = {s for s in S3 if site2prot[s] in c6_prots}

sites_1 = set(df["site_id"])
sites_2 = set(df.loc[c2, "site_id"])
sites_3 = S3
A4 = sites_3 & c4_sites
A5 = A4 & c5_sites
A6 = A5 & c6_sites
B5 = sites_3 & c5_sites
B6 = B5 & c6_sites

all_rows = pd.Series(True, index=df.index)
STEPS = [
    ("A", "1", "全部", all_rows, sites_1),
    ("A", "2", "行满足 c2", c2, sites_2),
    ("A", "3", "行满足 c2&c3；位点 ∈ S3", c23, sites_3),
    ("A", "4", "行满足 c2&c3；位点 ∈ S3 ∩ c4", c23, A4),
    ("A", "5", "行满足 c2&c3；位点 ∈ ④剩余 ∩ c5", c23, A5),
    ("A", "6", "行满足 c2&c3；位点 ∈ ⑤剩余 ∩ c6", c23, A6),
    ("A", "ctrl_c", "宽口径：任一 log2occ_* 非空", any_occ, None),
    ("A", "ctrl_d", "全表中 10 个数值列（log2occ_*、log2int_*）全空", all_num_empty, None),
    ("B", "1", "全部", all_rows, sites_1),
    ("B", "2", "行满足 c2", c2, sites_2),
    ("B", "3", "行满足 c2&c3；位点 ∈ S3", c23, sites_3),
    ("B", "5", "行满足 c2&c3；位点 ∈ ③剩余 ∩ c5", c23, B5),
    ("B", "6", "行满足 c2&c3；位点 ∈ ⑤剩余 ∩ c6", c23, B6),
]


def fmt_median(vals):
    if len(vals) == 0:
        return NONE_MARK
    return f"{statistics.median(vals):.4f}"


records = []
r4_mismatch = []  # 用 4 位小数丰度值取中位数时与全精度结果不同的格
for ver, step, cond, rowmask, siteset in STEPS:
    mask = rowmask.copy()
    if siteset is not None:
        mask &= df["site_id"].isin(siteset)
    for scope in ["all"] + FRS:
        m = mask if scope == "all" else mask & (df["fraction"] == scope)
        sub = df.loc[m]
        prots = sorted(set(sub["protein"]))
        prots_ab = [p for p in prots if p in prot2abd]
        ab_full = fmt_median([prot2abd[p] for p in prots_ab])
        ab_r4 = fmt_median([prot2abd_r4[p] for p in prots_ab])
        if ab_full != ab_r4:
            r4_mismatch.append((ver, step, scope, ab_full, ab_r4))
        ctrl_vals = [float(x) for x in sub["log2int_CTRL"] if x != ""]
        records.append({
            "version": ver, "step": step, "condition": cond, "scope": scope,
            "n_rows": int(len(sub)),
            "n_sites": int(sub["site_id"].nunique()),
            "n_proteins": len(prots),
            "n_proteins_with_abundance": len(prots_ab),
            "protein_abundance_median": ab_full,
            "site_log2int_CTRL_median": fmt_median(ctrl_vals),
        })
res = pd.DataFrame(records, columns=COLUMNS)
res.to_csv(OUT_TSV, sep="\t", index=False, lineterminator="\n")

# ---------------------------------------------------------------- 自检
def get(ver, step, scope, col):
    r = res[(res["version"] == ver) & (res["step"] == step) & (res["scope"] == scope)]
    assert len(r) == 1, (ver, step, scope)
    return r[col].iloc[0]


checks = []


def chk(name, ok, detail=""):
    checks.append((name, bool(ok), detail))


s1 = tuple(int(get("A", "1", "all", c)) for c in ("n_rows", "n_sites", "n_proteins"))
chk("A step1/all = 39813/18268/4996", s1 == (39813, 18268, 4996), "/".join(map(str, s1)))
v = int(get("A", "2", "all", "n_rows"))
chk("A step2/all n_rows = 14963（强度不随分母变）", v == 14963, f"实得 {v}")
old_s2 = int(old[(old["step"] == "2") & (old["scope"] == "all")]["n_rows"].iloc[0])
chk("A step2/all n_rows = 06_pilot 旧漏斗 step2/all n_rows", v == old_s2, f"{v} vs {old_s2}")
bad = []
for step in ("1", "2", "3"):
    for scope in ["all"] + FRS:
        for col in COLUMNS[4:]:
            if get("A", step, scope, col) != get("B", step, scope, col):
                bad.append(f"{step}/{scope}/{col}")
chk("B step1–3 与 A 相同（7 个 scope × 6 个数值列逐格）", not bad, "; ".join(bad) or "全部相同")
for ver, steps in (("A", ["1", "2", "3", "4", "5", "6"]), ("B", ["1", "2", "3", "5", "6"])):
    for col in ("n_rows", "n_sites", "n_proteins", "n_proteins_with_abundance"):
        seq = [int(get(ver, s, "all", col)) for s in steps]
        chk(f"{ver} all 口径 {col} 单调不增", all(x >= y for x, y in zip(seq, seq[1:])),
            "→".join(map(str, seq)))
a5, b5 = int(get("A", "5", "all", "n_sites")), int(get("B", "5", "all", "n_sites"))
chk("B step5 n_sites ≥ A step5 n_sites", b5 >= a5, f"B {b5} ≥ A {a5}")
bad = []
for ver, step, *_ in STEPS:
    tot = sum(int(get(ver, step, fr, "n_rows")) for fr in FRS)
    if tot != int(get(ver, step, "all", "n_rows")):
        bad.append(f"{ver}{step}: {tot}≠{get(ver, step, 'all', 'n_rows')}")
chk("每步 6 个 fraction n_rows 之和 = all（13 个 version×step）", not bad, "; ".join(bad) or "全部相等")
chk("n_proteins_with_abundance ≤ n_proteins（全部 91 行）",
    (res["n_proteins_with_abundance"] <= res["n_proteins"]).all())
chk("丰度表蛋白数 = 3356（advisor 参考）且 ≤ 4996", len(abd_out) == 3356 and len(abd_out) <= 4996,
    f"{len(abd_out)}")
chk("丰度表蛋白 ⊆ site_traj 蛋白", set(abd_out["protein"]) <= set(site2prot.values))
for prot, n_ref, v_ref in (("Q8NB90", 1, "28.9955"), ("Q6PGN9", 3, "24.2320"), ("Q99733", 14, "27.1755")):
    r = abd_out[abd_out["protein"] == prot]
    got = (int(r["n_ctrl_runs_with_denom"].iloc[0]), r["abundance_median_log2"].iloc[0]) if len(r) else None
    chk(f"丰度抽查 {prot} n={n_ref} 值 {v_ref}", got == (n_ref, v_ref), f"实得 {got}")
chk("CTRL 中同一 (protein, fraction, rep) 的 denom_protein 唯一", (g_ctrl == 1).all(),
    f"{len(g_ctrl)} 组")
chk("全部 timepoint 同一 protein×run 的 denom_protein 唯一（advisor：125260 组）",
    n_prun_multi == 0 and n_prun_groups == 125260, f"{n_prun_groups} 组，不唯一 {n_prun_multi}")
chk("v2 denom_tier 全为 protein", tier_counts == {"protein": 39813}, str(tier_counts))
chk("v2 protein 列 = site_traj 映射", n_v2_prot_diff == 0, f"不同 {n_v2_prot_diff}")
chk("log2int_* 有无 v2 与 site_traj 逐格一致", all(r["有无不一致格"] == 0 for r in int_cmp))
nan_free = not res.isna().any().any() and not res.astype(str).apply(
    lambda s: s.str.lower().isin(["nan", "none", ""])).any().any()
chk("tsv 无 NaN/空格", nan_free)
grp = res.groupby(["version", "step", "scope"]).size()
chk("tsv 91 行（A 8 step × 7 + B 5 step × 7），每 version×step×scope 恰 1 行",
    len(res) == 91 and (grp == 1).all(), f"{len(res)} 行")

aux = [
    ("|S3|（位点）", len(S3)),
    ("S3 中满足 c4 的位点", len(c4_sites)),
    ("满足 c5 的蛋白（S3 中 ≥2 位点）", len(c5_prots)),
    ("S3 中满足 c5 的位点（= 版本 B ⑤ 位点）", len(c5_sites)),
    ("occupancy 响应位点（全表，任一行 class∈{both,occupancy_only}）", len(resp_sites_all)),
    ("其中 ∈ S3", len(resp_in_s3)),
    ("其中 ∉ S3", len(resp_sites_all - S3)),
    ("非响应位点（∈ S3 且无响应行）", len(nonresp_sites)),
    ("满足 c6 的蛋白", len(c6_prots)),
    ("S3 中满足 c6 的位点", len(c6_sites)),
    ("版本 A ⑤ 位点（S3 ∩ c4 ∩ c5）", len(A5)),
    ("版本 A ⑥ 位点", len(A6)),
    ("版本 B ⑤ 位点（S3 ∩ c5）", len(B5)),
    ("版本 B ⑥ 位点", len(B6)),
    ("class∈{both,occupancy_only} 但不满足 c2&c3 的行", int((df["class"].isin(RESP_CLASSES) & ~c23).sum())),
]

# ---------------------------------------------------------------- md
def md_table(frame):
    cols = list(frame.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in frame.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


DEFS = """沿用 `analysis/07_pilot2/00_design.md` 通用约定与 `analysis/06_pilot/00_design.md` 任务 1 定义，全部基于 v2 表（`01_site_traj_v2.tsv`）。

### 蛋白丰度（07_pilot2 任务 2）
`protein_abundance = median over CTRL runs of log2(denom_protein)`：对每个蛋白，取长表中 `timepoint == CTRL` 且 `denom_protein` 非空的 (fraction, rep) 组合（最多 24 个 run），**按 (protein, fraction, rep) 去重**（同一蛋白多个位点的 `denom_protein` 是同一个数），取 `log2(denom_protein)` 的中位数；用原值，不做 run 中心化。protein 由 site_id 经 `site_traj.tsv` 的 protein 列映射。每层漏斗的 `protein_abundance_median` = 剩余行的位点所在蛋白（去重）中有丰度值者的 `abundance_median_log2` 的中位数（4 位小数）；无值的蛋白不参与，参与蛋白数报在 `n_proteins_with_abundance`；全无值则填 `无`。

### 行级条件（在一个 位点×fraction 行上判断）
- `c2`（②）：`log2int_CTRL != ''` 且 `log2int_{tp} != ''` 对某个 tp ∈ EGF_TPS 成立。
- `c3`（③）：`log2occ_CTRL != ''` 且 `log2occ_{tp} != ''` 对某个 tp ∈ EGF_TPS 成立（v2 里分母只有 protein 档）。
- 对照行：ctrl_c = 任一 `log2occ_*` 非空（不要求 CTRL+EGF 配对）；ctrl_d = 全表中 10 个数值列全空。旧的 ctrl_a/ctrl_b（按 denom_tier 拆 ③）不再适用（v2 全为 protein 档），不报。

### 位点级 / 蛋白级条件
- `S3` = 在 ≥1 个 fraction 上满足 `c2 & c3` 的位点集合。
- `c4`（④）：位点满足 `c2 & c3` 的 fraction 数 ≥ 2。
- `c5`（⑤）：位点所在蛋白在 `S3` 中有 ≥2 个位点。
- 位点"occupancy 响应" := 该位点任一行（v2）`class ∈ {both, occupancy_only}`；"非响应" := 位点 ∈ S3 且无此类行。
- `c6`（⑥）：位点所在蛋白在 `S3` 中既有 ≥1 个响应位点又有 ≥1 个非响应位点。

### 漏斗（累积：每步在上一步剩余集合上再加条件）
| version | step | 剩余行 | 剩余位点 |
|---|---|---|---|
| A, B | ① 全部 | 全表 | 全表 |
| A, B | ② | 行满足 c2 | 有 ≥1 行满足 c2 的位点 |
| A, B | ③ | 行满足 c2 & c3 | S3 |
| A | ④ | 行满足 c2&c3 且位点满足 c4 | S3 ∩ c4 |
| A | ⑤ | 行满足 c2&c3 且位点 ∈ (④剩余 ∩ c5) | S3 ∩ c4 ∩ c5 |
| A | ⑥ | 行满足 c2&c3 且位点 ∈ (⑤剩余 ∩ c6) | S3 ∩ c4 ∩ c5 ∩ c6 |
| B | ⑤ | 行满足 c2&c3 且位点 ∈ (③剩余 ∩ c5) | S3 ∩ c5 |
| B | ⑥ | 行满足 c2&c3 且位点 ∈ (⑤剩余 ∩ c6) | S3 ∩ c5 ∩ c6 |

c4/c5/c6 的位点集合均按上述定义在 S3 上算一次，两版共用；两版只差是否取 ∩ c4。

每步报：`n_rows`、`n_sites`（distinct site_id）、`n_proteins`（distinct protein）、`n_proteins_with_abundance`、`protein_abundance_median`、`site_log2int_CTRL_median` = 剩余行中非空 `log2int_CTRL`（v2 表的值）的中位数（位点级磷酸肽强度，非蛋白丰度）。
**按 fraction 再报一次**：scope = FR1..FR6，行 = 该 fraction 内满足该步行级条件且位点属于该步剩余集合的行；n_sites / n_proteins / 丰度在这些行上算（丰度是蛋白级、跨 fraction 的单一值，不按 fraction 重算）。

常量：`EGF_TPS = [2min, 8min, 20min, 90min]`，`FRS = FR1..FR6`，`FC_THRESH = 1.0`，`MIN_REPS = 2`。读表：`pd.read_csv(path, sep='\\t', dtype=str, keep_default_na=False)`，空串 = 缺失；中位数 `statistics.median`，4 位小数。"""

md = []
md.append("# 02 漏斗重跑（v2）+ 蛋白丰度\n")
md.append("## ① 输入\n")
md.append(f"- `{REL[P_V2]}`：{len(v2)} 行（不含表头），{v2['site_id'].nunique()} 个 site_id，"
          f"{v2['protein'].nunique()} 个 protein，{v2.shape[1]} 列；denom_tier 计数 "
          + "，".join(f"{k} {v}" for k, v in tier_counts.items()) + "；class 计数 "
          + "，".join(f"{k} {class_counts.get(k, 0)}" for k in ["none", "both", "occupancy_only", "intensity_only"]) + "。")
md.append(f"- `{REL[P_LONG]}`：{n_long} 行；timepoint==CTRL {n_ctrl_rows_all} 行，其中 denom_protein 非空 "
          f"{n_ctrl_rows_denom} 行；按 (protein, fraction, rep) 去重后 {len(ded)} 个 protein×CTRL run。")
md.append(f"- `{REL[P_ST]}`：{len(st)} 行；只用 site_id→protein 映射（{len(site2prot)} 个 site_id → "
          f"{site2prot.nunique()} 个 protein）和 log2int 对照。v2 自带 protein 列与该映射不同的行：{n_v2_prot_diff}。")
md.append(f"- `{REL[P_OLD]}`：{len(old)} 行（旧漏斗，直接抄进 ⑥）。")
md.append("- 脚本：`analysis/07_pilot2/02_funnel_v2.py`；输出：`02_funnel_v2.tsv`（91 行）、`02_funnel_v2.md`、"
          f"`02_protein_abundance.tsv`（{len(abd_out)} 行）。\n")
md.append("### log2int_* v2 vs site_traj 核对\n")
md.append(md_table(pd.DataFrame(int_cmp)) + "\n")

md.append("## ② 定义\n")
md.append(DEFS + "\n")
md.append("### 定义执行时的中间计数\n")
md.append(md_table(pd.DataFrame([{"量": k, "值": v} for k, v in aux])) + "\n")

md.append("## ③ 丰度表摘要（`02_protein_abundance.tsv`）\n")
vals = abd["abundance_full"].tolist()
md.append(f"- 有丰度值的蛋白：{len(abd_out)}（site_traj 共 {site2prot.nunique()} 个蛋白；无值 {site2prot.nunique() - len(abd_out)} 个，不列入丰度表）。")
md.append(f"- `abundance_median_log2`：最小 {min(vals):.4f}，中位数 {statistics.median(vals):.4f}，最大 {max(vals):.4f}。")
md.append(f"- 一致性：长表全部 timepoint 下 protein×run 组 {n_prun_groups} 个，denom_protein 取值不唯一的组 {n_prun_multi} 个；"
          f"CTRL 下 (protein, fraction, rep) 组 {len(g_ctrl)} 个，全部唯一。")
md.append("- advisor 抽查：" + "；".join(
    f"{p} n={int(abd_out.loc[abd_out.protein == p, 'n_ctrl_runs_with_denom'].iloc[0])} 值 "
    f"{abd_out.loc[abd_out.protein == p, 'abundance_median_log2'].iloc[0]}" for p in ("Q8NB90", "Q6PGN9", "Q99733")) + "。\n")
dist = abd_out["n_ctrl_runs_with_denom"].value_counts().reindex(range(1, 25), fill_value=0)
md.append("`n_ctrl_runs_with_denom` 分布：\n")
md.append(md_table(pd.DataFrame({"n_ctrl_runs_with_denom": dist.index, "n_proteins": dist.values,
                                 "累计": dist.values.cumsum()})) + "\n")
md.append(f"合计 {int(dist.sum())}。\n")

show_cols = [c for c in COLUMNS if c != "version"]
md.append("## ④ 版本 A（① → ② → ③ → ④ → ⑤ → ⑥ + ctrl_c、ctrl_d）\n")
md.append(md_table(res[res["version"] == "A"][show_cols]) + "\n")
md.append("## ⑤ 版本 B（① → ② → ③ → ⑤ → ⑥，⑤ 在 ③ 剩余上取交）\n")
md.append(md_table(res[res["version"] == "B"][show_cols]) + "\n")

md.append("## ⑥ 旧漏斗（06_pilot，all 口径）vs 本轮版本 A（all 口径）\n")
md.append("旧数直接抄自 `analysis/06_pilot/01_funnel.tsv`（site_traj，旧分类）；新 = 本轮版本 A（v2）。Δ = 新 − 旧。\n")
cmp_rows = []
old_all = old[old["scope"] == "all"].set_index("step")
new_all = res[(res["version"] == "A") & (res["scope"] == "all")].set_index("step")
for step in ["1", "2", "3", "4", "5", "6", "ctrl_a", "ctrl_b", "ctrl_c", "ctrl_d"]:
    o = old_all.loc[step]
    row = {"step": step}
    if step in new_all.index:
        n = new_all.loc[step]
        for c in ("n_rows", "n_sites", "n_proteins"):
            row[f"旧 {c}"] = int(o[c])
            row[f"新 {c}"] = int(n[c])
            row[f"Δ {c}"] = int(n[c]) - int(o[c])
        row["旧 protein_abundance_median"] = o["protein_abundance_median"]
        row["新 protein_abundance_median"] = n["protein_abundance_median"]
        row["旧 site_log2int_CTRL_median"] = o["site_log2int_CTRL_median"]
        row["新 site_log2int_CTRL_median"] = n["site_log2int_CTRL_median"]
    else:
        for c in ("n_rows", "n_sites", "n_proteins"):
            row[f"旧 {c}"] = int(o[c])
            row[f"新 {c}"] = "不适用"
            row[f"Δ {c}"] = "不适用"
        row["旧 protein_abundance_median"] = o["protein_abundance_median"]
        row["新 protein_abundance_median"] = "不适用"
        row["旧 site_log2int_CTRL_median"] = o["site_log2int_CTRL_median"]
        row["新 site_log2int_CTRL_median"] = "不适用"
    cmp_rows.append(row)
cmp_df = pd.DataFrame(cmp_rows)
md.append(md_table(cmp_df) + "\n")
md.append("旧 ctrl_a / ctrl_b（③ 按 denom_tier=protein / residue 拆分）在 v2 中不适用（v2 全为 protein 档）。\n")

md.append("## ⑦ 跳过/缺失项\n")
md.append(f"- 丰度：site_traj 中 {site2prot.nunique() - len(abd_out)} 个蛋白在 CTRL run 中无 `denom_protein`，无丰度值；"
          "不列入 `02_protein_abundance.tsv`，不参与 `protein_abundance_median`，体现在 `n_proteins_with_abundance` < `n_proteins`。")
emp = res[res["protein_abundance_median"] == NONE_MARK][["version", "step", "scope"]]
md.append("- `protein_abundance_median` 填 `无` 的格："
          + ("，".join(f"{r.version}{r.step}/{r.scope}" for r in emp.itertuples()) if len(emp) else "无") + "。")
emp = res[res["site_log2int_CTRL_median"] == NONE_MARK][["version", "step", "scope"]]
md.append("- `site_log2int_CTRL_median` 填 `无` 的格（剩余行中非空 `log2int_CTRL` 为 0 个）："
          + ("，".join(f"{r.version}{r.step}/{r.scope}" for r in emp.itertuples()) if len(emp) else "无") + "。")
n_int_diff = sum(r["数值字符串不同格"] for r in int_cmp)
md.append(f"- v2 的 `log2int_*` 与 site_traj 有无逐格一致，数值有 {n_int_diff} 格相差 ≤ 0.001"
          f"（其中 `log2int_CTRL` {int_cmp[0]['数值字符串不同格']} 格）；本表 `site_log2int_CTRL_median` 用 v2 值，"
          "故与旧漏斗该列可能在第 4 位小数上不同，行/位点/蛋白计数不受影响。")
md.append("- 旧漏斗的 ctrl_a、ctrl_b 不适用（见 ⑥）。")
md.append("- 口径说明（方案未指定）：`protein_abundance_median` 用全精度 `abundance_median_log2` 取中位数后保留 4 位；"
          "若改用 `02_protein_abundance.tsv` 中已取 4 位的值，91 格中有 "
          f"{len(r4_mismatch)} 格第 4 位小数相差 0.0001："
          + ("，".join(f"{a}{b}/{c} {d}→{e}" for a, b, c, d, e in r4_mismatch) if r4_mismatch else "无") + "。")
md.append("- 未使用：`data/pilot/run_medians.tsv`、`site_features.tsv`、`site_profiles.tsv`、`kinase_name_mapping.csv`（本任务不需要）。\n")

md.append("## 附：自检\n")
md.append(md_table(pd.DataFrame(
    [{"项": n, "结果": "通过" if ok else "未通过", "说明": d} for n, ok, d in checks])) + "\n")

with open(OUT_MD, "w", encoding="utf-8") as fh:
    fh.write("\n".join(md) + "\n")

# ---------------------------------------------------------------- 终端汇报
pd.set_option("display.width", 250)
print(res[res["scope"] == "all"].to_string(index=False))
print()
print(cmp_df.to_string(index=False))
print()
for k, v in aux:
    print(f"{k}\t{v}")
print()
for n, ok, d in checks:
    print(("PASS" if ok else "FAIL") + "\t" + n + "\t" + d)
