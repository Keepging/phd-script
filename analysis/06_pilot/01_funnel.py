#!/usr/bin/env python3
"""06_pilot 任务 1：漏斗表。

按 analysis/06_pilot/00_design.md 的 "0. 通用约定" 与 "任务 1" 执行。
输入 : data/pilot/site_traj.tsv（只读）
输出 : analysis/06_pilot/01_funnel.tsv, analysis/06_pilot/01_funnel.md
复跑 : python3 analysis/06_pilot/01_funnel.py   （在仓库根目录或任意目录均可）
"""
import os
import statistics

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
IN_PATH = os.path.join(ROOT, "data", "pilot", "site_traj.tsv")
IN_REL = "data/pilot/site_traj.tsv"
OUT_TSV = os.path.join(HERE, "01_funnel.tsv")
OUT_MD = os.path.join(HERE, "01_funnel.md")

# 脚本已有常量（07 :22-25），不改
EGF_TPS = ["2min", "8min", "20min", "90min"]
TPS = ["CTRL"] + EGF_TPS
FRS = [f"FR{i}" for i in range(1, 7)]
OCC_COLS = [f"log2occ_{tp}" for tp in TPS]
INT_COLS = [f"log2int_{tp}" for tp in TPS]
NUM_COLS = OCC_COLS + INT_COLS
RESP_CLASSES = {"both", "occupancy_only"}
NONE_MARK = "无"

COLUMNS = ["step", "condition", "scope", "n_rows", "n_sites", "n_proteins",
           "protein_abundance_median", "site_log2int_CTRL_median"]

# ---------------------------------------------------------------- 读表
df = pd.read_csv(IN_PATH, sep="\t", dtype=str, keep_default_na=False)
assert df.isna().sum().sum() == 0
assert not df.duplicated(["site_id", "fraction"]).any(), "site_id×fraction 不唯一"
assert (df.groupby("site_id")["protein"].nunique() == 1).all(), "site_id 对应多个 protein"
site2prot = df.drop_duplicates("site_id").set_index("site_id")["protein"]


def nonempty(col):
    return df[col] != ""


# ---------------------------------------------------------------- 行级条件
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

s3_prot = site2prot.loc[sorted(S3)]
prot_n_s3 = s3_prot.value_counts()
c5_prots = set(prot_n_s3[prot_n_s3 >= 2].index)
c5_sites = {s for s in S3 if site2prot[s] in c5_prots}

resp_sites_all = set(df.loc[df["class"].isin(RESP_CLASSES), "site_id"])
resp_in_s3 = resp_sites_all & S3
nonresp_sites = S3 - resp_sites_all
prots_with_resp = {site2prot[s] for s in resp_in_s3}
prots_with_nonresp = {site2prot[s] for s in nonresp_sites}
c6_prots = prots_with_resp & prots_with_nonresp
c6_sites = {s for s in S3 if site2prot[s] in c6_prots}

# ---------------------------------------------------------------- 漏斗（累积）
sites_1 = set(df["site_id"])
sites_2 = set(df.loc[c2, "site_id"])
sites_3 = S3
sites_4 = sites_3 & c4_sites
sites_5 = sites_4 & c5_sites
sites_6 = sites_5 & c6_sites

all_rows = pd.Series(True, index=df.index)
STEPS = [
    ("1", "全部", all_rows, sites_1),
    ("2", "行满足 c2", c2, sites_2),
    ("3", "行满足 c2&c3；位点 ∈ S3", c23, sites_3),
    ("4", "行满足 c2&c3；位点 ∈ S3 ∩ c4", c23, sites_4),
    ("5", "行满足 c2&c3；位点 ∈ ④剩余 ∩ c5", c23, sites_5),
    ("6", "行满足 c2&c3；位点 ∈ ⑤剩余 ∩ c6", c23, sites_6),
    ("ctrl_a", "③ 剩余行中 denom_tier=protein", c23 & (df["denom_tier"] == "protein"), None),
    ("ctrl_b", "③ 剩余行中 denom_tier=residue", c23 & (df["denom_tier"] == "residue"), None),
    ("ctrl_c", "宽口径：任一 log2occ_* 非空", any_occ, None),
    ("ctrl_d", "全表中 10 个数值列（log2occ_*、log2int_*）全空", all_num_empty, None),
]


def fmt_median(vals):
    if len(vals) == 0:
        return NONE_MARK
    return f"{statistics.median(vals):.4f}"


records = []
for step, cond, rowmask, siteset in STEPS:
    mask = rowmask.copy()
    if siteset is not None:
        mask &= df["site_id"].isin(siteset)
    for scope in ["all"] + FRS:
        m = mask if scope == "all" else mask & (df["fraction"] == scope)
        sub = df.loc[m]
        ctrl_vals = [float(x) for x in sub["log2int_CTRL"] if x != ""]
        records.append({
            "step": step,
            "condition": cond,
            "scope": scope,
            "n_rows": int(len(sub)),
            "n_sites": int(sub["site_id"].nunique()),
            "n_proteins": int(sub["protein"].nunique()),
            "protein_abundance_median": NONE_MARK,
            "site_log2int_CTRL_median": fmt_median(ctrl_vals),
        })
res = pd.DataFrame(records, columns=COLUMNS)
res.to_csv(OUT_TSV, sep="\t", index=False, lineterminator="\n")

# ---------------------------------------------------------------- 自检（验收标准）
def get(step, scope, col):
    r = res[(res["step"] == step) & (res["scope"] == scope)]
    assert len(r) == 1
    return int(r[col].iloc[0])


checks = []


def chk(name, ok, detail):
    checks.append((name, bool(ok), detail))


s1 = (get("1", "all", "n_rows"), get("1", "all", "n_sites"), get("1", "all", "n_proteins"))
chk("1. step1/all = 39813/18268/4996", s1 == (39813, 18268, 4996), f"实得 {s1[0]}/{s1[1]}/{s1[2]}")
v = get("2", "all", "n_rows")
chk("2a. step2/all n_rows = 14963", v == 14963, f"实得 {v}")
v = get("3", "all", "n_rows")
chk("2b. step3/all n_rows = 13106", v == 13106, f"实得 {v}")
a, b = get("ctrl_a", "all", "n_rows"), get("ctrl_b", "all", "n_rows")
chk("2c. ctrl_a = 9037", a == 9037, f"实得 {a}")
chk("2d. ctrl_b = 4069", b == 4069, f"实得 {b}")
chk("2e. ctrl_a + ctrl_b = step3/all n_rows", a + b == get("3", "all", "n_rows"), f"{a}+{b}={a + b}")
v = get("ctrl_d", "all", "n_rows")
chk("2f. ctrl_d = 16854", v == 16854, f"实得 {v}")
for col in ["n_rows", "n_sites", "n_proteins"]:
    seq = [get(str(i), "all", col) for i in range(1, 7)]
    chk(f"3. all 口径 {col} 单调不增", all(x >= y for x, y in zip(seq, seq[1:])),
        "→".join(map(str, seq)))
bad = []
for step, *_ in STEPS:
    tot = sum(get(step, fr, "n_rows") for fr in FRS)
    if tot != get(step, "all", "n_rows"):
        bad.append(f"{step}: {tot}≠{get(step, 'all', 'n_rows')}")
chk("4. 每步 6 个 fraction n_rows 之和 = all n_rows（10 个 step）", not bad, "; ".join(bad) or "全部相等")
nan_free = not res.isna().any().any() and not res.astype(str).apply(
    lambda s: s.str.lower().isin(["nan", "none", ""])).any().any()
chk("5a. 无 NaN/空值", nan_free, "")
grp = res.groupby(["step", "scope"]).size()
chk("5b. 70 行，每 step×scope 恰 1 行", len(res) == 70 and (grp == 1).all() and len(grp) == 70,
    f"行数 {len(res)}")

# 0 节基线数（附带核对）
occ_resp = df["class"].isin(RESP_CLASSES)
base = {
    "class": df["class"].value_counts().to_dict(),
    "denom_tier": df["denom_tier"].value_counts().to_dict(),
    "occ_resp": (int(occ_resp.sum()), df.loc[occ_resp, "site_id"].nunique(),
                 df.loc[occ_resp, "protein"].nunique()),
    "occ_resp_90_empty": int((occ_resp & (df["log2occ_90min"] == "")).sum()),
    "occ_resp_ctrl_empty": int((occ_resp & (df["log2occ_CTRL"] == "")).sum()),
}
chk("0. 基线 class 计数 none 32092/both 3202/occupancy_only 2447/intensity_only 2072",
    base["class"] == {"none": 32092, "both": 3202, "occupancy_only": 2447, "intensity_only": 2072},
    str(base["class"]))
chk("0. 基线 denom_tier protein 29687/residue 10126",
    base["denom_tier"] == {"protein": 29687, "residue": 10126}, str(base["denom_tier"]))
chk("0. 基线 occupancy 响应行 5649 行/3787 位点/1837 蛋白", base["occ_resp"] == (5649, 3787, 1837),
    "/".join(map(str, base["occ_resp"])))
chk("0. 基线 响应行中 log2occ_90min 空 801、log2occ_CTRL 空 0",
    (base["occ_resp_90_empty"], base["occ_resp_ctrl_empty"]) == (801, 0),
    f"{base['occ_resp_90_empty']}/{base['occ_resp_ctrl_empty']}")

# 附带计数（定义执行时的中间量，供 ④ 节）
aux = {
    "S3": len(S3),
    "c4_sites_in_S3": len(c4_sites),
    "c5_prots": len(c5_prots),
    "c5_sites_in_S3": len(c5_sites),
    "resp_sites_all": len(resp_sites_all),
    "resp_sites_in_S3": len(resp_in_s3),
    "resp_sites_not_in_S3": len(resp_sites_all - S3),
    "nonresp_sites": len(nonresp_sites),
    "c6_prots": len(c6_prots),
    "c6_sites_in_S3": len(c6_sites),
    "occresp_rows_not_c23": int((occ_resp & ~c23).sum()),
}

# ---------------------------------------------------------------- md
DEFS = """### 行级条件（在一个 位点×fraction 行上判断）
- `c2`（②）：`log2int_CTRL != ''` 且 `log2int_{tp} != ''` 对某个 tp ∈ EGF_TPS 成立。（非空 ⇔ 该时间点 ≥2 重复，07 `:295`。）
- `c3`（③"有蛋白档分母"）：`log2occ_CTRL != ''` 且 `log2occ_{tp} != ''` 对某个 tp ∈ EGF_TPS 成立。解释：occupancy 值存在 ⇔ 07 为该轨迹找到了 proteome 分母（residue 档或 protein 档均算，07 `:270-272`；residue 档存在时 protein 档必然存在，07 `:230-234`）。按 `denom_tier` 拆开报 protein / residue 两个子计数，residue 子计数即"残基档分母覆盖"对照行。
- 对照行（宽口径）：任一 `log2occ_*` 非空（不要求 CTRL+EGF 配对）。

### 位点级 / 蛋白级条件
- 记 `S3` = 在 ≥1 个 fraction 上满足 `c2 & c3` 的位点集合（③ 的位点级口径）。
- `c4`（④）：位点满足 `c2 & c3` 的 fraction 数 ≥ 2。
- `c5`（⑤）：位点所在蛋白在 `S3` 中有 ≥2 个位点。
- 位点"occupancy 响应" := 该位点任一行 `class ∈ {both, occupancy_only}`；"非响应" := 位点 ∈ S3 且无此类行。
- `c6`（⑥）：位点所在蛋白在 `S3` 中既有 ≥1 个响应位点又有 ≥1 个非响应位点。

### 漏斗（累积：每步在上一步剩余集合上再加条件）
| step | 剩余行 | 剩余位点 |
|---|---|---|
| ① 全部 | 全表 | 全表 |
| ② | 行满足 c2 | 有 ≥1 行满足 c2 的位点 |
| ③ | 行满足 c2 & c3 | S3 |
| ④ | 行满足 c2&c3 且位点满足 c4 | S3 ∩ c4 |
| ⑤ | 行满足 c2&c3 且位点 ∈ (④剩余 ∩ c5) | ④剩余 ∩ c5 |
| ⑥ | 行满足 c2&c3 且位点 ∈ (⑤剩余 ∩ c6) | ⑤剩余 ∩ c6 |

每步报：`n_rows`、`n_sites`（distinct site_id）、`n_proteins`（distinct protein）、`protein_abundance_median` = `无`、`site_log2int_CTRL_median` = 剩余行中非空 `log2int_CTRL` 的中位数（标注"位点级磷酸肽强度，非蛋白丰度"）。
**按 fraction 再报一次**：scope = FR1..FR6，行 = 该 fraction 内满足该步行级条件且位点属于该步剩余集合的行；n_sites/n_proteins 在这些行上数。
**对照行**：(a) ③ 的 protein 档子计数；(b) ③ 的 residue 档子计数；(c) 宽口径"任一 log2occ 非空"的行/位点/蛋白数；(d) 全表中 10 个数值列全空的行数。

常量：`EGF_TPS = [2min, 8min, 20min, 90min]`，`FRS = FR1..FR6`（07 `:22-25`）。读表：`pd.read_csv(path, sep='\\t', dtype=str, keep_default_na=False)`，空串 = 缺失。"""


def md_table(frame):
    lines = ["| " + " | ".join(frame.columns) + " |",
             "|" + "|".join(["---"] * len(frame.columns)) + "|"]
    for _, r in frame.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in frame.columns) + " |")
    return "\n".join(lines)


n_in = len(df)
fr_counts = df["fraction"].value_counts().reindex(FRS).fillna(0).astype(int)
md = []
md.append("# 01 漏斗表\n")
md.append("## ① 输入\n")
md.append(f"- `{IN_REL}`：{n_in} 行（不含表头），{df['site_id'].nunique()} 个 site_id，"
          f"{df['protein'].nunique()} 个 protein，{df.shape[1]} 列。")
md.append("- 各 fraction 行数：" + "，".join(f"{fr} {fr_counts[fr]}" for fr in FRS) + "。")
md.append("- 读表方式：`pd.read_csv(path, sep='\\t', dtype=str, keep_default_na=False)`。")
md.append("- 脚本：`analysis/06_pilot/01_funnel.py`；输出：`analysis/06_pilot/01_funnel.tsv`、`analysis/06_pilot/01_funnel.md`。\n")
md.append("## ② 定义（照抄 00_design.md 任务 1）\n")
md.append(DEFS + "\n")
md.append("列说明：`site_log2int_CTRL_median` = 位点级磷酸肽强度，非蛋白丰度；"
          "`protein_abundance_median` 固定为 `无`（输入中无蛋白级强度列）。\n")
md.append("## ③ 结果表（与 01_funnel.tsv 相同）\n")
md.append(md_table(res) + "\n")
md.append("### 中间计数\n")
md.append(md_table(pd.DataFrame(
    [{"量": k, "值": v} for k, v in [
        ("|S3|（位点）", aux["S3"]),
        ("S3 中满足 c4 的位点", aux["c4_sites_in_S3"]),
        ("满足 c5 的蛋白（S3 中 ≥2 位点）", aux["c5_prots"]),
        ("S3 中满足 c5 的位点", aux["c5_sites_in_S3"]),
        ("occupancy 响应位点（全表，任一行 class∈{both,occupancy_only}）", aux["resp_sites_all"]),
        ("其中 ∈ S3", aux["resp_sites_in_S3"]),
        ("其中 ∉ S3", aux["resp_sites_not_in_S3"]),
        ("非响应位点（∈ S3 且无响应行）", aux["nonresp_sites"]),
        ("满足 c6 的蛋白", aux["c6_prots"]),
        ("S3 中满足 c6 的位点", aux["c6_sites_in_S3"]),
        ("class∈{both,occupancy_only} 但不满足 c2&c3 的行", aux["occresp_rows_not_c23"]),
    ]])) + "\n")
md.append("### 自检\n")
md.append(md_table(pd.DataFrame(
    [{"项": n, "结果": "通过" if ok else "未通过", "说明": d} for n, ok, d in checks])) + "\n")
md.append("## ④ 跳过/缺失项\n")
md.append("- 蛋白级丰度：输入中无蛋白级强度列，`protein_abundance_median` 全部填 `无`。")
empty_med = res[res["site_log2int_CTRL_median"] == NONE_MARK][["step", "scope"]]
if len(empty_med):
    md.append("- `site_log2int_CTRL_median` 填 `无` 的格（该格剩余行中非空 `log2int_CTRL` 为 0 个，中位数无定义）："
              + "，".join(f"{r.step}/{r.scope}" for r in empty_med.itertuples()) + "。")
else:
    md.append("- `site_log2int_CTRL_median`：无空集格。")
md.append("- 中位数取 `statistics.median`，保留 4 位小数（输入值为 3 位小数，偶数个时取两数均值）。")
md.append("- 未使用的输入：`site_profiles.tsv`、`run_medians.tsv`、`site_features.tsv`、`kinase_name_mapping.csv`（本任务不需要）。")
with open(OUT_MD, "w", encoding="utf-8") as fh:
    fh.write("\n".join(md) + "\n")

# ---------------------------------------------------------------- 终端汇报
print(res.to_string(index=False))
print()
for k, v in aux.items():
    print(f"{k}\t{v}")
print()
for n, ok, d in checks:
    print(("PASS" if ok else "FAIL") + "\t" + n + "\t" + d)
