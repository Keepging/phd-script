#!/usr/bin/env python3
"""01_shares.py — 08_gate1 任务 1：GRB2 / SHC1 / CBL / EGFR 的区室份额表。

按 analysis/08_gate1/00_design.md §0、§2 与「任务 1」执行。
输入（只读）:
  data/pilot/gate1_proteins.tsv            列 run, genes, protein_group, quantity
  analysis/07_pilot2/01_run_factors.tsv    取 layer=proteome 的 120 行的 a 列
输出（analysis/08_gate1/）:
  01_groups.tsv, 01_shares.tsv, 01_shares_secondary.tsv, 01_shares.md

定义（照方案 §2）:
  run 解析：07 parse_design 三个正则 _(2min|8min|20min|90min|CTRL)_ (re.I)、_(FR\\d)_、_(Rep\\d)；
           timepoint ∈ {CTRL,2min,8min,20min,90min}，fraction ∈ FR1..FR6，rep ∈ Rep1..Rep4。
  基因匹配：genes 按 ';' 拆开，任一 token 等于目标基因即匹配；记录原始 genes 字符串。
  group：一个基因的全部 protein_group 都保留；primary = 测到 run 数最多者
         （并列取 protein_group 字符串排序第一，并在 md 记录并列），其余 secondary。
  quantity：转 float；空串 / 非数视为缺失。
  两版：raw = quantity；centered = quantity × 2^a（a 按 (timepoint, fraction, rep) 匹配 proteome 行）。
  份额：每个 (gene, protein_group, timepoint, rep, normalization)：
         share_f = q_f / Σ_{present} q_f；Cyt = FR1+FR2，Mem = FR3+FR4，Nuc = FR5+FR6；
         某区室两个 fraction 都缺 → 留空；6 个 fraction 全缺 → 不出行（md 记录）；不填补。

用法: python3 analysis/08_gate1/01_shares.py      （任意目录可复跑，输出确定）
      可选 --proteins / --factors / --outdir 覆盖默认路径（仅用于代码测试）。
"""
import argparse
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEF_PROT = REPO / "data" / "pilot" / "gate1_proteins.tsv"
DEF_FACT = REPO / "analysis" / "07_pilot2" / "01_run_factors.tsv"
DEF_OUT = HERE

GENES = ["GRB2", "SHC1", "CBL", "EGFR"]
TPS = ["CTRL", "2min", "8min", "20min", "90min"]
FRS = [f"FR{i}" for i in range(1, 7)]
REPS = [f"Rep{i}" for i in range(1, 5)]
COMPS = [("Cyt", ("FR1", "FR2")), ("Mem", ("FR3", "FR4")), ("Nuc", ("FR5", "FR6"))]
NORMS = ["raw", "centered"]
N_EXPECTED_RUNS = len(TPS) * len(FRS) * len(REPS)  # 120
EXPECTED_TRIPS = {(t, f, r) for t in TPS for f in FRS for r in REPS}
N_EXPECTED_ROWS = len(GENES) * len(TPS) * len(REPS) * len(NORMS)  # 160
TOL = 1e-9

IN_COLS = ["run", "genes", "protein_group", "quantity"]
FACT_COLS = ["layer", "run", "timepoint", "fraction", "rep", "a"]
GROUP_COLS = ["gene", "protein_group", "genes_raw", "n_runs_detected", "is_primary"]
SHARE_COLS = ["gene", "protein_group", "timepoint", "rep", "normalization",
              "n_fractions_present", "Cyt", "Mem", "Nuc"]
NODATA = "缺数据：`gate1_proteins.tsv` 0 行，跳过"


# ---- 照抄 code/Protein contour/Zhihan/Test/07_main_analysis.py parse_design (:32-38) ----
def parse_design(run):
    tp = re.search(r"_(2min|8min|20min|90min|CTRL)_", run, re.I)
    fr = re.search(r"_(FR\d)_", run)
    rep = re.search(r"_(Rep\d)", run)
    if not (tp and fr and rep):
        return None
    return (tp.group(1), fr.group(1), rep.group(1))


def parse_run(run):
    """返回 ((tp, fr, rep), None) 或 (None, 失败原因)。tp 统一成 TPS 中的写法（正则带 re.I）。"""
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
    s = s.strip()
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


def fmt_share(x):
    return "" if x is None else f"{x:.12f}"


def fmt4(s):
    return "" if s == "" else f"{float(s):.4f}"


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    for r in rows:
        out.append("| " + " | ".join(str(v) for v in r) + " |")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--proteins", default=str(DEF_PROT))
    ap.add_argument("--factors", default=str(DEF_FACT))
    ap.add_argument("--outdir", default=str(DEF_OUT))
    args = ap.parse_args()
    p_prot, p_fact, outdir = Path(args.proteins), Path(args.factors), Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    o_groups = outdir / "01_groups.tsv"
    o_shares = outdir / "01_shares.tsv"
    o_sec = outdir / "01_shares_secondary.tsv"
    o_md = outdir / "01_shares.md"

    # ------------------------------------------------------------ 读表
    P = read(p_prot)
    miss_cols = [c for c in IN_COLS if c not in P.columns]
    if miss_cols:
        sys.exit(f"输入缺列 {miss_cols}：{p_prot}")
    n_input = len(P)
    extra_cols = [c for c in P.columns if c not in IN_COLS]

    F = read(p_fact)
    miss_fcols = [c for c in FACT_COLS if c not in F.columns]
    if miss_fcols:
        sys.exit(f"run_factors 缺列 {miss_fcols}：{p_fact}")
    n_fact = len(F)
    FP = F[F["layer"] == "proteome"]
    a_map, fact_bad, fact_dup = {}, [], []
    for r in FP.itertuples(index=False):
        key = (r.timepoint, r.fraction, r.rep)
        a = to_float(r.a)
        if key not in EXPECTED_TRIPS or a is None:
            fact_bad.append((r.run, r.timepoint, r.fraction, r.rep, r.a))
            continue
        if key in a_map:
            fact_dup.append(key)
            continue
        a_map[key] = a

    # ------------------------------------------------------------ run 解析
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
    dup_trips = {t: rs for t, rs in trip_runs.items() if len(rs) > 1}
    n_runs_ok = len(run_trip)
    n_trips = len(trip_runs)
    trips_no_a = sorted(t for t in trip_runs if t not in a_map)

    # ------------------------------------------------------------ 基因匹配
    E = P.reset_index(drop=False).rename(columns={"index": "x_row"})
    E = E.assign(x_tok=E["genes"].str.split(";")).explode("x_tok")
    E["x_tok"] = E["x_tok"].fillna("").astype(str).str.strip()
    M = E[E["x_tok"].isin(GENES)].drop_duplicates(["x_row", "x_tok"])
    n_rows_matched = M["x_row"].nunique()
    n_rows_matched_by_gene = {g: int((M["x_tok"] == g).sum()) for g in GENES}

    groups_raw = defaultdict(set)                    # (gene, group) -> {genes 原串}
    cell_vals = defaultdict(list)                    # (gene, group, tp, fr, rep) -> [(run, q)]
    n_unparsed_matched, n_q_missing, nonpos = 0, 0, []
    for r in M.itertuples(index=False):
        gene, grp = r.x_tok, r.protein_group
        groups_raw[(gene, grp)].add(r.genes)
        q = to_float(r.quantity)
        if q is None:
            n_q_missing += 1
        elif q <= 0:
            nonpos.append((gene, grp, r.run, r.quantity))
        trip = run_trip.get(r.run)
        if trip is None:
            n_unparsed_matched += 1
            continue
        cell_vals[(gene, grp) + trip].append((r.run, q))

    # 同一 (gene, group, tp, fr, rep) 多条有效 quantity：值全同取该值，否则冲突→缺失（方案未定义）
    cell_q, conflicts = {}, []
    detected = defaultdict(set)                      # (gene, group) -> {trip}
    for key, lst in cell_vals.items():
        qs = [q for _, q in lst if q is not None]
        if not qs:
            continue
        detected[key[:2]].add(key[2:])
        if len(set(qs)) == 1:
            cell_q[key] = qs[0]
        else:
            conflicts.append((key, lst))
    gene_detected = {g: set().union(*[detected[k] for k in groups_raw if k[0] == g] or [set()])
                     for g in GENES}

    # ------------------------------------------------------------ group 表
    group_rows, primary, ties = [], {}, []
    for g in GENES:
        grps = sorted(k[1] for k in groups_raw if k[0] == g)
        if not grps:
            if n_input > 0:  # 0 行输入时表只留表头（方案验收），见 md ②
                group_rows.append([g, "", "", 0, ""])
            continue
        nmax = max(len(detected[(g, x)]) for x in grps)
        tops = [x for x in grps if len(detected[(g, x)]) == nmax]
        primary[g] = tops[0]
        if len(tops) > 1:
            ties.append((g, tops, nmax))
        for x in grps:
            group_rows.append([g, x, " | ".join(sorted(groups_raw[(g, x)])),
                               len(detected[(g, x)]), "yes" if x == primary[g] else "no"])
    G = pd.DataFrame(group_rows, columns=GROUP_COLS)

    # ------------------------------------------------------------ 份额
    def value(g, x, tp, fr, rep, norm):
        q = cell_q.get((g, x, tp, fr, rep))
        if q is None:
            return None
        if norm == "raw":
            return q
        a = a_map.get((tp, fr, rep))
        return None if a is None else q * 2.0 ** a

    share_rows = {"primary": [], "secondary": []}
    missing_list, empty_list, zero_total, out_of_01 = [], [], [], []
    for g in GENES:
        grps = sorted(k[1] for k in groups_raw if k[0] == g)
        for x in grps:
            role = "primary" if x == primary.get(g) else "secondary"
            if not detected[(g, x)]:
                empty_list.append((g, x, role, "全部", "全部", "全部", "该 group 0 run 测到"))
                continue
            for tp in TPS:
                for rep in REPS:
                    for norm in NORMS:
                        vals = {fr: value(g, x, tp, fr, rep, norm) for fr in FRS}
                        present = {fr: v for fr, v in vals.items() if v is not None}
                        n = len(present)
                        if n == 0:
                            empty_list.append((g, x, role, tp, rep, norm, "6 个 fraction 全缺，不出行"))
                            continue
                        tot = sum(present.values())
                        comp = {}
                        for cname, frs in COMPS:
                            have = [present[f] for f in frs if f in present]
                            comp[cname] = None if (not have or tot == 0) else sum(have) / tot
                        if tot == 0:
                            zero_total.append((g, x, tp, rep, norm))
                        if any(v is not None and not (0 <= v <= 1) for v in comp.values()):
                            out_of_01.append((g, x, tp, rep, norm))
                        if n < 6:
                            miss = [fr for fr in FRS if fr not in present]
                            empty_c = [c for c, frs in COMPS if all(f not in present for f in frs)]
                            missing_list.append((g, x, role, tp, rep, norm, n,
                                                 ",".join(miss), ",".join(empty_c)))
                        share_rows[role].append([g, x, tp, rep, norm, n,
                                                 fmt_share(comp["Cyt"]), fmt_share(comp["Mem"]),
                                                 fmt_share(comp["Nuc"])])
    S1 = pd.DataFrame(share_rows["primary"], columns=SHARE_COLS)
    S2 = pd.DataFrame(share_rows["secondary"], columns=SHARE_COLS)

    G.to_csv(o_groups, sep="\t", index=False)
    S1.to_csv(o_shares, sep="\t", index=False)
    S2.to_csv(o_sec, sep="\t", index=False)

    # ------------------------------------------------------------ 自检（读回写出的文件）
    checks = []

    def chk(name, ok, val=""):
        status = ok if isinstance(ok, str) else ("PASS" if ok else "FAIL")
        checks.append((name, status, val))

    R_G, R_1, R_2 = read(o_groups), read(o_shares), read(o_sec)
    chk("01_groups.tsv 列 = 方案列", list(R_G.columns) == GROUP_COLS)
    chk("01_shares.tsv / 01_shares_secondary.tsv 列 = 方案列",
        list(R_1.columns) == SHARE_COLS and list(R_2.columns) == SHARE_COLS)
    chk("run_factors proteome 行 = 120，(tp, fr, rep) 唯一且 a 全为数值",
        len(FP) == N_EXPECTED_RUNS and len(a_map) == N_EXPECTED_RUNS and not fact_bad and not fact_dup,
        f"{len(FP)} 行 / 有效 a {len(a_map)}")
    for tag, R in (("01_shares", R_1), ("01_shares_secondary", R_2)):
        if len(R) == 0:
            chk(f"{tag}：Cyt+Mem+Nuc = 1 ± 1e-9（三格非空的行）", "不适用", "0 行")
            chk(f"{tag}：raw 与 centered 行数相同", (R["normalization"] == "raw").sum()
                == (R["normalization"] == "centered").sum(), "0 / 0")
            chk(f"{tag}：n_fractions_present ∈ 1..6", "不适用", "0 行")
            continue
        full = R[(R.Cyt != "") & (R.Mem != "") & (R.Nuc != "")]
        dev = [abs(float(a) + float(b) + float(c) - 1) for a, b, c in zip(full.Cyt, full.Mem, full.Nuc)]
        chk(f"{tag}：Cyt+Mem+Nuc = 1 ± 1e-9（三格非空的行）", all(d <= TOL for d in dev),
            f"{len(full)} 行，max|偏差| = {max(dev) if dev else 0:.2e}")
        part = R[[any(v == "" for v in t) and any(v != "" for v in t)
                  for t in zip(R.Cyt, R.Mem, R.Nuc)]]
        devp = [abs(sum(float(v) for v in t if v != "") - 1) for t in zip(part.Cyt, part.Mem, part.Nuc)]
        chk(f"{tag}：有区室留空的行，非空区室之和 = 1 ± 1e-9", all(d <= TOL for d in devp),
            f"{len(part)} 行")
        nr, nc = (R.normalization == "raw").sum(), (R.normalization == "centered").sum()
        chk(f"{tag}：raw 与 centered 行数相同", nr == nc, f"{nr} / {nc}")
        nf = R.n_fractions_present.astype(int)
        chk(f"{tag}：n_fractions_present ∈ 1..6", bool(((nf >= 1) & (nf <= 6)).all()),
            f"{nf.min()}..{nf.max()}")
    if n_input == 0:
        chk("0 行输入：01_groups / 01_shares / 01_shares_secondary 只有表头",
            len(R_G) == 0 and len(R_1) == 0 and len(R_2) == 0, f"{len(R_G)} / {len(R_1)} / {len(R_2)}")
    else:
        chk("0 行输入：所有表只有表头", "不适用", f"输入 {n_input} 行")
        chk("01_shares.tsv 行数 = 160（数据齐全时）", "PASS" if len(R_1) == N_EXPECTED_ROWS else "不齐",
            f"{len(R_1)}")
        chk("每个目标基因在 01_groups.tsv 至少一行", set(R_G.gene) == set(GENES))
        chk("每个测到的基因恰一个 primary",
            all((R_G[(R_G.gene == g)].is_primary == "yes").sum() == (1 if g in primary else 0)
                for g in GENES))

    # ------------------------------------------------------------ md
    L = []
    L.append("# 08 任务 1：GRB2 / SHC1 / CBL / EGFR 区室份额表\n")
    L.append(f"脚本：`analysis/08_gate1/01_shares.py`（`python3`，任意目录可复跑，输出确定）。"
             f"只报计数和数值。\n")

    # ① 输入
    L.append("## ① 输入\n")
    L.append(f"- `{rel(p_prot)}`：{n_input} 行（不含表头），列 `{', '.join(P.columns)}`；"
             f"读法 `pd.read_csv(sep='\\t', dtype=str, keep_default_na=False)`，空串 = 缺失。"
             + (f" 方案外的列（不使用）：`{', '.join(extra_cols)}`。" if extra_cols else ""))
    L.append(f"- `{rel(p_fact)}`：{n_fact} 行，layer=proteome {len(FP)} 行，"
             f"有效 (timepoint, fraction, rep) → a：{len(a_map)} 个"
             + (f"；无效行 {len(fact_bad)}、重复三元组 {len(fact_dup)}" if (fact_bad or fact_dup) else "")
             + "。")
    L.append(f"- distinct run：{len(runs)}；解析成功 {n_runs_ok}/{N_EXPECTED_RUNS}；"
             f"解析失败 {len(parse_fail)}；覆盖的 (timepoint, fraction, rep) 三元组 "
             f"{n_trips}/{N_EXPECTED_RUNS}；多个 run 映射到同一三元组：{len(dup_trips)}；"
             f"在 run_factors 中找不到 a 的三元组：{len(trips_no_a)}。")
    if parse_fail:
        L.append("\n解析失败的 run：\n")
        L.append(md_table(["run", "原因"], [(f"`{r}`", w) for r, w in parse_fail]))
        L.append("")
    if dup_trips:
        L.append("\n多个 run 映射到同一三元组：\n")
        L.append(md_table(["timepoint", "fraction", "rep", "runs"],
                          [t + ("; ".join(f"`{x}`" for x in rs),) for t, rs in sorted(dup_trips.items())]))
        L.append("")
    if trips_no_a:
        L.append("\n找不到 a 的三元组（centered 版该 run 当缺失）：" +
                 "、".join(f"{t[0]}/{t[1]}/{t[2]}" for t in trips_no_a) + "\n")
    L.append(f"- 匹配目标基因的行：{n_rows_matched}（按基因："
             + "、".join(f"{g} {n_rows_matched_by_gene[g]}" for g in GENES)
             + f"）；其中 run 解析失败 {n_unparsed_matched}、quantity 缺失（空串/非数）{n_q_missing}、"
               f"quantity ≤ 0 {len(nonpos)}。")
    L.append("- **四个蛋白各在多少 run 测到（/120，任一 protein_group 有有效 quantity 的三元组数）**："
             + "，".join(f"{g} {len(gene_detected[g])}/{N_EXPECTED_RUNS}" for g in GENES) + "。")
    L.append("\ngroup 表（`01_groups.tsv`）：\n")
    if n_input == 0:
        L.append(NODATA + "（`01_groups.tsv` 只有表头）。\n")
    else:
        L.append(md_table(GROUP_COLS, G.values.tolist()))
        L.append("")
        if ties:
            L.append("primary 并列（取 protein_group 字符串排序第一）：" +
                     "；".join(f"{g}：{', '.join(t)}（各 {n} run）" for g, t, n in ties) + "\n")
        else:
            L.append("primary 并列：0 个基因。\n")

    # ② 定义
    L.append("## ② 定义\n")
    L.append("照抄 `analysis/08_gate1/00_design.md` §2 / 任务 1：\n")
    L.append("- 蛋白集：GRB2、SHC1、CBL、EGFR。匹配 `genes` 列（分号分隔的多基因按 `;` 拆开后任一等于目标基因即匹配，记录原始 `genes` 字符串）。")
    L.append("- 一个基因对应多个 `protein_group`：全部保留；`primary` = 出现 run 数最多的 group（并列取 `protein_group` 字符串排序第一，并记录并列）；其余为 `secondary`，写附表。")
    L.append("- run 解析：`timepoint ∈ {CTRL,2min,8min,20min,90min}`，`fraction ∈ FR1..FR6`，`rep ∈ Rep1..Rep4`；期望 120 个 proteome run；解析失败的 run 计数并列出。正则为 07 `parse_design`（`code/Protein contour/Zhihan/Test/07_main_analysis.py:32-38`）的 `_(2min|8min|20min|90min|CTRL)_`（re.I）、`_(FR\\d)_`、`_(Rep\\d)`。")
    L.append("- `quantity` 转 float；空串/非数视为缺失。")
    L.append("- 两版归一化：`raw` = quantity；`centered` = quantity × 2^a（a 取 `01_run_factors.tsv` 中 layer=proteome、同 run 的 `a`）。run 名匹配用 (timepoint, fraction, rep) 三元组。")
    L.append("- 份额：对每个 (gene, protein_group, timepoint, rep)：`n_fractions_present` = 6 个 fraction 中有值的个数；`share_f = q_f / Σ_{present} q_f`；`Cyt = share_FR1 + share_FR2`、`Mem = share_FR3 + share_FR4`、`Nuc = share_FR5 + share_FR6`；某区室两个 fraction 都缺 → 该区室格留空；`n_fractions_present < 6` 的行单列\"缺 fraction 清单\"；不填补。raw、centered 各算一次（`normalization` 列）。")
    L.append("- `01_shares.tsv`（primary）/ `01_shares_secondary.tsv`（其余 group）列：`" + ", ".join(SHARE_COLS) + "`；某 (tp, rep) 一个 fraction 都没有则不出行，并在 md 记录。")
    L.append("- `01_groups.tsv` 列：`" + ", ".join(GROUP_COLS) + "`（每个 gene×group 一行；gene 没测到则一行 `n_runs_detected = 0`）。")
    L.append("\n脚本补充规则（方案未写明处，非替换定义）：\n")
    L.append("- timepoint 正则带 re.I，匹配到的写法统一成 `CTRL/2min/8min/20min/90min`；`FR\\d`、`Rep\\d` 匹配到但不在 FR1..FR6 / Rep1..Rep4 的 run 记为解析失败（原因\"超出设计范围\"）。")
    L.append("- 基因 token 去首尾空白后与目标基因做区分大小写的完全相等比较；同一行 genes 中重复的同一目标基因只计一次。")
    L.append("- `nan`/`inf` 字符串按\"非数\"处理为缺失；quantity ≤ 0 保留原值参与计算，计数见 ①，若有份额落在 [0,1] 外在 ⑥ 列出；某 (tp, rep) 有值 fraction 之和 = 0 时三格留空并在 ⑥ 列出。")
    L.append("- `n_runs_detected` = 该 gene×group 有 ≥1 条有效 quantity 的 (timepoint, fraction, rep) 三元组数（只计解析成功的 run）；基因级\"测到 run 数\"= 该基因任一 group 有有效 quantity 的三元组数。")
    L.append("- `genes_raw`：该 gene×group 出现过的全部原始 `genes` 字符串，排序后以 ` | ` 连接；`is_primary` 取 `yes`/`no`（gene 没测到的那一行留空）。")
    L.append("- 同一 (gene, protein_group, timepoint, fraction, rep) 有多条有效 quantity：值全相同取该值；值不同记为冲突，该格当缺失并在 ⑥ 列出。")
    L.append("- 方案两处对 0 行输入的要求不一致（`01_groups.tsv`\"gene 没测到则一行 n_runs_detected = 0\" vs 验收\"0 行输入时所有表只有表头\"）：本脚本在输入 0 行时按验收写只有表头的 `01_groups.tsv`，四个基因的 0/120 写在 ①；输入非 0 行时，没测到的基因各写一行 `n_runs_detected = 0`。")
    L.append("- tsv 中份额保留 12 位小数；md 表中保留 4 位。行序：gene（GRB2, SHC1, CBL, EGFR）→ protein_group → timepoint（CTRL, 2min, 8min, 20min, 90min）→ rep → normalization（raw, centered）。\n")

    # ③ 缺 fraction 清单
    L.append("## ③ 缺 fraction 清单\n")
    if n_input == 0:
        L.append(NODATA + "。\n")
    else:
        L.append(f"`n_fractions_present < 6` 的行：{len(missing_list)}（primary "
                 f"{sum(1 for m in missing_list if m[2] == 'primary')}、secondary "
                 f"{sum(1 for m in missing_list if m[2] == 'secondary')}）。\n")
        if missing_list:
            L.append(md_table(["gene", "protein_group", "role", "timepoint", "rep", "normalization",
                               "n_fractions_present", "缺的 fraction", "留空区室"], missing_list))
            L.append("")
        L.append(f"\n一个 fraction 都没有、因而不出行的 (gene, group, tp, rep, normalization)：{len(empty_list)} 条。\n")
        if empty_list:
            L.append(md_table(["gene", "protein_group", "role", "timepoint", "rep", "normalization", "说明"],
                              empty_list))
            L.append("")

    # ④ 份额表（primary，两版）
    L.append("## ④ 份额表（primary，raw 与 centered 并列）\n")
    if n_input == 0:
        L.append(NODATA + "（`01_shares.tsv` 只有表头）。\n")
    else:
        L.append(f"`01_shares.tsv`：{len(S1)} 行（raw {int((S1.normalization == 'raw').sum())}、"
                 f"centered {int((S1.normalization == 'centered').sum())}；数据齐全时应为 {N_EXPECTED_ROWS}）。"
                 "按基因：" + "、".join(f"{g} {int((S1.gene == g).sum())}" for g in GENES) + "。\n")
        nod = [g for g in GENES if g not in primary or not detected[(g, primary[g])]]
        if nod:
            L.append("无份额行的基因（0/120 run 测到）：" + "、".join(nod) + "。\n")
        if len(S1):
            idx = {}
            for r in S1.itertuples(index=False):
                idx[(r.gene, r.protein_group, r.timepoint, r.rep, r.normalization)] = r
            hdr = ["gene", "protein_group", "timepoint", "rep",
                   "n_fr raw", "Cyt raw", "Mem raw", "Nuc raw",
                   "n_fr centered", "Cyt centered", "Mem centered", "Nuc centered"]
            rows = []
            seen = []
            for r in S1.itertuples(index=False):
                k = (r.gene, r.protein_group, r.timepoint, r.rep)
                if k in seen:
                    continue
                seen.append(k)
                row = list(k)
                for norm in NORMS:
                    x = idx.get(k + (norm,))
                    row += (["", "", "", ""] if x is None else
                            [x.n_fractions_present, fmt4(x.Cyt), fmt4(x.Mem), fmt4(x.Nuc)])
                rows.append(row)
            L.append(md_table(hdr, rows))
            L.append("")

    # ⑤ 附表说明
    L.append("## ⑤ 附表说明\n")
    if n_input == 0:
        L.append(NODATA + "（`01_shares_secondary.tsv`、`01_groups.tsv` 只有表头）。\n")
    else:
        L.append("- `01_groups.tsv`：每个 gene×protein_group 一行，"
                 f"{len(G)} 行；列 `{', '.join(GROUP_COLS)}`。")
        L.append(f"- `01_shares_secondary.tsv`：非 primary group 的份额，列同 `01_shares.tsv`，{len(S2)} 行。")
        sec = [(g, x) for g, x in sorted(groups_raw, key=lambda k: (GENES.index(k[0]), k[1]))
               if x != primary.get(g)]
        if sec:
            L.append("\n| gene | protein_group | n_runs_detected | rows raw | rows centered |")
            L.append("|---|---|---|---|---|")
            for g, x in sec:
                sub = S2[(S2.gene == g) & (S2.protein_group == x)]
                L.append(f"| {g} | {x} | {len(detected[(g, x)])} | "
                         f"{int((sub.normalization == 'raw').sum())} | "
                         f"{int((sub.normalization == 'centered').sum())} |")
            L.append("")
        else:
            L.append("- secondary group：0 个。\n")

    # ⑥ 跳过项
    L.append("## ⑥ 跳过项 / 记录\n")
    if n_input == 0:
        L.append(f"- `{rel(p_prot)}` 0 行数据：③ 缺 fraction 清单、④ 份额表、⑤ 附表全部跳过；"
                 "三张 tsv 照常生成，只含表头。数据到位后直接复跑本脚本即可。")
    L.append(f"- 解析失败 run：{len(parse_fail)}；多 run 同三元组：{len(dup_trips)}；"
             f"找不到 a 的三元组：{len(trips_no_a)}。")
    L.append(f"- 值冲突的 (gene, group, tp, fraction, rep) 格：{len(conflicts)}"
             + ("：" + "；".join(f"{k[0]}/{k[1]}/{k[2]}/{k[3]}/{k[4]}" for k, _ in conflicts) if conflicts else "")
             + "。")
    L.append(f"- 有值 fraction 之和 = 0（三格留空）的行：{len(zero_total)}"
             + ("：" + "；".join("/".join(z) for z in zero_total) if zero_total else "") + "。")
    L.append(f"- 份额落在 [0,1] 外的行：{len(out_of_01)}"
             + ("：" + "；".join("/".join(z) for z in out_of_01) if out_of_01 else "") + "。")
    L.append("- 不复现 DAPAR 的 LOESS 归一化与 KNN 填补（方案 §1），只用 raw 与 centered 两版；不填补缺失。")
    L.append("- 代码冒烟测试：脚本在 scratchpad 合成夹具上跑通（夹具数据与结果均未写入本目录）。\n")
    L.append("### 自检（方案验收项）\n")
    L.append(md_table(["检查", "结果", "值"], checks))
    L.append("")

    o_md.write_text("\n".join(L) + "\n", encoding="utf-8")

    # ------------------------------------------------------------ stdout
    print(f"input rows {n_input}; runs parsed {n_runs_ok}/{N_EXPECTED_RUNS}; parse failures {len(parse_fail)}")
    print("detected runs: " + ", ".join(f"{g} {len(gene_detected[g])}/{N_EXPECTED_RUNS}" for g in GENES))
    print(f"rows: 01_groups {len(G)}, 01_shares {len(S1)}, 01_shares_secondary {len(S2)}")
    for name, st, val in checks:
        print(f"  [{st}] {name} {val}")
    print("wrote:", ", ".join(rel(p) for p in (o_groups, o_shares, o_sec, o_md)))
    return 0 if all(st != "FAIL" for _, st, _ in checks) else 1


if __name__ == "__main__":
    sys.exit(main())
