#!/usr/bin/env python3
"""03_egfr_total.py — 08_gate1 任务 3：EGFR 总量。

按 analysis/08_gate1/00_design.md 的 §0（输入现状）、§2（共用定义）与「任务 3：EGFR 总量」执行。

输入（只读）:
  data/pilot/gate1_proteins.tsv                 列 run, genes, protein_group, quantity
  analysis/07_pilot2/01_run_factors.tsv         取 layer=proteome 的 120 行的 a 列
输出（analysis/08_gate1/）:
  03_egfr_total.tsv   gene, protein_group, is_primary, timepoint, rep, normalization,
                      n_fractions_present, total_quantity
  03_egfr_ratio.tsv   gene, protein_group, normalization, rep, total_90min, total_CTRL, ratio
  03_egfr_total.md    ① 输入 ② 定义 ③ 总量表 ④ 90 min / CTRL 比值表 ⑤ 跳过项

定义（照方案）：
  - genes 按 ";" 拆开（去首尾空白）后任一项 == "EGFR" 即匹配，记录原始 genes 字符串。
  - 一个基因多个 protein_group 全部保留；primary = 测到 run 数最多的 group
    （并列取 protein_group 字符串排序第一，并记录并列）；其余 secondary（附表）。
  - run 解析：07 parse_design 三个正则（照抄，:32-38）；timepoint 正则为大小写不敏感，
    匹配结果统一为 CTRL/2min/8min/20min/90min；fraction 须 ∈ FR1..FR6、rep 须 ∈ Rep1..Rep4，
    否则计为解析失败并列出。期望 120 个 proteome run。
  - quantity 转 float；空串 / 非数 / 非有限值视为缺失。
  - raw = quantity；centered = quantity × 2^a（a 取 01_run_factors.tsv 中 layer=proteome、
    同 (timepoint, fraction, rep) 的 a）。
  - 每个 (protein_group, timepoint, rep, normalization)：
      n_fractions_present = 6 个 fraction 中有值的个数；total = Σ_present q_f；不填补；
      n_fractions_present = 0 时 total 留空。
  - ratio_r = total_90min,r / total_CTRL,r，按 rep 编号配对；任一侧 total 缺则留空。

用法:
  python3 analysis/08_gate1/03_egfr_total.py            （任意目录可复跑，输出确定）
  可选参数（仅供代码测试，默认即仓库路径）:
    --gate1 PATH  --factors PATH  --outdir DIR
"""
import argparse
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
DEF_GATE1 = REPO / "data" / "pilot" / "gate1_proteins.tsv"
DEF_FACTORS = REPO / "analysis" / "07_pilot2" / "01_run_factors.tsv"
DEF_OUTDIR = REPO / "analysis" / "08_gate1"

GENE = "EGFR"
TPS = ["CTRL", "2min", "8min", "20min", "90min"]
FRS = [f"FR{i}" for i in range(1, 7)]
REPS = [f"Rep{i}" for i in range(1, 5)]
NORMS = ["raw", "centered"]
N_DESIGN = len(TPS) * len(FRS) * len(REPS)  # 120

IN_COLS = ["run", "genes", "protein_group", "quantity"]
TOTAL_COLS = ["gene", "protein_group", "is_primary", "timepoint", "rep", "normalization",
              "n_fractions_present", "total_quantity"]
RATIO_COLS = ["gene", "protein_group", "normalization", "rep", "total_90min", "total_CTRL", "ratio"]

SKIP_MSG = "缺数据：`gate1_proteins.tsv` 0 行，跳过"


# ---- 照抄 code/Protein contour/Zhihan/Test/07_main_analysis.py parse_design (:32-38) ----
def parse_design(run):
    tp = re.search(r"_(2min|8min|20min|90min|CTRL)_", run, re.I)
    fr = re.search(r"_(FR\d)_", run)
    rep = re.search(r"_(Rep\d)", run)
    if not (tp and fr and rep):
        return None
    return (tp.group(1), fr.group(1), rep.group(1))
# ---- 照抄结束 ----


TP_CANON = {t.lower(): t for t in TPS}


def parse_run(run):
    """返回 ((tp, fr, rep), None) 或 (None, 失败原因)。"""
    p = parse_design(run)
    if p is None:
        return None, "三个正则未全部匹配"
    tp, fr, rep = p
    tp = TP_CANON[tp.lower()]
    if fr not in FRS:
        return None, f"fraction {fr} 不在 FR1..FR6"
    if rep not in REPS:
        return None, f"rep {rep} 不在 Rep1..Rep4"
    return (tp, fr, rep), None


def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def to_float(s):
    """空串/非数/非有限值 -> None。"""
    s = s.strip()
    if s == "":
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    if not math.isfinite(v):
        return None
    return v


def gene_match(genes_str, target=GENE):
    return any(g.strip() == target for g in genes_str.split(";"))


def fnum(v):
    """tsv 数值：最短可往返表示；None -> 空串。"""
    return "" if v is None else repr(float(v))


def g6(v):
    return "" if v is None else f"{v:.6g}"


def f4(v):
    return "" if v is None else f"{v:.4f}"


def md_cell(s):
    return str(s).replace("|", "\\|")


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    for r in rows:
        out.append("| " + " | ".join(md_cell(x) for x in r) + " |")
    return out


def write_tsv(path, cols, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")


def rel(p):
    p = Path(p).resolve()
    try:
        return str(p.relative_to(REPO))
    except ValueError:
        return str(p)


def main():
    ap = argparse.ArgumentParser(description="08_gate1 任务 3：EGFR 总量")
    ap.add_argument("--gate1", default=str(DEF_GATE1))
    ap.add_argument("--factors", default=str(DEF_FACTORS))
    ap.add_argument("--outdir", default=str(DEF_OUTDIR))
    args = ap.parse_args()
    gate1_path, factors_path, outdir = Path(args.gate1), Path(args.factors), Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_total = outdir / "03_egfr_total.tsv"
    out_ratio = outdir / "03_egfr_ratio.tsv"
    out_md = outdir / "03_egfr_total.md"

    # ------------------------------------------------------------------ 读入
    G = read_tsv(gate1_path)
    miss_cols = [c for c in IN_COLS if c not in G.columns]
    if miss_cols:
        sys.exit(f"缺列 {miss_cols}：{gate1_path}")
    n_in = len(G)
    in_cols = list(G.columns)

    F = read_tsv(factors_path)
    Fp = F[F["layer"] == "proteome"]
    factor = {}
    factor_bad = []
    for _, r in Fp.iterrows():
        key = (r["timepoint"], r["fraction"], r["rep"])
        a = to_float(r["a"])
        if a is None or key in factor:
            factor_bad.append(key)
            continue
        factor[key] = a
    n_factor_rows = len(Fp)

    # ------------------------------------------------------------------ run 解析
    run_strings = sorted(set(G["run"]))
    run_key, parse_fail = {}, []
    for rs in run_strings:
        k, why = parse_run(rs)
        if k is None:
            parse_fail.append((rs, why))
        else:
            run_key[rs] = k
    key_to_runs = defaultdict(list)
    for rs, k in run_key.items():
        key_to_runs[k].append(rs)
    n_design_cells = len(key_to_runs)
    multi_run_cells = {k: v for k, v in key_to_runs.items() if len(v) > 1}

    # ------------------------------------------------------------------ EGFR 行
    G = G.assign(_match=G["genes"].map(gene_match))
    E = G[G["_match"]]
    n_egfr_rows = len(E)
    egfr_parse_fail_rows = 0
    egfr_missing_q_rows = 0
    egfr_nonpos_rows = 0
    genes_raw = defaultdict(set)
    # obs[(group, tp, fr, rep)] -> list of (quantity, run)
    obs = defaultdict(list)
    for _, r in E.iterrows():
        grp = r["protein_group"]
        genes_raw[grp].add(r["genes"])
        k = run_key.get(r["run"])
        if k is None:
            egfr_parse_fail_rows += 1
            continue
        q = to_float(r["quantity"])
        if q is None:
            egfr_missing_q_rows += 1
            continue
        if q <= 0:
            egfr_nonpos_rows += 1
        obs[(grp,) + k].append((q, r["run"]))

    # 同一 (group, tp, fr, rep) 多条有值记录：值相同计一次；值不同 -> 冲突，该 fraction 不计入
    qval, conflicts, dup_same = {}, [], 0
    for key, lst in obs.items():
        vals = sorted(set(q for q, _ in lst))
        if len(vals) == 1:
            qval[key] = vals[0]
            if len(lst) > 1:
                dup_same += 1
        else:
            conflicts.append((key, lst))

    # group 统计：测到 run 数（解析成功、quantity 有值、无冲突的 (tp, fr, rep) 数）
    groups = sorted(genes_raw)
    n_runs = {g: 0 for g in groups}
    for (grp, tp, fr, rep) in qval:
        n_runs[grp] += 1
    egfr_cells = {(tp, fr, rep) for (grp, tp, fr, rep) in qval}
    n_egfr_runs = len(egfr_cells)

    primary, tie = None, []
    if groups:
        best = max(n_runs.values())
        tied = sorted(g for g in groups if n_runs[g] == best)
        primary = tied[0]
        tie = tied if len(tied) > 1 else []
    order = ([primary] if primary is not None else []) + sorted(
        [g for g in groups if g != primary], key=lambda g: (-n_runs[g], g))

    # ------------------------------------------------------------------ 总量
    def qn(grp, tp, fr, rep, norm):
        q = qval.get((grp, tp, fr, rep))
        if q is None:
            return None
        if norm == "raw":
            return q
        a = factor.get((tp, fr, rep))
        return None if a is None else q * (2.0 ** a)

    total_rows = []
    totals = {}       # (grp, norm, tp, rep) -> (n_present, total or None, missing frs)
    centered_no_factor = set()
    for grp in order:
        for norm in NORMS:
            for tp in TPS:
                for rep in REPS:
                    vals, missing = [], []
                    for fr in FRS:
                        v = qn(grp, tp, fr, rep, norm)
                        if v is None:
                            missing.append(fr)
                            if norm == "centered" and (grp, tp, fr, rep) in qval:
                                centered_no_factor.add((tp, fr, rep))
                        else:
                            vals.append(v)
                    n = len(vals)
                    tot = math.fsum(vals) if n > 0 else None
                    totals[(grp, norm, tp, rep)] = (n, tot, missing)
                    total_rows.append({
                        "gene": GENE, "protein_group": grp,
                        "is_primary": "yes" if grp == primary else "no",
                        "timepoint": tp, "rep": rep, "normalization": norm,
                        "n_fractions_present": n, "total_quantity": fnum(tot)})

    ratio_rows = []
    ratios = {}       # (grp, norm, rep) -> (t90, tctrl, ratio, note)
    for grp in order:
        for norm in NORMS:
            for rep in REPS:
                t90 = totals[(grp, norm, "90min", rep)][1]
                tct = totals[(grp, norm, "CTRL", rep)][1]
                note = ""
                if t90 is None or tct is None:
                    ratio = None
                    note = "缺 " + "、".join(
                        s for s, t in (("90min", t90), ("CTRL", tct)) if t is None)
                elif tct == 0:
                    ratio = None
                    note = "CTRL total = 0"
                else:
                    ratio = t90 / tct
                ratios[(grp, norm, rep)] = (t90, tct, ratio, note)
                ratio_rows.append({
                    "gene": GENE, "protein_group": grp, "normalization": norm, "rep": rep,
                    "total_90min": fnum(t90), "total_CTRL": fnum(tct), "ratio": fnum(ratio)})

    write_tsv(out_total, TOTAL_COLS, total_rows)
    write_tsv(out_ratio, RATIO_COLS, ratio_rows)

    # ------------------------------------------------------------------ 自检
    checks = []
    # total = Σ 有值 fraction（从 qval / factor 独立重算）
    ok = True
    for (grp, norm, tp, rep), (n, tot, missing) in totals.items():
        vs = []
        for fr in FRS:
            q = qval.get((grp, tp, fr, rep))
            if q is None:
                continue
            if norm == "centered":
                a = factor.get((tp, fr, rep))
                if a is None:
                    continue
                q = q * 2.0 ** a
            vs.append(q)
        if len(vs) != n or (n and abs(sum(vs) - tot) > 1e-9 * max(1.0, abs(tot))) or \
                (n == 0 and tot is not None) or len(missing) != 6 - n:
            ok = False
    checks.append(("total = Σ 有值 fraction；n_fractions_present + 缺 fraction 数 = 6", ok))
    checks.append(("n_fractions_present ∈ 0..6",
                   all(0 <= r["n_fractions_present"] <= 6 for r in total_rows)))
    nr = sum(r["normalization"] == "raw" for r in total_rows)
    nc = sum(r["normalization"] == "centered" for r in total_rows)
    checks.append(("raw 与 centered 行数相同", nr == nc))
    checks.append(("总量表每个 group 40 行（5 tp × 4 rep × 2 版）",
                   len(total_rows) == 40 * len(order)))
    checks.append(("比值表每个 group 8 行（4 rep × 2 版）", len(ratio_rows) == 8 * len(order)))
    checks.append(("primary 恰 1 个（有 EGFR group 时）",
                   sum(1 for g in order if g == primary) == (1 if order else 0)))
    for name, res in checks:
        print(f"[check] {'PASS' if res else 'FAIL'}  {name}")

    # ------------------------------------------------------------------ md
    L = []
    L.append("# 08_gate1 任务 3：EGFR 总量")
    L.append("")
    L.append("生成脚本：`analysis/08_gate1/03_egfr_total.py`（按 `analysis/08_gate1/00_design.md` "
             "§0、§2、任务 3 执行）。")
    L.append("")
    L.append("## ① 输入")
    L.append("")
    L.append(f"- `{rel(gate1_path)}`：数据行 {n_in}（不含表头）；列 `{', '.join(G.columns[:-1])}`。")
    L.append(f"- `{rel(factors_path)}`：layer=proteome {n_factor_rows} 行，"
             f"可用 (timepoint, fraction, rep) → a {len(factor)} 个；异常（a 非数或三元组重复）"
             f"{len(factor_bad)} 个。")
    L.append(f"- distinct run 字符串 {len(run_strings)}；解析成功 {len(run_key)}，"
             f"落入 (timepoint, fraction, rep) 设计格 {n_design_cells}/{N_DESIGN}；"
             f"解析失败 {len(parse_fail)}。")
    if parse_fail:
        L.append("  - 解析失败 run：")
        for rs, why in parse_fail:
            L.append(f"    - `{rs}`（{why}）")
    if multi_run_cells:
        L.append(f"  - 多个 run 字符串映射到同一设计格：{len(multi_run_cells)} 格")
        for k in sorted(multi_run_cells):
            L.append(f"    - {k}: " + ", ".join(f"`{x}`" for x in sorted(multi_run_cells[k])))
    L.append(f"- `genes` 含 EGFR 的行 {n_egfr_rows}；其中 run 解析失败 {egfr_parse_fail_rows}、"
             f"quantity 缺失 {egfr_missing_q_rows}、quantity ≤ 0 {egfr_nonpos_rows}；"
             f"同一 (protein_group, tp, fraction, rep) 多条同值记录 {dup_same} 格（计一次），"
             f"不同值冲突 {len(conflicts)} 格（不计入，见下）。")
    if conflicts:
        for key, lst in sorted(conflicts):
            L.append(f"  - 冲突 {key}: " + "; ".join(f"{g6(q)} (`{rs}`)" for q, rs in lst))
    L.append(f"- **EGFR 在 {n_egfr_runs}/{N_DESIGN} 个 run 测到**"
             f"（任一 EGFR protein_group 的 quantity 有值）。")
    L.append("")
    L.append("EGFR protein_group 表：")
    L.append("")
    if order:
        L += md_table(["gene", "protein_group", "genes_raw", "n_runs_detected", "is_primary"],
                      [[GENE, g, " / ".join(sorted(genes_raw[g])), f"{n_runs[g]}/{N_DESIGN}",
                        "yes" if g == primary else "no"] for g in order])
        if tie:
            L.append("")
            L.append(f"primary 并列（n_runs_detected 相同）：{', '.join(tie)}；按字符串排序取 `{primary}`。")
    else:
        L.append((SKIP_MSG if n_in == 0 else f"缺数据：EGFR 在 0/{N_DESIGN} 个 run 测到，跳过")
                 + "（无 EGFR protein_group）。")
    L.append("")

    L.append("## ② 定义")
    L.append("")
    L += [
        "- 匹配：`genes` 按 `;` 拆开（去首尾空白）后任一项等于 `EGFR` 即匹配，记录原始 `genes` 字符串。",
        "- 多个 `protein_group` 全部保留；primary = 测到 run 数最多的 group（测到 = run 解析成功、"
        "落入 120 设计格且 quantity 有值），并列取 `protein_group` 字符串排序第一并记录并列；"
        "其余为 secondary，列附表。",
        "- run 解析：07 `parse_design`（`code/Protein contour/Zhihan/Test/07_main_analysis.py:32-38`）"
        "三个正则 `_(2min|8min|20min|90min|CTRL)_`（忽略大小写，结果统一写成 CTRL/2min/8min/20min/90min）、"
        "`_(FR\\d)_`、`_(Rep\\d)`；fraction ∉ FR1..FR6 或 rep ∉ Rep1..Rep4 计为解析失败。",
        "- `quantity` 转 float；空串 / 非数 / 非有限值视为缺失。",
        "- 两版归一化：`raw` = quantity；`centered` = quantity × 2^a，a 取 "
        "`analysis/07_pilot2/01_run_factors.tsv` 中 layer=proteome、同 (timepoint, fraction, rep) 的 `a`"
        "（等价于 quantity / 2^(m_07 − anchor)，anchor 为常数）。",
        "- 总量：对每个 (protein_group, timepoint, rep, normalization)，"
        "`n_fractions_present` = FR1..FR6 中有值的个数；`total_quantity = Σ_present q_f`；不填补；"
        "`n_fractions_present = 0` 时 `total_quantity` 留空。",
        "- 比值：`ratio_r = total_90min,r / total_CTRL,r`，按 rep 编号配对，4 个 rep 各一行；"
        "任一侧 total 缺则 ratio 留空（CTRL total = 0 时也留空并注明）。",
    ]
    L.append("")

    have = bool(order)
    skip_line = SKIP_MSG if n_in == 0 else f"缺数据：EGFR 在 0/{N_DESIGN} 个 run 测到，跳过"

    def total_table(grp, norm):
        rows = []
        for tp in TPS:
            for rep in REPS:
                n, tot, missing = totals[(grp, norm, tp, rep)]
                rows.append([tp, rep, n, g6(tot),
                             "—" if not missing else ",".join(missing),
                             "是" if n < 6 else ""])
        return md_table(["timepoint", "rep", "n_fractions_present", "total_quantity",
                         "缺的 fraction", "有缺 fraction"], rows)

    def ratio_table(grp, norm):
        rows = []
        for rep in REPS:
            t90, tct, ratio, note = ratios[(grp, norm, rep)]
            rows.append([rep, g6(t90), totals[(grp, norm, "90min", rep)][0],
                         g6(tct), totals[(grp, norm, "CTRL", rep)][0], f4(ratio), note])
        return md_table(["rep", "total_90min", "n_fr_90min", "total_CTRL", "n_fr_CTRL",
                         "ratio (90min/CTRL)", "备注"], rows)

    L.append("## ③ EGFR 总量表")
    L.append("")
    L.append(f"EGFR 在 {n_egfr_runs}/{N_DESIGN} 个 run 测到。")
    L.append("")
    if not have:
        L.append(skip_line)
        L.append("")
    else:
        L.append(f"gene = EGFR；primary protein_group = `{primary}`"
                 f"（genes = `{' / '.join(sorted(genes_raw[primary]))}`，"
                 f"测到 {n_runs[primary]}/{N_DESIGN} run）。"
                 "total_quantity = 该 (timepoint, rep) 下有值 fraction 的 quantity 之和，不填补；"
                 "raw = quantity，centered = quantity × 2^a（a 来自 07 `01_run_factors.tsv`, proteome）。"
                 "`有缺 fraction` = 是 表示 n_fractions_present < 6。数值为 6 位有效数字，"
                 "完整精度见 `03_egfr_total.tsv`。")
        L.append("")
        for norm in NORMS:
            L.append(f"**{norm}**")
            L.append("")
            L += total_table(primary, norm)
            L.append("")
        miss_list = [(norm, tp, rep, m) for norm in NORMS for tp in TPS for rep in REPS
                     for (n, _, m) in [totals[(primary, norm, tp, rep)]] if n < 6]
        L.append(f"缺 fraction 清单（primary，n_fractions_present < 6）：{len(miss_list)} 行"
                 + ("。" if miss_list else "（无）。"))
        if miss_list:
            L.append("")
            L += md_table(["normalization", "timepoint", "rep", "缺的 fraction"],
                          [[a, b, c, ",".join(m)] for a, b, c, m in miss_list])
        L.append("")
        if centered_no_factor:
            L.append("centered 缺因子 a 的 (timepoint, fraction, rep)：" +
                     "; ".join(str(k) for k in sorted(centered_no_factor)) + "。")
            L.append("")
        secs = [g for g in order if g != primary]
        L.append(f"附表（secondary protein_group）：{len(secs)} 个" + ("。" if secs else "（无）。"))
        L.append("")
        for g in secs:
            for norm in NORMS:
                L.append(f"secondary `{g}`（测到 {n_runs[g]}/{N_DESIGN} run），**{norm}**")
                L.append("")
                L += total_table(g, norm)
                L.append("")

    L.append("## ④ EGFR 90 min / CTRL 比值表")
    L.append("")
    L.append(f"EGFR 在 {n_egfr_runs}/{N_DESIGN} 个 run 测到。")
    L.append("")
    if not have:
        L.append(skip_line)
        L.append("")
    else:
        L.append(f"gene = EGFR；primary protein_group = `{primary}`。"
                 "ratio_r = total_90min,r / total_CTRL,r，按 rep 编号配对；任一侧 total 缺则留空。"
                 "n_fr_* = 对应 total 的 n_fractions_present（< 6 即该侧有缺 fraction）。"
                 "total 为 6 位有效数字，ratio 为 4 位小数，完整精度见 `03_egfr_ratio.tsv`。")
        L.append("")
        for norm in NORMS:
            L.append(f"**{norm}**")
            L.append("")
            L += ratio_table(primary, norm)
            L.append("")
        secs = [g for g in order if g != primary]
        L.append(f"附表（secondary protein_group）：{len(secs)} 个" + ("。" if secs else "（无）。"))
        L.append("")
        for g in secs:
            for norm in NORMS:
                L.append(f"secondary `{g}`，**{norm}**")
                L.append("")
                L += ratio_table(g, norm)
                L.append("")

    L.append("## ⑤ 跳过项")
    L.append("")
    if n_in == 0:
        L.append(f"- {SKIP_MSG}。③ 总量表、④ 比值表均无数值；`03_egfr_total.tsv`、"
                 "`03_egfr_ratio.tsv` 只有表头。数据到位后直接复跑本脚本即可。")
    elif not have:
        L.append(f"- 缺数据：EGFR 在 0/{N_DESIGN} 个 run 测到，③、④ 跳过；两个 tsv 只有表头。")
    else:
        L.append("- 无整节跳过。")
    if conflicts:
        L.append(f"- 同一 (protein_group, tp, fraction, rep) 不同值冲突 {len(conflicts)} 格，未计入总量（见 ①）。")
    L.append("- 不复现 DAPAR 的 LOESS 归一化与 KNN/MEC/detQuant 填补（00_design.md §1），只给 raw 与 centered 两版。")
    L.append("- 脚本在 scratchpad 合成夹具上跑通。")
    L.append("")

    out_md.write_text("\n".join(L), encoding="utf-8")

    print(f"input rows: {n_in}; distinct runs: {len(run_strings)}; parsed: {len(run_key)}; "
          f"design cells: {n_design_cells}/{N_DESIGN}; parse failures: {len(parse_fail)}")
    print(f"EGFR rows: {n_egfr_rows}; EGFR detected in {n_egfr_runs}/{N_DESIGN} runs; "
          f"groups: {len(order)}; primary: {primary}")
    print(f"wrote {out_total} ({len(total_rows)} rows), {out_ratio} ({len(ratio_rows)} rows), {out_md}")


if __name__ == "__main__":
    main()
