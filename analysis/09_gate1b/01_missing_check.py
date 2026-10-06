#!/usr/bin/env python3
"""01_missing_check.py — 09_gate1b 任务 1：缺失核查（Q2b）。

按 analysis/09_gate1b/00_design.md §0 与「任务 1：缺失核查（Q2b）」执行。
输入（只读）:
  data/pilot/gate1_missing_check.tsv        列 run, gene, n_precursors, min_PG_qvalue, min_precursor_qvalue（480 行）
  analysis/gate1_proteins.tsv               列 run, genes, protein_group, quantity（344 行）
  data/pilot/proteome_pg_matrix_120.tsv     交叉核对用（Protein.Group, Genes + 120 run 列）
  code/Protein contour/Zhihan/Test/07_main_analysis.py   只读取文本，核对下面照抄的正则与 :32-38 一致
输出（analysis/09_gate1b/）:
  01_missing_classes.tsv, 01_missing_classes.md

定义（照方案「任务 1」）:
  对每个 (gene ∈ {GRB2, SHC1, CBL, EGFR}, run ∈ 120)：
    A = analysis/gate1_proteins.tsv 中该 (gene, run) 有 quantity（非空串）；
    B = 非 A 且 n_precursors > 0；
    C = 非 A 且 n_precursors == 0。
  另报：A 类中 n_precursors == 0 的行数；B 类按 min_PG_qvalue <= 0.01 / > 0.01 拆分。
  读表：sep='\\t', dtype=str, keep_default_na=False，空串 = 缺失。
  run 解析：07 parse_design 三个正则 _(2min|8min|20min|90min|CTRL)_ (re.I)、_(FR\\d)_、_(Rep\\d)。
  gene 匹配（gate1_proteins）：genes 按 ';' 拆开，任一 token 等于目标基因即匹配。
  输出行序：gene（GRB2, SHC1, CBL, EGFR）→ timepoint（CTRL, 2min, 8min, 20min, 90min）→ fraction → rep。
  pg_quantity：gate1_proteins 的 quantity 原字符串（非 A 为空串）；q 值列照抄原字符串。

用法: python3 analysis/09_gate1b/01_missing_check.py   （任意目录可复跑，输出确定）
"""
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
F_MISS = ROOT / "data/pilot/gate1_missing_check.tsv"
F_PROT = ROOT / "analysis/gate1_proteins.tsv"
F_PROT_OLD = ROOT / "data/pilot/gate1_proteins.tsv"
F_PGM = ROOT / "data/pilot/proteome_pg_matrix_120.tsv"
F_07 = ROOT / "code/Protein contour/Zhihan/Test/07_main_analysis.py"
OUTDIR = ROOT / "analysis/09_gate1b"
OUT_TSV = OUTDIR / "01_missing_classes.tsv"
OUT_MD = OUTDIR / "01_missing_classes.md"

GENES = ["GRB2", "SHC1", "CBL", "EGFR"]
TPS = ["CTRL", "2min", "8min", "20min", "90min"]
FRS = [f"FR{i}" for i in range(1, 7)]
REPS = [f"Rep{i}" for i in range(1, 5)]
CLASSES = ["A", "B", "C"]
Q_TH = 0.01
OUT_COLS = ["gene", "run", "timepoint", "fraction", "rep", "class", "n_precursors",
            "min_PG_qvalue", "min_precursor_qvalue", "pg_quantity"]

# advisor 参考值（00_design.md「任务 1」验收），仅用于自检比对，不参与分类
REF_ABC = {"GRB2": (90, 3, 27), "SHC1": (115, 1, 4), "CBL": (59, 12, 49), "EGFR": (80, 7, 33)}
REF_B_QLE = {"GRB2": 0, "SHC1": 0, "CBL": 1, "EGFR": 0}


# ---- 照抄 code/Protein contour/Zhihan/Test/07_main_analysis.py parse_design (:32-38) ----
def parse_design(run):
    tp = re.search(r"_(2min|8min|20min|90min|CTRL)_", run, re.I)
    fr = re.search(r"_(FR\d)_", run)
    rep = re.search(r"_(Rep\d)", run)
    if not (tp and fr and rep):
        return None
    return (tp.group(1), fr.group(1), rep.group(1))
# ---- 照抄结束 ----

REGEX_LINES = [
    'tp = re.search(r"_(2min|8min|20min|90min|CTRL)_", run, re.I)',
    'fr = re.search(r"_(FR\\d)_", run)',
    'rep = re.search(r"_(Rep\\d)", run)',
]


def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x).replace("|", "\\|") for x in r) + " |")
    return out


def main():
    checks = []  # (name, ok, detail)

    def check(name, ok, detail=""):
        checks.append((name, bool(ok), detail))

    # ---- 07 正则核对（只读文本）----
    src_lines = F_07.read_text(encoding="utf-8").splitlines()
    seg = [l.strip() for l in src_lines[31:38]]  # 1-based :32-38
    check("07_main_analysis.py:32-38 为 parse_design，三个正则与本脚本照抄一致",
          seg[0].startswith("def parse_design(run):") and all(rl in seg for rl in REGEX_LINES),
          "; ".join(seg[1:4]))

    # ---- 读入 ----
    miss = read_tsv(F_MISS)
    prot = read_tsv(F_PROT)
    pgm = read_tsv(F_PGM)
    prot_old = read_tsv(F_PROT_OLD) if F_PROT_OLD.exists() else None

    check("gate1_missing_check.tsv 列名", list(miss.columns) ==
          ["run", "gene", "n_precursors", "min_PG_qvalue", "min_precursor_qvalue"], str(list(miss.columns)))
    check("gate1_proteins.tsv 列名", list(prot.columns) == ["run", "genes", "protein_group", "quantity"],
          str(list(prot.columns)))
    check("gate1_missing_check.tsv 480 行", len(miss) == 480, str(len(miss)))
    check("analysis/gate1_proteins.tsv 344 行", len(prot) == 344, str(len(prot)))

    # ---- run 解析 ----
    runs = sorted(set(miss["run"]))
    design = {}
    bad = []
    for r in runs:
        p = parse_design(r)
        if p is None:
            bad.append(r)
        else:
            design[r] = p
    check("missing_check 120 个 run 且全部可解析", len(runs) == 120 and not bad, f"n_run={len(runs)}, 解析失败={len(bad)}")
    cells = set(design.values())
    check("run 解析为 5 tp × 6 FR × 4 rep 且一一对应",
          cells == {(t, f, r) for t in TPS for f in FRS for r in REPS} and len(cells) == len(runs),
          f"不同 (tp,FR,rep) = {len(cells)}")

    # ---- missing_check 结构 ----
    check("missing_check 的 gene 集合 = 四蛋白，每基因 120 行",
          set(miss["gene"]) == set(GENES) and all((miss["gene"] == g).sum() == 120 for g in GENES),
          str(miss["gene"].value_counts().to_dict()))
    check("missing_check (run, gene) 无重复", not miss.duplicated(["run", "gene"]).any(),
          str(int(miss.duplicated(["run", "gene"]).sum())))
    npre_ok = miss["n_precursors"].str.fullmatch(r"\d+").all()
    check("n_precursors 全为非负整数", npre_ok, "")
    miss["n_pre_int"] = miss["n_precursors"].astype(int)
    z = miss["n_pre_int"] == 0
    q_empty_z = ((miss.loc[z, "min_PG_qvalue"] == "") & (miss.loc[z, "min_precursor_qvalue"] == "")).all()
    q_full_nz = ((miss.loc[~z, "min_PG_qvalue"] != "") & (miss.loc[~z, "min_precursor_qvalue"] != "")).all()
    check("n_precursors == 0 的行两列 q 值均为空串；> 0 的行均非空", q_empty_z and q_full_nz,
          f"n_precursors==0 行数 {int(z.sum())}")

    # ---- gate1_proteins → (gene, run) 的 quantity ----
    qty = {}
    multi = []
    unmatched = 0
    for _, r in prot.iterrows():
        toks = [t for t in r["genes"].split(";") if t != ""]
        hit = [g for g in GENES if g in toks]
        if not hit:
            unmatched += 1
            continue
        for g in hit:
            key = (g, r["run"])
            if key in qty:
                multi.append(key)
            qty[key] = (r["quantity"], r["protein_group"])
    check("gate1_proteins 每行都匹配到四蛋白之一，(gene, run) 无重复", unmatched == 0 and not multi,
          f"未匹配 {unmatched} 行，重复 {len(multi)} 个")
    check("gate1_proteins 的 run 都在 missing_check 的 120 run 内",
          set(prot["run"]) <= set(runs), f"多出 {len(set(prot['run']) - set(runs))}")
    check("gate1_proteins quantity 全部非空且可转为正数",
          all(v[0] != "" for v in qty.values()) and all(float(v[0]) > 0 for v in qty.values()), "")
    pg_by_gene = {g: sorted({v[1] for k, v in qty.items() if k[0] == g}) for g in GENES}

    # ---- 分类 ----
    rows = []
    for _, r in miss.iterrows():
        g, run = r["gene"], r["run"]
        tp, fr, rep = design[run]
        q = qty.get((g, run), ("", ""))[0]
        n = int(r["n_pre_int"])
        if q != "":
            cls = "A"
        elif n > 0:
            cls = "B"
        else:
            cls = "C"
        rows.append({"gene": g, "run": run, "timepoint": tp, "fraction": fr, "rep": rep, "class": cls,
                     "n_precursors": r["n_precursors"], "min_PG_qvalue": r["min_PG_qvalue"],
                     "min_precursor_qvalue": r["min_precursor_qvalue"], "pg_quantity": q})
    out = pd.DataFrame(rows, columns=OUT_COLS)
    out["_g"] = out["gene"].map({g: i for i, g in enumerate(GENES)})
    out["_t"] = out["timepoint"].map({t: i for i, t in enumerate(TPS)})
    out = out.sort_values(["_g", "_t", "fraction", "rep"]).drop(columns=["_g", "_t"]).reset_index(drop=True)
    check("所有 gate1_proteins 的 (gene, run) 都落入 A", (out["class"] == "A").sum() == len(qty),
          f"A={int((out['class'] == 'A').sum())}, gate1_proteins={len(qty)}")
    check("输出 480 行、列顺序正确", len(out) == 480 and list(out.columns) == OUT_COLS, str(len(out)))

    # ---- 计数 ----
    abc = {g: {c: int(((out["gene"] == g) & (out["class"] == c)).sum()) for c in CLASSES} for g in GENES}
    for g in GENES:
        tot = sum(abc[g].values())
        check(f"{g} A+B+C = 120", tot == 120, str(tot))
        got = (abc[g]["A"], abc[g]["B"], abc[g]["C"])
        check(f"{g} A/B/C 与 advisor 参考一致", got == REF_ABC[g], f"本轮 {got} vs 参考 {REF_ABC[g]}")

    a_zero = out[(out["class"] == "A") & (out["n_precursors"].astype(int) == 0)]
    check("A 类中 n_precursors == 0 的行数 = 0", len(a_zero) == 0, str(len(a_zero)))

    bdf = out[out["class"] == "B"].copy()
    bdf["_pgq"] = bdf["min_PG_qvalue"].astype(float)
    b_le = {g: int(((bdf["gene"] == g) & (bdf["_pgq"] <= Q_TH)).sum()) for g in GENES}
    b_gt = {g: int(((bdf["gene"] == g) & (bdf["_pgq"] > Q_TH)).sum()) for g in GENES}
    check("B 类 min_PG_qvalue ≤ 0.01 计数与 advisor 参考一致（CBL 1，其余 0）", b_le == REF_B_QLE, str(b_le))
    b_eq = int((bdf["_pgq"] == Q_TH).sum())

    frtab = {}
    for g in GENES:
        for f in FRS:
            sub = out[(out["gene"] == g) & (out["fraction"] == f)]
            frtab[(g, f)] = {c: int((sub["class"] == c).sum()) for c in CLASSES}
            frtab[(g, f)]["tot"] = len(sub)
    check("每 (蛋白, fraction) 行合计 = 20", all(v["tot"] == 20 for v in frtab.values()),
          str(sorted({v["tot"] for v in frtab.values()})))

    # ---- 与 pg_matrix 交叉核对 ----
    pg_info = {}
    pgm_runs = [c for c in pgm.columns if c not in ("Protein.Group", "Genes")]
    check("pg_matrix 的 120 个 run 列与 missing_check 的 run 集合一致", set(pgm_runs) == set(runs),
          f"pg_matrix run 列 {len(pgm_runs)}")
    for g in GENES:
        hit = pgm[pgm["Genes"].apply(lambda s: g in [t for t in s.split(";") if t != ""])]
        a_runs = set(out.loc[(out["gene"] == g) & (out["class"] == "A"), "run"])
        if len(hit) != 1:
            pg_info[g] = (len(hit), "", None, False, None)
            check(f"pg_matrix 中 {g} 恰好 1 行", False, str(len(hit)))
            continue
        row = hit.iloc[0]
        nonempty = {c for c in pgm_runs if row[c] != ""}
        same = nonempty == a_runs
        # 数值一致性：gate1_proteins 的 quantity 与 pg_matrix 值的最大相对差
        rel = 0.0
        for run in a_runs & nonempty:
            a = float(qty[(g, run)][0])
            b = float(row[run])
            rel = max(rel, abs(a - b) / abs(b))
        pg_info[g] = (1, row["Protein.Group"], len(nonempty), same, rel)
        check(f"pg_matrix {g} 行（{row['Protein.Group']}）非空 run 集合 = A 类 run 集合", same,
              f"pg_matrix 非空 {len(nonempty)}，A {len(a_runs)}，对称差 {len(nonempty ^ a_runs)}")
        check(f"pg_matrix {g} 行与 gate1_proteins 的 protein_group 一致",
              pg_by_gene[g] == [row["Protein.Group"]], f"{pg_by_gene[g]} vs {row['Protein.Group']}")

    old_same = None
    if prot_old is not None:
        old_same = F_PROT_OLD.read_bytes() == F_PROT.read_bytes()

    # ---- 写 tsv ----
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    # ---- 写 md ----
    L = []
    L.append("# 09_gate1b 任务 1：缺失核查（Q2b）")
    L.append("")
    L.append("脚本：`analysis/09_gate1b/01_missing_check.py`（`python3 analysis/09_gate1b/01_missing_check.py` 可复跑，输出确定）。"
             f"pandas {pd.__version__}。逐行结果：`analysis/09_gate1b/01_missing_classes.tsv`（{len(out)} 行）。")
    L.append("")
    L.append("## ① 输入")
    L.append("")
    L += md_table(["文件", "用途", "行数 / 说明"], [
        ["`data/pilot/gate1_missing_check.tsv`", "n_precursors、min_PG_qvalue、min_precursor_qvalue",
         f"{len(miss)} 行 = 4 蛋白 × {len(runs)} run；n_precursors == 0 的行 {int(z.sum())}（"
         + " / ".join(f"{g} {int(((miss['gene'] == g) & z).sum())}" for g in GENES) + "），这些行两列 q 值均为空串"],
        ["`analysis/gate1_proteins.tsv`", "A 类判定（有 quantity）",
         f"{len(prot)} 行；" + " / ".join(f"{g} {abc[g]['A']}（{', '.join(pg_by_gene[g])}）" for g in GENES)],
        ["`data/pilot/proteome_pg_matrix_120.tsv`", "交叉核对",
         f"{len(pgm)} 行 × {len(pgm.columns)} 列；四基因各 1 行"],
        ["`code/Protein contour/Zhihan/Test/07_main_analysis.py:32-38`", "run 名解析",
         "`parse_design` 三个正则 `_(2min|8min|20min|90min|CTRL)_`（re.I）、`_(FR\\d)_`、`_(Rep\\d)`；脚本照抄，并核对源文件该段文本一致"],
    ])
    L.append("")
    L.append(f"- run 解析：{len(runs)} 个 run 全部可解析，解析失败 {len(bad)}；(timepoint, fraction, rep) 共 {len(cells)} 种 = 5 × 6 × 4，与 run 一一对应。")
    L.append("- 与 pg_matrix 交叉核对（四基因行的非空 run 集合 vs 本表 A 类 run 集合；数值为 gate1_proteins quantity 与 pg_matrix 值的最大相对差）：")
    L.append("")
    L += md_table(["gene", "pg_matrix Protein.Group", "pg_matrix 非空 run 数", "A 类 run 数", "集合一致", "最大相对差"],
                  [[g, pg_info[g][1], pg_info[g][2], abc[g]["A"], "是" if pg_info[g][3] else "否",
                    f"{pg_info[g][4]:.2e}" if pg_info[g][4] is not None else ""] for g in GENES])
    L.append("")
    if old_same is not None:
        L.append(f"- `data/pilot/gate1_proteins.tsv`：现为 {len(prot_old)} 行，与 `analysis/gate1_proteins.tsv` "
                 + ("逐字节相同" if old_same else "内容不同")
                 + "；00_design §0 记为 0 行，与现状不符。本轮按任务书只用 `analysis/gate1_proteins.tsv`。")
        L.append("")
    L.append("## ② 定义")
    L.append("")
    L.append("对每个 (gene ∈ {GRB2, SHC1, CBL, EGFR}, run ∈ 120)：")
    L.append("")
    L.append("- **A**：`analysis/gate1_proteins.tsv` 中该 (gene, run) 的 `quantity` 非空（gene 匹配：`genes` 按 `;` 拆开，任一 token 等于目标基因）。")
    L.append("- **B**：非 A 且 `gate1_missing_check.tsv` 的 `n_precursors > 0`。")
    L.append("- **C**：非 A 且 `n_precursors == 0`。")
    L.append(f"- B 类按 `min_PG_qvalue` 拆分：≤ {Q_TH} / > {Q_TH}（`min_PG_qvalue` 转 float 比较；恰等于 {Q_TH} 的 B 行 {b_eq} 个）。")
    L.append("- 读表 `sep='\\t', dtype=str, keep_default_na=False`，空串 = 缺失；tsv 中 `n_precursors`、q 值、`pg_quantity` 均为原字符串（非 A 的 `pg_quantity` 为空串）。")
    L.append("- tsv 行序：gene（GRB2, SHC1, CBL, EGFR）→ timepoint（CTRL, 2min, 8min, 20min, 90min）→ fraction → rep。")
    L.append("")
    L.append("## ③ 每蛋白 A/B/C 计数")
    L.append("")
    L += md_table(["gene", "A", "B", "C", "合计", "advisor 参考 A/B/C", "一致"],
                  [[g, abc[g]["A"], abc[g]["B"], abc[g]["C"], sum(abc[g].values()),
                    "/".join(str(x) for x in REF_ABC[g]),
                    "是" if (abc[g]["A"], abc[g]["B"], abc[g]["C"]) == REF_ABC[g] else "否"] for g in GENES]
                  + [["合计", sum(abc[g]["A"] for g in GENES), sum(abc[g]["B"] for g in GENES),
                      sum(abc[g]["C"] for g in GENES), len(out), "", ""]])
    L.append("")
    L.append(f"- A 类中 `n_precursors == 0` 的行数：{len(a_zero)}（A 类 `n_precursors` 最小值 "
             f"{int(out.loc[out['class'] == 'A', 'n_precursors'].astype(int).min())}）。")
    L.append("")
    L.append("## ④ 每蛋白 × fraction 的 A/B/C")
    L.append("")
    for g in GENES:
        L.append(f"**{g}**")
        L.append("")
        L += md_table(["fraction", "A", "B", "C", "合计"],
                      [[f, frtab[(g, f)]["A"], frtab[(g, f)]["B"], frtab[(g, f)]["C"], frtab[(g, f)]["tot"]] for f in FRS]
                      + [["合计", abc[g]["A"], abc[g]["B"], abc[g]["C"], sum(abc[g].values())]])
        L.append("")
    L.append("## ⑤ B 类清单")
    L.append("")
    L.append(f"B 类共 {len(bdf)} 行。按 `min_PG_qvalue` 拆分：")
    L.append("")
    L += md_table(["gene", "B", f"min_PG_qvalue ≤ {Q_TH}", f"min_PG_qvalue > {Q_TH}"],
                  [[g, abc[g]["B"], b_le[g], b_gt[g]] for g in GENES]
                  + [["合计", len(bdf), sum(b_le.values()), sum(b_gt.values())]])
    L.append("")
    L.append("清单（行序同 tsv）：")
    L.append("")
    L += md_table(["gene", "run", "n_precursors", "min_PG_qvalue", "min_precursor_qvalue"],
                  [[r["gene"], r["run"], r["n_precursors"], r["min_PG_qvalue"], r["min_precursor_qvalue"]]
                   for _, r in bdf.iterrows()])
    L.append("")
    L.append("## ⑥ 跳过项")
    L.append("")
    L.append("- 无跳过项：三个输入文件均存在且结构符合方案。")
    L.append("- 自检（脚本内全部执行）：")
    L.append("")
    L += md_table(["检查", "结果", "细节"], [[n, "通过" if ok else "**未通过**", d] for n, ok, d in checks])
    L.append("")
    OUT_MD.write_text("\n".join(L) + "\n", encoding="utf-8")

    # ---- 控制台 ----
    n_fail = sum(1 for _, ok, _ in checks if not ok)
    for n, ok, d in checks:
        print(("PASS " if ok else "FAIL ") + n + (f"  [{d}]" if d else ""))
    print(f"wrote {OUT_TSV} ({len(out)} rows), {OUT_MD}")
    print(f"checks: {len(checks) - n_fail}/{len(checks)} passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
