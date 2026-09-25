#!/usr/bin/env python3
"""04_qc_summary.py — 汇总 240 个 run 的质量指标, 标红异常项。

用法: python3 04_qc_summary.py <phospho|proteome> [结果根目录]
默认根目录: $VSC_DATA/rerun_2026-08
输出: qc_<line>.tsv (当前目录) + 屏幕摘要

每 run 指标:
  run, timepoint, fraction, rep, precursors, proteins,
  phospho 线额外: phos_precursors, phos_pct, loc075_pass (PTM.Site.Confidence >= 0.75)
红旗规则:
  - stats/parquet 缺失
  - precursors < 全体中位数的 50%
  - phospho 线 phos_pct < 50%
  - 5 时间点 x 6 组分 x 4 重复 的网格有缺口
"""
import sys, os, re, csv, glob, statistics

LOC_CUTOFF = 0.75

def parse_design(run):
    tp = re.search(r"_(2min|8min|20min|90min|CTRL)_", run, re.I)
    fr = re.search(r"_(FR\d)_", run)
    rep = re.search(r"_(Rep\d)", run)
    return (tp.group(1) if tp else "?", fr.group(1) if fr else "?", rep.group(1) if rep else "?")

def read_stats(d):
    """DIA-NN report.stats.tsv: 取 Precursors.Identified / Proteins.Identified"""
    p = os.path.join(d, "report.stats.tsv")
    if not os.path.isfile(p):
        return None, None
    with open(p) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        return None, None
    r = rows[0]
    def num(key):
        for k in r:
            if k.strip() == key:
                try:
                    return int(float(r[k]))
                except (ValueError, TypeError):
                    return None
        return None
    return num("Precursors.Identified"), num("Proteins.Identified")

def scan_phospho(d):
    """扫 report.parquet: 总行/磷酸化行/定位>=0.75 行"""
    p = os.path.join(d, "report.parquet")
    if not os.path.isfile(p):
        return None
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(p)
    cols = pf.schema_arrow.names
    seq = next((c for c in ("Modified.Sequence", "ModifiedPeptide") if c in cols), None)
    conf = next((c for c in cols if "Site.Confidence" in c), None)
    if not seq:
        return None
    total = phos = loc = 0
    use = [c for c in (seq, conf) if c]
    for b in pf.iter_batches(columns=use, batch_size=100_000):
        total += b.num_rows
        s = b.column(seq).to_pylist()
        cf = b.column(conf).to_pylist() if conf else [None] * b.num_rows
        for i in range(b.num_rows):
            if s[i] and "UniMod:21" in s[i]:
                phos += 1
                if cf[i] is not None and cf[i] >= LOC_CUTOFF:
                    loc += 1
    return total, phos, loc, (conf is not None)

def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ("phospho", "proteome"):
        print(__doc__); sys.exit(2)
    line = sys.argv[1]
    base = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.environ["VSC_DATA"], "rerun_2026-08")
    dirs = sorted(glob.glob(os.path.join(base, line, "*/")))
    print(f"{line}: 发现 {len(dirs)} 个 run 目录")

    rows, flags = [], []
    have_conf_col = True
    for d in dirs:
        run = os.path.basename(d.rstrip("/"))
        tp, fr, rep = parse_design(run)
        prec, prot = read_stats(d)
        row = {"run": run, "timepoint": tp, "fraction": fr, "rep": rep,
               "precursors": prec, "proteins": prot}
        if prec is None:
            flags.append(f"[缺文件] {run}: report.stats.tsv 缺失或空")
        if line == "phospho":
            r = scan_phospho(d)
            if r is None:
                flags.append(f"[缺文件] {run}: report.parquet 缺失或无序列列")
                row.update({"phos_precursors": None, "phos_pct": None, "loc075_pass": None})
            else:
                total, phos, loc, has_conf = r
                have_conf_col &= has_conf
                pct = round(phos / total * 100, 1) if total else 0.0
                row.update({"phos_precursors": phos, "phos_pct": pct, "loc075_pass": loc})
                if pct < 50:
                    flags.append(f"[低磷酸化占比] {run}: {pct}%")
        rows.append(row)

    # 离群: precursors < 中位数一半
    vals = [r["precursors"] for r in rows if r["precursors"]]
    if vals:
        med = statistics.median(vals)
        for r in rows:
            if r["precursors"] is not None and r["precursors"] < med * 0.5:
                flags.append(f"[低鉴定数] {r['run']}: {r['precursors']} (中位数 {int(med)})")

    # 设计网格完整性
    grid = {(r["timepoint"], r["fraction"], r["rep"]) for r in rows}
    expect = {(t, f"FR{i}", f"Rep{j}") for t in ("2min", "8min", "20min", "90min", "CTRL")
              for i in range(1, 7) for j in range(1, 5)}
    missing = expect - grid
    for m in sorted(missing):
        flags.append(f"[网格缺口] 缺 {m[0]} {m[1]} {m[2]}")

    out = f"qc_{line}.tsv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    print(f"\n== 摘要 ({line}) ==")
    if vals:
        print(f"precursors: 中位 {int(statistics.median(vals))} | 最小 {min(vals)} | 最大 {max(vals)}")
    if line == "phospho":
        pcts = [r["phos_pct"] for r in rows if r.get("phos_pct") is not None]
        locs = [r["loc075_pass"] for r in rows if r.get("loc075_pass") is not None]
        if pcts:
            print(f"phos_pct: 中位 {statistics.median(pcts)}% | 最小 {min(pcts)}%")
        if locs:
            print(f"loc>= {LOC_CUTOFF} 通过数: 中位 {int(statistics.median(locs))}")
        if not have_conf_col:
            print("WARN: 部分 run 无 PTM.Site.Confidence 列 — 0.75 过滤依据需人工核实")
    print(f"\n红旗 {len(flags)} 条:")
    for fl in flags[:40]:
        print("  " + fl)
    if len(flags) > 40:
        print(f"  ... 共 {len(flags)} 条, 其余见输出文件对照")
    print(f"\n明细已写入 {out}")
    sys.exit(1 if flags else 0)

if __name__ == "__main__":
    main()
