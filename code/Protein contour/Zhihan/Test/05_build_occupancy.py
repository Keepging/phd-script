#!/usr/bin/env python3
"""05_build_occupancy.py — 位点级 occupancy 输入表 (方案 A)

逻辑:
  1. 读 FASTA (UniProt accession -> 序列)
  2. 扫 120 个 phospho run 的 report.parquet:
     Q.Value<=0.01 且含 UniMod:21 且 PTM.Site.Confidence>=0.75 的前体
     -> 解析每个磷酸位点 -> 映射到蛋白绝对残基位置
     -> 分子: 该位点在该 run 的 Precursor.Quantity 之和
  3. 扫 120 个 proteome run: 覆盖同一残基的所有肽段 (任意修饰形式)
     -> 分母: Precursor.Quantity 之和
  4. 按 (timepoint, fraction, rep) 严格配对, 输出长表

已知方法学取舍 (记入论文方法节):
  - 多位点前体的强度全额记给它的每个位点 (不做拆分)
  - 共享肽段的位点归属取 Protein.Ids 中第一个能定位的 accession, 其余记入 ambiguity 列
  - occupancy_index = phospho 富集 run 强度 / proteome run 强度, 是相对指数而非绝对化学计量

用法: python3 05_build_occupancy.py [结果根目录] [FASTA路径]
默认: $VSC_DATA/rerun_2026-08 和 $VSC_SCRATCH/rerun/fasta/ 下最新的 human_*_UP*.fasta
输出: <根目录>/occupancy/{sites.tsv, occupancy_long.tsv, summary.txt}
"""
import sys, os, re, csv, glob
from collections import defaultdict

LOC_CUTOFF = 0.75
QVAL_CUTOFF = 0.01

def log(msg):
    print(msg, flush=True)

def parse_design(run):
    tp = re.search(r"_(2min|8min|20min|90min|CTRL)_", run, re.I)
    fr = re.search(r"_(FR\d)_", run)
    rep = re.search(r"_(Rep\d)", run)
    if not (tp and fr and rep):
        return None
    return (tp.group(1), fr.group(1), rep.group(1))

def load_fasta(path):
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

def phospho_sites_in(modseq):
    """返回 [(肽内位置1-based, 残基)]; 解析 (UniMod:21) 标注"""
    sites = []
    pos = 0
    last = ""
    i = 0
    while i < len(modseq):
        c = modseq[i]
        if c == "(":
            j = modseq.index(")", i)
            if modseq[i + 1:j] == "UniMod:21" and pos > 0:
                sites.append((pos, last))
            i = j + 1
        else:
            if c.isalpha():
                pos += 1
                last = c
            i += 1
    return sites

def pick_col(cols, exact, fuzzy):
    for c in exact:
        if c in cols:
            return c
    for c in cols:
        if fuzzy(c):
            return c
    return None

def open_report(d):
    import pyarrow.parquet as pq
    p = os.path.join(d, "report.parquet")
    if not os.path.isfile(p):
        return None, None
    pf = pq.ParquetFile(p)
    return pf, pf.schema_arrow.names

def resolve_cols(cols, need_conf):
    c = {}
    c["mod"] = pick_col(cols, ["Modified.Sequence"], lambda x: "Modified" in x and "Sequence" in x)
    c["strip"] = pick_col(cols, ["Stripped.Sequence"], lambda x: "Stripped" in x)
    c["prot"] = pick_col(cols, ["Protein.Ids"], lambda x: x == "Protein.Group")
    c["qty"] = pick_col(cols, ["Precursor.Quantity"], lambda x: False)
    c["qval"] = pick_col(cols, ["Q.Value"], lambda x: x.endswith("Q.Value") and "PG" not in x and "PTM" not in x and "Lib" not in x and "Global" not in x)
    c["conf"] = pick_col(cols, ["PTM.Site.Confidence"], lambda x: "Site.Confidence" in x) if need_conf else None
    return c

def main():
    base = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ["VSC_DATA"], "rerun_2026-08")
    if len(sys.argv) > 2:
        fasta_path = sys.argv[2]
    else:
        cands = sorted(glob.glob(os.path.join(os.environ["VSC_SCRATCH"], "rerun/fasta/human_*_UP*.fasta")))
        fasta_path = cands[-1] if cands else ""
    if not os.path.isfile(fasta_path):
        log(f"FAIL: FASTA 不存在: {fasta_path}"); sys.exit(1)
    outdir = os.path.join(base, "occupancy")
    os.makedirs(outdir, exist_ok=True)

    seqs = load_fasta(fasta_path)
    log(f"FASTA: {fasta_path} ({len(seqs)} 条)")

    phos_dirs = sorted(glob.glob(os.path.join(base, "phospho", "*/")))
    prot_dirs = sorted(glob.glob(os.path.join(base, "proteome", "*/")))
    log(f"phospho run: {len(phos_dirs)} | proteome run: {len(prot_dirs)}")
    if not phos_dirs or not prot_dirs:
        log("FAIL: run 目录缺失"); sys.exit(1)

    # ---------- 第一遍: phospho -> 位点 + 分子 ----------
    find_cache = {}          # (acc, pep) -> start(0-based) or -1
    site_info = {}           # site_id -> (acc, pos, residue)
    site_ambig = defaultdict(set)
    numer = defaultdict(float)   # (site_id, designkey) -> intensity
    skipped_unmapped = 0
    checked_cols = False

    for d in phos_dirs:
        run = os.path.basename(d.rstrip("/"))
        key = parse_design(run)
        if key is None:
            log(f"WARN: 无法解析设计标签, 跳过 {run}"); continue
        pf, cols = open_report(d)
        if pf is None:
            log(f"WARN: 缺 report.parquet, 跳过 {run}"); continue
        c = resolve_cols(cols, need_conf=True)
        missing = [k for k in ("mod", "strip", "prot", "qty", "conf") if not c[k]]
        if missing:
            log(f"FAIL: {run} 缺列 {missing}; 实际列: {cols}"); sys.exit(1)
        if not checked_cols:
            log(f"phospho 列映射: {c}")
            checked_cols = True
        use = [c[k] for k in ("mod", "strip", "prot", "qty", "conf") if c[k]]
        if c["qval"]:
            use.append(c["qval"])
        for b in pf.iter_batches(columns=use, batch_size=100_000):
            cd = {name: b.column(name).to_pylist() for name in use}
            for i in range(b.num_rows):
                if c["qval"] and (cd[c["qval"]][i] or 1) > QVAL_CUTOFF:
                    continue
                mod = cd[c["mod"]][i] or ""
                if "UniMod:21" not in mod:
                    continue
                conf = cd[c["conf"]][i]
                if conf is None or conf < LOC_CUTOFF:
                    continue
                qty = cd[c["qty"]][i]
                if not qty or qty <= 0:
                    continue
                pep = cd[c["strip"]][i] or ""
                accs = [a for a in re.split(r"[;,]", cd[c["prot"]][i] or "") if a]
                psites = phospho_sites_in(mod)
                if not pep or not accs or not psites:
                    continue
                placed = None
                for acc in accs:
                    seq = seqs.get(acc)
                    if seq is None:
                        continue
                    ck = (acc, pep)
                    if ck not in find_cache:
                        find_cache[ck] = seq.find(pep)
                    start = find_cache[ck]
                    if start >= 0:
                        placed = (acc, start)
                        break
                if placed is None:
                    skipped_unmapped += 1
                    continue
                acc, start = placed
                for ppos, res in psites:
                    abspos = start + ppos      # 1-based 蛋白位置
                    sid = f"{acc}_{res}{abspos}"
                    site_info[sid] = (acc, abspos, res)
                    if len(accs) > 1:
                        site_ambig[sid].update(a for a in accs if a != acc)
                    numer[(sid, key)] += qty
        log(f"  [分子] {run} 完成 (累计位点 {len(site_info)})")

    log(f"位点总数: {len(site_info)} | 未能映射到 FASTA 的磷酸化前体行: {skipped_unmapped}")
    if not site_info:
        log("FAIL: 零位点, 检查列映射与 FASTA"); sys.exit(1)

    # 蛋白 -> 位点索引 (供分母)
    sites_by_acc = defaultdict(list)
    for sid, (acc, pos, _res) in site_info.items():
        sites_by_acc[acc].append((pos, sid))

    # ---------- 第二遍: proteome -> 分母 ----------
    denom = defaultdict(float)
    checked_cols = False
    for d in prot_dirs:
        run = os.path.basename(d.rstrip("/"))
        key = parse_design(run)
        if key is None:
            log(f"WARN: 无法解析设计标签, 跳过 {run}"); continue
        pf, cols = open_report(d)
        if pf is None:
            log(f"WARN: 缺 report.parquet, 跳过 {run}"); continue
        c = resolve_cols(cols, need_conf=False)
        missing = [k for k in ("strip", "prot", "qty") if not c[k]]
        if missing:
            log(f"FAIL: {run} 缺列 {missing}; 实际列: {cols}"); sys.exit(1)
        if not checked_cols:
            log(f"proteome 列映射: {c}")
            checked_cols = True
        use = [c[k] for k in ("strip", "prot", "qty") if c[k]]
        if c["qval"]:
            use.append(c["qval"])
        for b in pf.iter_batches(columns=use, batch_size=100_000):
            cd = {name: b.column(name).to_pylist() for name in use}
            for i in range(b.num_rows):
                if c["qval"] and (cd[c["qval"]][i] or 1) > QVAL_CUTOFF:
                    continue
                qty = cd[c["qty"]][i]
                if not qty or qty <= 0:
                    continue
                pep = cd[c["strip"]][i] or ""
                accs = [a for a in re.split(r"[;,]", cd[c["prot"]][i] or "") if a]
                for acc in accs:
                    if acc not in sites_by_acc:
                        continue
                    seq = seqs.get(acc)
                    if seq is None:
                        continue
                    ck = (acc, pep)
                    if ck not in find_cache:
                        find_cache[ck] = seq.find(pep)
                    start = find_cache[ck]
                    if start < 0:
                        continue
                    end = start + len(pep)
                    for pos, sid in sites_by_acc[acc]:
                        if start < pos <= end:      # pos 为 1-based
                            denom[(sid, key)] += qty
                    break   # 一个前体只按第一个可定位蛋白计一次, 避免重复计数
        log(f"  [分母] {run} 完成")

    # ---------- 输出 ----------
    with open(os.path.join(outdir, "sites.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["site_id", "protein", "position", "residue", "ambiguous_with"])
        for sid in sorted(site_info):
            acc, pos, res = site_info[sid]
            w.writerow([sid, acc, pos, res, ";".join(sorted(site_ambig.get(sid, [])))])

    n_pairs_with_denom = 0
    with open(os.path.join(outdir, "occupancy_long.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["site_id", "timepoint", "fraction", "rep",
                    "phos_intensity", "prot_intensity", "occupancy_index"])
        for (sid, key) in sorted(numer):
            tp, fr, rep = key
            p = numer[(sid, key)]
            t = denom.get((sid, key))
            occ = round(p / t, 6) if t else ""
            if t:
                n_pairs_with_denom += 1
            w.writerow([sid, tp, fr, rep, round(p, 1), round(t, 1) if t else "", occ])

    total_pairs = len(numer)
    lines = [
        f"位点数: {len(site_info)}",
        f"(位点 x 样本) 观测对: {total_pairs}",
        f"其中有 proteome 分母的: {n_pairs_with_denom} ({n_pairs_with_denom/total_pairs*100:.1f}%)",
        f"未映射 FASTA 的磷酸化前体行: {skipped_unmapped}",
        f"多蛋白歧义位点: {len(site_ambig)}",
    ]
    with open(os.path.join(outdir, "summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    log("\n== 汇总 ==")
    for l in lines:
        log(l)
    log(f"输出目录: {outdir}")

if __name__ == "__main__":
    main()
