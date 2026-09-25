#!/usr/bin/env python3
"""06_occupancy_variants.py — 三种分母方案同跑, 产出对比表供方法学讨论。

方案:
  v1  严格残基级: 同 (timepoint, fraction, rep) 样本中覆盖该残基的肽段强度和 (= 05 的定义)
  vA  残基级 + 蛋白级回退: v1 缺失时改用同样本中该蛋白全部肽段强度和, fallback_type 标注
  vB  跨重复合并残基级: 分母 = 同 (timepoint, fraction) 下各重复残基级分母的均值 (仅对有值的重复取均值)

附加输出: 残基级与蛋白级分母同时存在的行上, log10(occ_res) vs log10(occ_prot) 的 Pearson r
          —— 衡量回退方案与严格方案的一致性, 供 Paddy 讨论。

用法: python3 06_occupancy_variants.py [结果根目录] [FASTA路径]   (默认同 05)
输出: <根目录>/occupancy/{occupancy_variants_long.tsv, summary_variants.txt}
"""
import sys, os, re, csv, glob, math
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
    seqs, acc, buf = {}, None, []
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
    sites, pos, last, i = [], 0, "", 0
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
    c["qval"] = pick_col(cols, ["Q.Value"], lambda x: x.endswith("Q.Value") and all(t not in x for t in ("PG", "PTM", "Lib", "Global")))
    c["conf"] = pick_col(cols, ["PTM.Site.Confidence"], lambda x: "Site.Confidence" in x) if need_conf else None
    return c

def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)

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
    phos_dirs = sorted(glob.glob(os.path.join(base, "phospho", "*/")))
    prot_dirs = sorted(glob.glob(os.path.join(base, "proteome", "*/")))
    log(f"FASTA {len(seqs)} 条 | phospho {len(phos_dirs)} | proteome {len(prot_dirs)}")
    if not phos_dirs or not prot_dirs:
        log("FAIL: run 目录缺失"); sys.exit(1)

    find_cache = {}
    site_info = {}
    numer = defaultdict(float)

    # ---- 分子 (同 05) ----
    for d in phos_dirs:
        run = os.path.basename(d.rstrip("/"))
        key = parse_design(run)
        if key is None:
            continue
        pf, cols = open_report(d)
        if pf is None:
            log(f"WARN: 缺 parquet {run}"); continue
        c = resolve_cols(cols, need_conf=True)
        if any(not c[k] for k in ("mod", "strip", "prot", "qty", "conf")):
            log(f"FAIL: {run} 列缺失; 实际列: {cols}"); sys.exit(1)
        use = [v for v in c.values() if v]
        for b in pf.iter_batches(columns=use, batch_size=100_000):
            cd = {n: b.column(n).to_pylist() for n in use}
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
                    if find_cache[ck] >= 0:
                        placed = (acc, find_cache[ck]); break
                if placed is None:
                    continue
                acc, start = placed
                for ppos, res in psites:
                    sid = f"{acc}_{res}{start + ppos}"
                    site_info[sid] = (acc, start + ppos, res)
                    numer[(sid, key)] += qty
        log(f"  [分子] {run}")
    log(f"位点数: {len(site_info)}")
    if not site_info:
        log("FAIL: 零位点"); sys.exit(1)

    sites_by_acc = defaultdict(list)
    for sid, (acc, pos, _r) in site_info.items():
        sites_by_acc[acc].append((pos, sid))

    # ---- 分母: 残基级 + 蛋白级, 一遍扫完 ----
    denom_res = defaultdict(float)    # (sid, key)
    denom_prot = defaultdict(float)   # (acc, key)
    for d in prot_dirs:
        run = os.path.basename(d.rstrip("/"))
        key = parse_design(run)
        if key is None:
            continue
        pf, cols = open_report(d)
        if pf is None:
            log(f"WARN: 缺 parquet {run}"); continue
        c = resolve_cols(cols, need_conf=False)
        if any(not c[k] for k in ("strip", "prot", "qty")):
            log(f"FAIL: {run} 列缺失; 实际列: {cols}"); sys.exit(1)
        use = [v for k, v in c.items() if v and k != "mod"]
        for b in pf.iter_batches(columns=use, batch_size=100_000):
            cd = {n: b.column(n).to_pylist() for n in use}
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
                    denom_prot[(acc, key)] += qty
                    end = start + len(pep)
                    for pos, sid in sites_by_acc[acc]:
                        if start < pos <= end:
                            denom_res[(sid, key)] += qty
                    break
        log(f"  [分母] {run}")

    # vB: 残基级跨重复均值, keyed (sid, (tp, fr))
    res_by_tf = defaultdict(list)
    for (sid, (tp, fr, rep)), v in denom_res.items():
        res_by_tf[(sid, (tp, fr))].append(v)
    denom_res_pooled = {k: sum(v) / len(v) for k, v in res_by_tf.items()}

    # ---- 输出长表 ----
    n = {"total": 0, "v1": 0, "vA": 0, "vA_fb": 0, "vB": 0}
    log_pairs_x, log_pairs_y = [], []
    outp = os.path.join(outdir, "occupancy_variants_long.tsv")
    with open(outp, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["site_id", "timepoint", "fraction", "rep", "phos_intensity",
                    "denom_residue", "denom_protein", "denom_residue_repmean",
                    "occ_v1_strict", "occ_vA_fallback", "vA_type", "occ_vB_repmean"])
        for (sid, key) in sorted(numer):
            tp, fr, rep = key
            p = numer[(sid, key)]
            acc = site_info[sid][0]
            dr = denom_res.get((sid, key))
            dp = denom_prot.get((acc, key))
            dm = denom_res_pooled.get((sid, (tp, fr)))
            occ1 = round(p / dr, 6) if dr else ""
            if dr:
                fb_type, occA = "residue", round(p / dr, 6)
            elif dp:
                fb_type, occA = "protein", round(p / dp, 6)
            else:
                fb_type, occA = "none", ""
            occB = round(p / dm, 6) if dm else ""
            n["total"] += 1
            n["v1"] += bool(dr)
            n["vA"] += fb_type != "none"
            n["vA_fb"] += fb_type == "protein"
            n["vB"] += bool(dm)
            if dr and dp:
                log_pairs_x.append(math.log10(p / dr))
                log_pairs_y.append(math.log10(p / dp))
            w.writerow([sid, tp, fr, rep, round(p, 1),
                        round(dr, 1) if dr else "", round(dp, 1) if dp else "",
                        round(dm, 1) if dm else "", occ1, occA, fb_type, occB])

    r = pearson(log_pairs_x, log_pairs_y)
    lines = [
        f"位点数: {len(site_info)} | 观测对: {n['total']}",
        f"v1 严格残基级覆盖:        {n['v1']} ({n['v1']/n['total']*100:.1f}%)",
        f"vA 残基+蛋白回退覆盖:     {n['vA']} ({n['vA']/n['total']*100:.1f}%)  其中回退行 {n['vA_fb']} ({n['vA_fb']/n['total']*100:.1f}%)",
        f"vB 跨重复均值残基级覆盖:  {n['vB']} ({n['vB']/n['total']*100:.1f}%)",
        f"一致性: 双分母共存行 n={len(log_pairs_x)}, log10 occupancy Pearson r = {r:.3f}" if r is not None else "一致性: 双分母共存行不足, 无法计算",
    ]
    with open(os.path.join(outdir, "summary_variants.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    log("\n== 方案对比 ==")
    for l in lines:
        log(l)
    log(f"长表: {outp}")

if __name__ == "__main__":
    main()
