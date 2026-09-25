#!/usr/bin/env python3
"""13_breadth_features.py — Task C: compartment breadth 分组 + b2bTools 特征比较 (C1+C2, C3 条件出图)。
C1: breadth = CTRL 样本中位点被检出(q<=0.01, UniMod:21, loc>=0.75)且 >=2/4 重复 的 fraction 数
    组: G1_specific(=1) / G2_moderate(2-5) / G3_ubiquitous(=6); 同时输出 2/3/4/5 细分计数
    丰度检查内嵌: 各组 CTRL 总强度(log2)中位数 + breadth x log2强度 Spearman + 强度十分位匹配子集复检
C2: 16 特征 (8 通道 x 单点/±5窗, biophys_json, seqpos 索引 + 残基三方互证)
统计: Kruskal-Wallis (df=2, p=exp(-H/2) 精确) + epsilon^2 + BH(16); 蛋白级聚合敏感性
C3: matplotlib 可用则出 4x4 分面 boxplot PNG, 否则只出 tsv
输入: phospho CTRL 24 run parquet + FASTA + $VSC_DATA/biophys_json
输出: analysis07/{breadth_features.tsv, breadth_summary.txt[, breadth_boxplots.png]}
预计: sbatch zen4 30 分钟档, 实际 ~5-10 分钟
"""
import os, re, csv, glob, json, math, random
from collections import defaultdict, Counter

LOC, QV, MIN_REPS = 0.75, 0.01, 2
CH = ["backbone","sidechain","ppII","coil","sheet","helix","earlyFolding","disoMine"]
FEATS = [f"single_{c}" for c in CH] + [f"win5_{c}" for c in CH]
FRS = [f"FR{i}" for i in range(1,7)]
FORBIDDEN = ["证明","验证了","breakthrough","突破","隐藏的信号","hidden signal"]

def parse_design(run):
    tp = re.search(r"_(2min|8min|20min|90min|CTRL)_", run, re.I)
    fr = re.search(r"_(FR\d)_", run); rep = re.search(r"_(Rep\d)", run)
    return (tp.group(1), fr.group(1), rep.group(1)) if (tp and fr and rep) else None

def load_fasta(path):
    seqs, acc, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if acc: seqs[acc] = "".join(buf)
            m = re.match(r">\w+\|([^|]+)\|", line); acc = m.group(1) if m else line[1:].split()[0]; buf = []
        else: buf.append(line)
    if acc: seqs[acc] = "".join(buf)
    return seqs

def psites(mod):
    out, pos, last, i = [], 0, "", 0
    while i < len(mod):
        c = mod[i]
        if c == "(":
            j = mod.index(")", i)
            if mod[i+1:j] == "UniMod:21" and pos > 0: out.append((pos, last))
            i = j + 1
        else:
            if c.isalpha(): pos += 1; last = c
            i += 1
    return out

def rankdata(v):
    idx = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0]*len(v); i = 0
    while i < len(v):
        j = i
        while j+1 < len(v) and v[idx[j+1]] == v[idx[i]]: j += 1
        for k in range(i, j+1): r[idx[k]] = (i+j)/2 + 1
        i = j + 1
    return r

def kw(groups):
    all_v = [x for g in groups for x in g]
    n = len(all_v)
    if n < 6 or any(len(g) < 2 for g in groups): return None, None
    ranks = rankdata(all_v)
    H, off = 0.0, 0
    for g in groups:
        rs = sum(ranks[off:off+len(g)]); off += len(g)
        H += rs*rs/len(g)
    H = 12/(n*(n+1))*H - 3*(n+1)
    tie = Counter(all_v); T = sum(t**3-t for t in tie.values())
    if T: H /= (1 - T/(n**3-n))
    p = math.exp(-H/2)                      # df=2 精确
    eps2 = H*(n+1)/(n*n-1)
    return p, eps2

def spearman(x, y):
    rx, ry = rankdata(x), rankdata(y)
    n = len(x); mx, my = sum(rx)/n, sum(ry)/n
    sxy = sum((a-mx)*(b-my) for a,b in zip(rx,ry))
    sxx = sum((a-mx)**2 for a in rx); syy = sum((b-my)**2 for b in ry)
    return sxy/math.sqrt(sxx*syy) if sxx*syy > 0 else float("nan")

def bh(ps):
    n = len(ps); order = sorted(range(n), key=lambda i: ps[i]); out = [0.0]*n; prev = 1.0
    for k, i in enumerate(reversed(order)):
        prev = min(prev, ps[i]*n/(n-k)); out[i] = prev
    return out

def main():
    import pyarrow.parquet as pq
    base = os.path.join(os.environ["VSC_DATA"], "rerun_2026-08")
    jdir = os.path.join(os.environ["VSC_DATA"], "biophys_json")
    fasta = sorted(glob.glob(os.path.join(os.environ["VSC_SCRATCH"], "rerun/fasta/human_*_UP*.fasta")))[-1]
    a7 = os.path.join(base, "occupancy/analysis07")
    seqs = load_fasta(fasta); find_cache = {}
    det = defaultdict(lambda: defaultdict(set))   # sid -> fr -> reps
    inten = defaultdict(float); site_acc = {}
    for d in sorted(glob.glob(os.path.join(base, "phospho", "*CTRL*/"))):
        run = os.path.basename(d.rstrip("/")); key = parse_design(run)
        if key is None or key[0] != "CTRL": continue
        _, fr, rep = key
        pf = pq.ParquetFile(os.path.join(d, "report.parquet"))
        cols = pf.schema_arrow.names
        use = [c for c in ("Modified.Sequence","Stripped.Sequence","Protein.Ids","Precursor.Quantity","Q.Value","PTM.Site.Confidence") if c in cols]
        for b in pf.iter_batches(columns=use, batch_size=100_000):
            cd = {n: b.column(n).to_pylist() for n in use}
            for i in range(b.num_rows):
                if (cd["Q.Value"][i] or 1) > QV: continue
                mod = cd["Modified.Sequence"][i] or ""
                if "UniMod:21" not in mod: continue
                conf = cd["PTM.Site.Confidence"][i]
                if conf is None or conf < LOC: continue
                qty = cd["Precursor.Quantity"][i] or 0
                pep = cd["Stripped.Sequence"][i] or ""
                accs = [a for a in re.split(r"[;,]", cd["Protein.Ids"][i] or "") if a]
                ps = psites(mod)
                if not pep or not accs or not ps: continue
                placed = None
                for acc in accs:
                    s = seqs.get(acc)
                    if s is None: continue
                    ck = (acc, pep)
                    if ck not in find_cache: find_cache[ck] = s.find(pep)
                    if find_cache[ck] >= 0: placed = (acc, find_cache[ck]); break
                if placed is None: continue
                acc, st = placed
                for pp, r in ps:
                    sid = f"{acc}_{r}{st+pp}"
                    site_acc[sid] = acc
                    det[sid][fr].add(rep)
                    if qty > 0: inten[sid] += qty
        print(f"  [scan] {run}", flush=True)

    breadth = {s: sum(1 for fr, reps in frs.items() if len(reps) >= MIN_REPS) for s, frs in det.items()}
    breadth = {s: b for s, b in breadth.items() if b >= 1}
    grp = {s: ("G1_specific" if b == 1 else "G3_ubiquitous" if b == 6 else "G2_moderate") for s, b in breadth.items()}
    print(f"CTRL 中 breadth>=1 的位点: {len(breadth)} | 分布: {sorted(Counter(breadth.values()).items())}")

    # 特征提取
    by_acc = defaultdict(list)
    for s in breadth:
        a, rp = s.rsplit("_", 1); by_acc[a].append((s, rp[0], int(rp[1:])))
    feats, nfail = {}, 0
    for a in sorted(by_acc):
        jp = os.path.join(jdir, a + ".json")
        try: res = json.load(open(jp))["residues"]
        except Exception: nfail += len(by_acc[a]); continue
        pm = {r["seqpos"]: r for r in res}
        if "".join(r["aa"] for r in res) != seqs.get(a, ""): nfail += len(by_acc[a]); continue
        for s, letter, pos in by_acc[a]:
            rr = pm.get(pos)
            if not rr or rr.get("aa") != letter: nfail += 1; continue
            row = {f"single_{c}": float(rr[c]) for c in CH}
            for c in CH:
                vv = [float(pm[p][c]) for p in range(pos-5, pos+6) if p in pm]
                row[f"win5_{c}"] = sum(vv)/len(vv)
            feats[s] = row
    print(f"特征成功 {len(feats)} | 失败 {nfail}")

    G = {g: [s for s in feats if grp[s] == g] for g in ("G1_specific","G2_moderate","G3_ubiquitous")}
    li = {s: math.log2(inten[s]) for s in feats if inten.get(s, 0) > 0}

    S = [f"各组位点数: " + " | ".join(f"{g}={len(v)}" for g, v in G.items()),
         f"breadth 细分: " + " ".join(f"{b}:{c}" for b, c in sorted(Counter(breadth[s] for s in feats).items()))]
    med = lambda v: sorted(v)[len(v)//2] if v else float("nan")
    S.append("丰度检查: " + " | ".join(f"{g} 中位log2强度={med([li[s] for s in v if s in li]):.2f}" for g, v in G.items()))
    xs = [breadth[s] for s in feats if s in li]; ys = [li[s] for s in feats if s in li]
    rho = spearman(xs, ys)
    S.append(f"breadth x log2强度 Spearman rho = {rho:.3f} (n={len(xs)})")

    # 十分位匹配子集
    rng = random.Random(20260825)
    vals = sorted(li.values()); dec = [vals[int(len(vals)*k/10)] for k in range(1, 10)]
    def decile(x):
        d = 0
        for t in dec:
            if x > t: d += 1
        return d
    bins = {g: defaultdict(list) for g in G}
    for g, ss in G.items():
        for s in ss:
            if s in li: bins[g][decile(li[s])].append(s)
    matched = {g: [] for g in G}
    for d in range(10):
        m = min(len(bins[g][d]) for g in G)
        for g in G: matched[g] += rng.sample(bins[g][d], m)
    S.append(f"十分位匹配子集: " + " | ".join(f"{g}={len(v)}" for g, v in matched.items()))

    def run_stats(groups_dict, label):
        ps, rows = [], []
        for feat in FEATS:
            gs = [[feats[s][feat] for s in groups_dict[g]] for g in ("G1_specific","G2_moderate","G3_ubiquitous")]
            p, e2 = kw(gs)
            if p is None: continue
            ps.append(p); rows.append([label, feat, round(e2, 4), p] + [round(med(g), 4) for g in gs])
        for pb, r in zip(bh(ps), rows): r.append(pb)
        return rows

    rows = run_stats(G, "full") + run_stats(matched, "abundance_matched")
    # 蛋白级聚合敏感性 (full)
    prot_g = {g: defaultdict(list) for g in G}
    for g, ss in G.items():
        for s in ss: prot_g[g][site_acc[s]].append(s)
    protG = {g: None for g in G}
    prows = []
    for feat in FEATS:
        gs = [[med([feats[s][feat] for s in sl]) for sl in prot_g[g].values()] for g in ("G1_specific","G2_moderate","G3_ubiquitous")]
        p, e2 = kw(gs)
        if p is not None: prows.append(["protein_level", feat, round(e2, 4), p] + [round(med(g), 4) for g in gs])
    pps = [r[3] for r in prows]
    for pb, r in zip(bh(pps), prows): r.append(pb)
    rows += prows

    with open(os.path.join(a7, "breadth_features.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["analysis","feature","epsilon2","p_KW","median_G1","median_G2","median_G3","p_BH"])
        for r in rows: w.writerow(r[:3] + [f"{r[3]:.3e}"] + r[4:7] + [f"{r[7]:.3e}"])
    with open(os.path.join(a7, "breadth_sites.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["site_id","protein","breadth","group","log2_intensity_ctrl"] + FEATS)
        for s in sorted(feats):
            w.writerow([s, site_acc[s], breadth[s], grp[s], round(li.get(s, float("nan")), 3)] + [round(feats[s][ft], 4) for ft in FEATS])

    for label in ("full","abundance_matched","protein_level"):
        sub = [r for r in rows if r[0] == label]
        nsig = sum(1 for r in sub if r[-1] < 0.05)
        top = max(sub, key=lambda r: r[2]) if sub else None
        S.append(f"{label}: BH<0.05 特征 {nsig}/{len(sub)}" + (f" | 最大 eps2={top[2]} ({top[1]})" if top else ""))
    txt = "\n".join(S)
    for wbad in FORBIDDEN: assert wbad not in txt, wbad
    with open(os.path.join(a7, "breadth_summary.txt"), "w") as f: f.write(txt + "\n")
    print("\n" + txt)

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, axes = plt.subplots(4, 4, figsize=(14, 12))
        for ax, feat in zip(axes.flat, FEATS):
            data = [[feats[s][feat] for s in G[g]] for g in ("G1_specific","G2_moderate","G3_ubiquitous")]
            ax.boxplot(data, tick_labels=["1fr","2-5fr","6fr"], showfliers=False)
            ax.set_title(feat, fontsize=8)
        fig.suptitle("b2bTools features by compartment breadth (CTRL, >=2/4 reps)")
        fig.tight_layout(); fig.savefig(os.path.join(a7, "breadth_boxplots.png"), dpi=150)
        print("PNG 已出")
    except ImportError:
        print("matplotlib 不可用 — 图交给 Codex (用 breadth_sites.tsv)")

if __name__ == "__main__":
    main()
