#!/usr/bin/env python3
"""09b_paired_patch.py — 补齐 09 缺失的蛋白内配对逐特征结果 (B vs C 同蛋白)。
每个可配对蛋白: 该蛋白 B 位点特征中位数 - C 位点特征中位数 -> Wilcoxon 符号秩 (正态近似) + BH。
输出: analysis07/paired_comparison.tsv + 屏幕摘要 (仅计数与数值)。
用法: python3 09b_paired_patch.py   (login 节点直跑, ~3 分钟)
"""
import sys, os, re, csv, json, glob, math
from collections import defaultdict

CHANNELS = ["backbone", "sidechain", "ppII", "coil", "sheet", "helix", "earlyFolding", "disoMine"]
FEATURES = [f"single_{c}" for c in CHANNELS] + [f"win5_{c}" for c in CHANNELS]
PAIR_MIN = int(os.environ.get("PAIR_MIN", "50"))
FORBIDDEN = ["证明", "验证了", "发现了机制", "breakthrough", "隐藏的信号", "hidden signal", "突破"]

def median(v):
    s = sorted(v); n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2

def rankdata(v):
    idx = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v); i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[idx[j + 1]] == v[idx[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[idx[k]] = avg
        i = j + 1
    return r

def wilcoxon_p(diffs):
    d = [x for x in diffs if x != 0]
    n = len(d)
    if n < 10:
        return None, n
    ranks = rankdata([abs(x) for x in d])
    wpos = sum(r for r, x in zip(ranks, d) if x > 0)
    mu = n * (n + 1) / 4
    var = n * (n + 1) * (2 * n + 1) / 24
    z = (wpos - mu) / math.sqrt(var)
    return math.erfc(abs(z) / math.sqrt(2)), n

def bh(pairs):
    n = len(pairs)
    o = sorted(pairs, key=lambda t: t[1])
    out, prev = {}, 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, o[i][1] * n / (i + 1))
        out[o[i][0]] = prev
    return out

def main():
    a7 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ["VSC_DATA"], "rerun_2026-08/occupancy/analysis07")
    jdir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.environ["VSC_DATA"], "biophys_json")
    fasta = sys.argv[3] if len(sys.argv) > 3 else sorted(glob.glob(os.path.join(os.environ["VSC_SCRATCH"], "rerun/fasta/human_*_UP*.fasta")))[-1]
    seqs, acc, buf = {}, None, []
    for line in open(fasta):
        line = line.rstrip()
        if line.startswith(">"):
            if acc: seqs[acc] = "".join(buf)
            m = re.match(r">\w+\|([^|]+)\|", line)
            acc = m.group(1) if m else line[1:].split()[0]; buf = []
        else: buf.append(line)
    if acc: seqs[acc] = "".join(buf)

    PRI = {"both": 3, "occupancy_only": 2, "intensity_only": 1, "none": 0}
    site_cls = {}
    with open(os.path.join(a7, "site_traj.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            p = PRI[r["class"]]
            if p > site_cls.get(r["site_id"], -1): site_cls[r["site_id"]] = p
    prof = {r["site_id"] for r in csv.DictReader(open(os.path.join(a7, "site_profiles.tsv")), delimiter="\t")}
    grp_B = {r["site_id"] for r in csv.DictReader(open(os.path.join(a7, "top_candidates.tsv")), delimiter="\t")}
    grp_C = {s for s, p in site_cls.items() if p == 0 and s in prof}

    need = grp_B | grp_C
    by_acc = defaultdict(list)
    for s in need:
        a, rp = s.rsplit("_", 1)
        by_acc[a].append((s, rp[0], int(rp[1:])))
    feats, nfail = {}, 0
    for a in sorted(by_acc):
        jp = os.path.join(jdir, a + ".json")
        if not os.path.isfile(jp): nfail += len(by_acc[a]); continue
        try:
            res = json.load(open(jp))["residues"]
        except Exception:
            nfail += len(by_acc[a]); continue
        pm = {r["seqpos"]: r for r in res}
        aa_str = "".join(r["aa"] for r in res)
        if a in seqs and aa_str != seqs[a]: nfail += len(by_acc[a]); continue
        for s, letter, pos in by_acc[a]:
            rr = pm.get(pos)
            if not rr or rr.get("aa") != letter: nfail += 1; continue
            row = {f"single_{c}": float(rr[c]) for c in CHANNELS}
            for c in CHANNELS:
                vv = [float(pm[p][c]) for p in range(pos - 5, pos + 6) if p in pm]
                row[f"win5_{c}"] = sum(vv) / len(vv)
            feats[s] = row
    print(f"位点 {len(need)} | 特征成功 {len(feats)} | 失败 {nfail}")

    pB, pC = defaultdict(list), defaultdict(list)
    for s in grp_B:
        if s in feats: pB[s.rsplit("_", 1)[0]].append(s)
    for s in grp_C:
        if s in feats: pC[s.rsplit("_", 1)[0]].append(s)
    common = sorted(set(pB) & set(pC))
    print(f"可配对蛋白 (B∩C): {len(common)}")
    if len(common) < PAIR_MIN:
        print(f"n < {PAIR_MIN}, 仅报计数"); sys.exit(0)

    rows, praw = [], []
    for feat in FEATURES:
        diffs = [median([feats[s][feat] for s in pB[a]]) - median([feats[s][feat] for s in pC[a]]) for a in common]
        p, n_nz = wilcoxon_p(diffs)
        npos = sum(1 for d in diffs if d > 0)
        rows.append([feat, len(common), n_nz, round(median(diffs), 5), npos, len(common) - npos,
                     f"{p:.3e}" if p is not None else ""])
        if p is not None: praw.append((feat, p))
    pbh = bh(praw)
    out = os.path.join(a7, "paired_comparison.tsv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["feature", "n_paired_proteins", "n_nonzero", "median_paired_diff_BminusC",
                    "n_positive", "n_negative", "wilcoxon_p", "p_BH"])
        for r in rows:
            w.writerow(r + [f"{pbh[r[0]]:.3e}" if r[0] in pbh else ""])
    nsig = sum(1 for f_, p_ in pbh.items() if p_ < 0.05)
    lines = [f"蛋白内配对 B vs C: 可配对蛋白 {len(common)} | BH<0.05 特征 {nsig}/16"]
    for r in sorted(rows, key=lambda r: abs(r[3]), reverse=True)[:4]:
        s_ = "+" if r[3] > 0 else ""
        lines.append(f"  {r[0]}: 中位配对差 {s_}{r[3]} ({r[4]}正/{r[5]}负) BH p={pbh.get(r[0], float('nan')):.2e}")
    txt = "\n".join(lines)
    for wbad in FORBIDDEN:
        assert wbad not in txt, f"禁用词: {wbad}"
    print(txt)
    with open(os.path.join(a7, "feature_summary.txt"), "a") as f:
        f.write("\n[09b 补丁] " + txt + "\n")

if __name__ == "__main__":
    main()
