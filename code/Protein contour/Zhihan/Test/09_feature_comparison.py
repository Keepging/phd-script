#!/usr/bin/env python3
"""09_feature_comparison.py — 候选位点 vs 背景位点的 per-residue 特征对比 (spec v2)。

通道 (8, 含 ppII): backbone sidechain ppII coil sheet helix earlyFolding disoMine
特征 16 = 8 单点 + 8 (±5 窗口均值)
分组: A(严格候选) B(宽松候选) C(位点级 neither 且有轮廓数据) ONI(occ_not_int)
      INO(int_not_occ) ALLexA(全体去 A)
对比: 1 A-C | 2 B-C | 3 A-ALLexA | 4 ONI-INO | 5 ONI-C  (MWU + 组内 BH + 蛋白级置换)
补充: Spearman 连续版 | 蛋白内配对 (B vs C 同蛋白, Wilcoxon, n<50 只报计数) | IDR Fisher
QC: 蛋白聚集度 | 残基组成 | 窗口重叠 | 映射失败四分类 (>5% 停机)
对齐: 按 JSON seqpos 建索引; 残基字母三方互证 (site_id / JSON aa / FASTA)

用法: python3 09_feature_comparison.py [analysis07目录] [b2b_json目录] [FASTA]
"""
import sys, os, csv, json, glob, math, random
from collections import defaultdict

CHANNELS = ["backbone", "sidechain", "ppII", "coil", "sheet", "helix", "earlyFolding", "disoMine"]
FEATURES = [f"single_{c}" for c in CHANNELS] + [f"win5_{c}" for c in CHANNELS]
N_PERM = int(os.environ.get("PERMS", "1000"))
MIN_GROUP = int(os.environ.get("MIN_GROUP", "20"))
MIN_PAIRED = 50
IDR_CH = "single_disoMine"
FORBIDDEN = ["证明", "验证了", "发现了机制", "breakthrough", "隐藏的信号", "hidden signal", "突破"]
random.seed(20260821)

def log(m):
    print(m, flush=True)

def load_fasta(path):
    seqs, acc, buf = {}, None, []
    import re
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

def median(v):
    s = sorted(v)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2

def bh(pvals):
    """Benjamini-Hochberg; 输入 list[(idx, p)], 返回 idx->p_bh"""
    items = sorted(pvals, key=lambda x: x[1])
    n = len(items)
    out, prev = {}, 1.0
    for rank in range(n, 0, -1):
        idx, p = items[rank - 1]
        val = min(prev, p * n / rank)
        out[idx] = val
        prev = val
    return out

def rankdata(vals):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = r
        i = j + 1
    return ranks

def mwu(g1, g2):
    """返回 U1, p(正态近似含 tie 校正), rbc"""
    n1, n2 = len(g1), len(g2)
    allv = g1 + g2
    ranks = rankdata(allv)
    r1 = sum(ranks[:n1])
    u1 = r1 - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2
    # tie 校正
    cnt = defaultdict(int)
    for v in allv:
        cnt[v] += 1
    n = n1 + n2
    tie = sum(c ** 3 - c for c in cnt.values())
    var = n1 * n2 / 12 * ((n + 1) - tie / (n * (n - 1))) if n > 1 else 0
    if var <= 0:
        return u1, 1.0, 0.0
    z = (u1 - mu - (0.5 if u1 > mu else -0.5)) / math.sqrt(var)
    p = math.erfc(abs(z) / math.sqrt(2))
    rbc = 2 * u1 / (n1 * n2) - 1
    return u1, p, rbc

def wilcoxon_signed(diffs):
    d = [x for x in diffs if x != 0]
    n = len(d)
    if n < 10:
        return None, None
    ranks = rankdata([abs(x) for x in d])
    wpos = sum(r for r, x in zip(ranks, d) if x > 0)
    mu = n * (n + 1) / 4
    var = n * (n + 1) * (2 * n + 1) / 24
    z = (wpos - mu) / math.sqrt(var)
    return wpos, math.erfc(abs(z) / math.sqrt(2))

def fisher(a, b, c, d):
    """2x2 精确检验 (a=g1_idr b=g1_non c=g2_idr d=g2_non), 超几何双侧"""
    def lchoose(n, k):
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    row1, row2, col1 = a + b, c + d, a + c
    n = row1 + row2
    def pmf(x):
        return math.exp(lchoose(row1, x) + lchoose(row2, col1 - x) - lchoose(n, col1))
    p_obs = pmf(a)
    p = sum(pmf(x) for x in range(max(0, col1 - row2), min(row1, col1) + 1) if pmf(x) <= p_obs * (1 + 1e-9))
    odds = (a * d) / (b * c) if b * c > 0 else float("inf")
    return odds, min(1.0, p)

def main():
    a7 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ["VSC_DATA"], "rerun_2026-08/occupancy/analysis07")
    jdir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.environ["VSC_DATA"], "biophys_json")
    if len(sys.argv) > 3:
        fasta_path = sys.argv[3]
    else:
        cands = sorted(glob.glob(os.path.join(os.environ["VSC_SCRATCH"], "rerun/fasta/human_*_UP*.fasta")))
        fasta_path = cands[-1] if cands else ""
    if not os.path.isdir(jdir) or not os.listdir(jdir):
        log(f"FAIL: JSON 目录不存在或为空: {jdir} — 请先上传"); sys.exit(1)
    if not os.path.isfile(fasta_path):
        log(f"FAIL: FASTA 不存在: {fasta_path}"); sys.exit(1)
    seqs = load_fasta(fasta_path)

    # ---- 位点级分组 ----
    PRI = {"both": 3, "occupancy_only": 2, "intensity_only": 1, "none": 0}
    NAME = {3: "occ_and_int", 2: "occ_not_int", 1: "int_not_occ", 0: "neither"}
    site_cls, site_occfc = {}, {}
    with open(os.path.join(a7, "site_traj.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sid = r["site_id"]
            pri = PRI[r["class"]]
            if pri > site_cls.get(sid, -1):
                site_cls[sid] = pri
            try:
                base = float(r["log2occ_CTRL"])
                for tp in ("2min", "8min", "20min", "90min"):
                    v = r.get(f"log2occ_{tp}", "")
                    if v != "":
                        fc = abs(float(v) - base)
                        if fc > site_occfc.get(sid, -1):
                            site_occfc[sid] = fc
            except ValueError:
                pass
    prof_sites, site_shift = set(), {}
    with open(os.path.join(a7, "site_profiles.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sid = r["site_id"]
            prof_sites.add(sid)
            if r["phos_shift_vs_ctrl"] != "":
                v = float(r["phos_shift_vs_ctrl"])
                if v > site_shift.get(sid, -1):
                    site_shift[sid] = v
    grp_A, grp_B = set(), set()
    with open(os.path.join(a7, "top_candidates.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            grp_B.add(r["site_id"])
            dcol = next((k for k in r if k.startswith("divergence_ge")), None)
            if r["residue_evidence"] == "1" and dcol and r[dcol] == "1":
                grp_A.add(r["site_id"])
    all_sites = set(site_cls)
    grp_C = {s for s in all_sites if NAME[site_cls[s]] == "neither" and s in prof_sites}
    grp_ONI = {s for s in all_sites if NAME[site_cls[s]] == "occ_not_int"}
    grp_INO = {s for s in all_sites if NAME[site_cls[s]] == "int_not_occ"}
    grp_ALLexA = all_sites - grp_A
    log(f"分组: A={len(grp_A)} B={len(grp_B)} C={len(grp_C)} ONI={len(grp_ONI)} INO={len(grp_INO)} ALLexA={len(grp_ALLexA)}")

    # ---- JSON 加载 + 三方对齐 + 特征提取 ----
    need_prots = {s.rsplit("_", 1)[0] for s in all_sites}
    feats = {}          # sid -> {feature: value}
    qc_rows = []
    prot_json = {}      # acc -> {seqpos: featdict}, 已验证
    bad_prot = {}       # acc -> 失败类型
    for acc in sorted(need_prots):
        jp = os.path.join(jdir, acc + ".json")
        if not os.path.isfile(jp):
            bad_prot[acc] = "missing_json"; continue
        try:
            data = json.load(open(jp))
            residues = data["residues"]
            aa_str = "".join(r["aa"] for r in residues)
        except Exception:
            bad_prot[acc] = "missing_json"; continue
        if acc in seqs and aa_str != seqs[acc]:
            bad_prot[acc] = "sequence_mismatch"; continue
        prot_json[acc] = {r["seqpos"]: r for r in residues}

    n_fail = 0
    for sid in sorted(all_sites):
        acc, tail = sid.rsplit("_", 1)
        letter, pos = tail[0], int(tail[1:])
        if acc in bad_prot:
            qc_rows.append([sid, bad_prot[acc]]); n_fail += 1; continue
        rd = prot_json.get(acc, {})
        if pos not in rd:
            qc_rows.append([sid, "position_out_of_range"]); n_fail += 1; continue
        if rd[pos]["aa"] != letter or (acc in seqs and seqs[acc][pos - 1] != letter):
            qc_rows.append([sid, "residue_mismatch"]); n_fail += 1; continue
        fv = {}
        for c in CHANNELS:
            fv[f"single_{c}"] = rd[pos][c]
        lo, hi = pos - 5, pos + 5
        wins = [rd[p] for p in range(lo, hi + 1) if p in rd]
        for c in CHANNELS:
            fv[f"win5_{c}"] = sum(w[c] for w in wins) / len(wins)
        feats[sid] = fv
    with open(os.path.join(a7, "mapping_qc.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["site_id", "failure_type"]); w.writerows(qc_rows)
    fail_rate = n_fail / len(all_sites) * 100
    log(f"映射: 成功 {len(feats)}/{len(all_sites)} | 失败 {n_fail} ({fail_rate:.1f}%)")
    if fail_rate > 5:
        log("FAIL: 映射失败率 > 5%, 按 spec 停止。失败明细见 mapping_qc.tsv"); sys.exit(1)

    def gv(grp, feat):
        return [feats[s][feat] for s in grp if s in feats]

    def gsites(grp):
        return [s for s in grp if s in feats]

    # ---- 主分析: MWU + 置换 ----
    comps = [("A_vs_C", grp_A, grp_C), ("B_vs_C", grp_B, grp_C),
             ("A_vs_ALLexA", grp_A, grp_ALLexA),
             ("ONI_vs_INO", grp_ONI, grp_INO), ("ONI_vs_C", grp_ONI, grp_C)]
    fc_rows = []
    summary_sig = {}
    top_by_rbc = []
    for cname, g1, g2 in comps:
        s1, s2 = gsites(g1 - g2), gsites(g2 - g1)
        dual_prots = {x.rsplit("_", 1)[0] for x in s1} & {x.rsplit("_", 1)[0] for x in s2}
        pvals = []
        crows = []
        for feat in FEATURES:
            v1, v2 = [feats[s][feat] for s in s1], [feats[s][feat] for s in s2]
            if len(v1) < MIN_GROUP or len(v2) < MIN_GROUP:
                crows.append([feat, cname, len(v1), len(v2), "", "", "", "", "", "", ""]); continue
            u1, p, rbc = mwu(v1, v2)
            # 蛋白级置换 (剔除双属蛋白)
            p1 = [s for s in s1 if s.rsplit("_", 1)[0] not in dual_prots]
            p2 = [s for s in s2 if s.rsplit("_", 1)[0] not in dual_prots]
            perm_p = ""
            if len(p1) >= MIN_GROUP and len(p2) >= MIN_GROUP:
                allv = [feats[s][feat] for s in p1] + [feats[s][feat] for s in p2]
                ranks = rankdata(allv)
                sidlist = p1 + p2
                by_prot = defaultdict(lambda: [0.0, 0])
                for i, s in enumerate(sidlist):
                    a = s.rsplit("_", 1)[0]
                    by_prot[a][0] += ranks[i]; by_prot[a][1] += 1
                prots = list(by_prot)
                nprot1 = len({s.rsplit("_", 1)[0] for s in p1})
                n1p, n2p = len(p1), len(p2)
                real_u = sum(ranks[:n1p]) - n1p * (n1p + 1) / 2
                real_rbc = 2 * real_u / (n1p * n2p) - 1
                hit = 0
                for _ in range(N_PERM):
                    random.shuffle(prots)
                    rsum = cnt = 0
                    for a in prots[:nprot1]:
                        rsum += by_prot[a][0]; cnt += by_prot[a][1]
                    if cnt == 0 or cnt == len(sidlist):
                        continue
                    m2 = len(sidlist) - cnt
                    pu = rsum - cnt * (cnt + 1) / 2
                    prbc = 2 * pu / (cnt * m2) - 1
                    if abs(prbc) >= abs(real_rbc) - 1e-12:
                        hit += 1
                perm_p = round((hit + 1) / (N_PERM + 1), 4)
            pvals.append((feat, p))
            crows.append([feat, cname, len(v1), len(v2), round(median(v1), 4), round(median(v2), 4),
                          round(u1, 1), f"{p:.3e}", None, round(rbc, 4), perm_p])
            if cname == "A_vs_C":
                top_by_rbc.append((abs(rbc), feat, rbc, v1, v2))
        pbh = bh(pvals) if pvals else {}
        nsig = nsig_perm = 0
        for row in crows:
            if row[7] != "":
                row[8] = f"{pbh[row[0]]:.3e}"
                if pbh[row[0]] < 0.05:
                    nsig += 1
                if row[10] != "" and row[10] < 0.05:
                    nsig_perm += 1
        summary_sig[cname] = (nsig, nsig_perm, len(dual_prots))
        fc_rows += crows
        log(f"  [对比] {cname}: BH<0.05 {nsig}/16, 置换p<0.05 {nsig_perm}")
    with open(os.path.join(a7, "feature_comparison.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["feature", "comparison", "n_group1", "n_group2", "median_g1", "median_g2",
                    "U1", "p_raw", "p_BH", "effect_size_rbc", "permutation_p"])
        w.writerows(fc_rows)

    # ---- Spearman 连续版 ----
    def spearman(xs, ys):
        rx, ry = rankdata(xs), rankdata(ys)
        n = len(xs)
        mx, my = sum(rx) / n, sum(ry) / n
        sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
        sxx = sum((a - mx) ** 2 for a in rx)
        syy = sum((b - my) ** 2 for b in ry)
        if sxx <= 0 or syy <= 0:
            return 0.0, 1.0
        rho = sxy / math.sqrt(sxx * syy)
        t = rho * math.sqrt((n - 2) / max(1e-12, 1 - rho ** 2))
        p = math.erfc(abs(t) / math.sqrt(2))   # 大样本正态近似
        return rho, p
    sp_rows = []
    n_sp_sig = 0
    for varname, mapping in (("max_abs_occ_log2fc", site_occfc), ("max_profile_shift", site_shift)):
        ss = [s for s in prof_sites if s in feats and s in mapping]
        for feat in FEATURES:
            xs = [feats[s][feat] for s in ss]
            ys = [mapping[s] for s in ss]
            rho, p = spearman(xs, ys)
            sp_rows.append([feat, varname, len(ss), round(rho, 4), f"{p:.3e}"])
            if p < 0.05 / 32:
                n_sp_sig += 1
    with open(os.path.join(a7, "spearman_continuous.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["feature", "continuous_var", "n", "rho", "p_value"]); w.writerows(sp_rows)

    # ---- 蛋白内配对: B vs C 同蛋白 ----
    by_prot_bc = defaultdict(lambda: ([], []))
    for s in gsites(grp_B - grp_C):
        by_prot_bc[s.rsplit("_", 1)[0]][0].append(s)
    for s in gsites(grp_C - grp_B):
        by_prot_bc[s.rsplit("_", 1)[0]][1].append(s)
    paired_prots = [a for a, (x, y) in by_prot_bc.items() if x and y]
    pr_rows = []
    for feat in FEATURES:
        diffs = []
        for a in paired_prots:
            x, y = by_prot_bc[a]
            diffs.append(sum(feats[s][feat] for s in x) / len(x) - sum(feats[s][feat] for s in y) / len(y))
        if len(paired_prots) >= MIN_PAIRED:
            wstat, p = wilcoxon_signed(diffs)
            pr_rows.append([feat, len(paired_prots), round(median(diffs), 4),
                            round(wstat, 1) if wstat is not None else "", f"{p:.3e}" if p is not None else ""])
        else:
            pr_rows.append([feat, len(paired_prots), round(median(diffs), 4) if diffs else "", "", ""])
    with open(os.path.join(a7, "paired_within_protein.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["feature", "n_paired_proteins", "median_paired_diff_BminusC", "W", "p"]); w.writerows(pr_rows)

    # ---- IDR Fisher ----
    idr_rows = []
    for gname, grp in (("A", grp_A), ("B", grp_B)):
        s1, s2 = gsites(grp), gsites(grp_C - grp)
        a = sum(1 for s in s1 if feats[s][IDR_CH] > 0.5); b = len(s1) - a
        c = sum(1 for s in s2 if feats[s][IDR_CH] > 0.5); d = len(s2) - c
        odds, p = fisher(a, b, c, d) if min(len(s1), len(s2)) else ("", "")
        idr_rows.append([gname, "C", round(a / len(s1), 4), round(c / len(s2), 4),
                         round(odds, 3) if odds != "" else "", f"{p:.3e}" if p != "" else ""])
    with open(os.path.join(a7, "idr_enrichment.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["group1", "group2", "idr_fraction_g1", "idr_fraction_g2", "odds_ratio", "fisher_p"])
        w.writerows(idr_rows)

    # ---- group QC ----
    gq = []
    comp_letters = {}
    for gname, grp in (("A", grp_A), ("B", grp_B), ("C", grp_C), ("ONI", grp_ONI), ("INO", grp_INO)):
        ss = gsites(grp)
        pp = defaultdict(int)
        letters = defaultdict(int)
        nb = defaultdict(int)
        for s in ss:
            acc, tail = s.rsplit("_", 1)
            pp[acc] += 1
            letters[tail[0]] += 1
            pos = int(tail[1:])
            if acc in seqs:
                sq = seqs[acc]
                if pos - 2 >= 0:
                    nb[sq[pos - 2]] += 1
                if pos < len(sq):
                    nb[sq[pos]] += 1
        cnts = sorted(pp.values())
        gq.append([gname, "n_sites", len(ss)])
        gq.append([gname, "n_proteins", len(pp)])
        gq.append([gname, "sites_per_protein_median", median(cnts) if cnts else ""])
        gq.append([gname, "sites_per_protein_max", max(cnts) if cnts else ""])
        tot = sum(letters.values()) or 1
        gq.append([gname, "STY_composition", " ".join(f"{k}:{v/tot*100:.1f}%" for k, v in sorted(letters.items()))])
        totn = sum(nb.values()) or 1
        top_nb = sorted(nb.items(), key=lambda x: -x[1])[:6]
        gq.append([gname, "neighbor_top6", " ".join(f"{k}:{v/totn*100:.1f}%" for k, v in top_nb)])
        comp_letters[gname] = {k: v / tot for k, v in letters.items()}
    # 窗口重叠 (B 组内, 同蛋白距离<=10)
    pos_by_prot = defaultdict(list)
    for s in gsites(grp_B):
        acc, tail = s.rsplit("_", 1)
        pos_by_prot[acc].append(int(tail[1:]))
    n_overlap = 0
    for acc, ps in pos_by_prot.items():
        ps.sort()
        for i in range(len(ps)):
            if (i > 0 and ps[i] - ps[i - 1] <= 10) or (i + 1 < len(ps) and ps[i + 1] - ps[i] <= 10):
                n_overlap += 1
    gq.append(["B", "window_overlap_sites(<=10aa)", n_overlap])
    with open(os.path.join(a7, "group_qc.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["group", "metric", "value"]); w.writerows(gq)

    # ---- summary ----
    top_by_rbc.sort(key=lambda x: -x[0])
    lines = [f"映射成功率: {len(feats)}/{len(all_sites)} ({100 - fail_rate:.1f}%)",
             f"分组规模(有特征): A={len(gsites(grp_A))} B={len(gsites(grp_B))} C={len(gsites(grp_C))} ONI={len(gsites(grp_ONI))} INO={len(gsites(grp_INO))}",
             ""]
    for cname, (nsig, nsig_perm, ndual) in summary_sig.items():
        lines.append(f"{cname}: BH<0.05 特征数 {nsig}/16 | 蛋白级置换 p<0.05 特征数 {nsig_perm}/16 | 置换中剔除双属蛋白 {ndual}")
        if nsig > 0 and nsig_perm == 0:
            lines.append(f"  注: {cname} 存在 MWU 显著但置换不显著的特征, 该显著性可能受蛋白内聚集影响")
        if nsig == 0:
            lines.append(f"  {cname}: 16 个特征在 BH 校正后均未达 0.05 显著水平")
    lines.append("")
    for absr, feat, rbc, v1, v2 in top_by_rbc[:3]:
        q = lambda v, p: sorted(v)[int((len(v) - 1) * p)]
        lines.append(f"A_vs_C 效应量 top: {feat} rbc={rbc:+.3f} "
                     f"(A 中位数{'高于' if rbc > 0 else '低于'}C) "
                     f"A[Q25/50/75]={q(v1,.25):.3f}/{q(v1,.5):.3f}/{q(v1,.75):.3f} "
                     f"C={q(v2,.25):.3f}/{q(v2,.5):.3f}/{q(v2,.75):.3f}")
    for row in idr_rows:
        lines.append(f"IDR(disoMine>0.5) 比例: 组{row[0]}={row[2]} vs C={row[3]}, OR={row[4]}, Fisher p={row[5]}")
    lines.append(f"Spearman 连续版: 32 组合中 p<0.05/32 的 {n_sp_sig} 个 (明细见 spearman_continuous.tsv)")
    lines.append(f"蛋白内配对 (B vs C 同蛋白): 可配对蛋白 {len(paired_prots)} 个"
                 + ("" if len(paired_prots) >= MIN_PAIRED else f" (< {MIN_PAIRED}, 仅报计数不做检验)"))
    a_comp = comp_letters.get("A", {})
    c_comp = comp_letters.get("C", {})
    drift = max((abs(a_comp.get(k, 0) - c_comp.get(k, 0)) for k in "STY"), default=0)
    lines.append(f"残基组成: A 与 C 的 S/T/Y 占比最大差 {drift*100:.1f} 个百分点 (明细见 group_qc.tsv)")
    txt = "\n".join(lines)
    for wbad in FORBIDDEN:
        assert wbad not in txt, f"禁用词出现: {wbad}"
    with open(os.path.join(a7, "feature_summary.txt"), "w") as f:
        f.write(txt + "\n")
    log("\n" + txt)

if __name__ == "__main__":
    main()
