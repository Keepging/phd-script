#!/usr/bin/env python3
"""07_main_analysis.py — 三层主分析。

层1 occupancy: 位点x组分轨迹, 分母档位按轨迹统一(残基级优先, 不足半数回退蛋白级);
               各 EGF 时间点 vs CTRL 的 log2FC (occupancy 与 phospho 强度并行);
               解离分类: both / intensity_only / occupancy_only(核心发现类) / none
层2 轮廓(分母无关): 位点磷酸化形式的 6 组分 CLR 轮廓, 相对 CTRL 的 Aitchison 位移;
               蛋白总量轮廓(proteome); 位点轮廓 vs 母蛋白轮廓的偏离度及其时间变化
层3 自检: 剔除 90min FR1/FR3 Rep3 两个离群 run 后解离分类计数的变化

强度均经 per-run log2 中位数中心化 (两条线各自), 消除进样量差异。
筛选标准: |log2FC|>=1, 每侧 >=2 个重复, 且方向一致 —— 明确为候选筛选, 非最终统计检验。

用法: python3 07_main_analysis.py [结果根目录] [FASTA路径]  (默认同 05/06)
输出: <根目录>/occupancy/analysis07/{site_traj.tsv, site_profiles.tsv, dissociation_summary.txt, summary07.txt}
"""
import sys, os, re, csv, glob, math
from collections import defaultdict

LOC_CUTOFF = 0.75
QVAL_CUTOFF = 0.01
FC_THRESH = 1.0
MIN_REPS = 2
TPS = ["CTRL", "2min", "8min", "20min", "90min"]
EGF_TPS = ["2min", "8min", "20min", "90min"]
FRS = [f"FR{i}" for i in range(1, 7)]
OUTLIER_OBS = {("90min", "FR1", "Rep3"), ("90min", "FR3", "Rep3")}

def log(m):
    print(m, flush=True)

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

def median(v):
    return sorted(v)[len(v) // 2] if len(v) % 2 else sum(sorted(v)[len(v) // 2 - 1:len(v) // 2 + 1]) / 2

def clr(props):
    logs = [math.log(p) for p in props]
    m = sum(logs) / len(logs)
    return [x - m for x in logs]

def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))

def main():
    base = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ["VSC_DATA"], "rerun_2026-08")
    if len(sys.argv) > 2:
        fasta_path = sys.argv[2]
    else:
        cands = sorted(glob.glob(os.path.join(os.environ["VSC_SCRATCH"], "rerun/fasta/human_*_UP*.fasta")))
        fasta_path = cands[-1] if cands else ""
    if not os.path.isfile(fasta_path):
        log(f"FAIL: FASTA 不存在: {fasta_path}"); sys.exit(1)
    outdir = os.path.join(base, "occupancy", "analysis07")
    os.makedirs(outdir, exist_ok=True)
    seqs = load_fasta(fasta_path)
    phos_dirs = sorted(glob.glob(os.path.join(base, "phospho", "*/")))
    prot_dirs = sorted(glob.glob(os.path.join(base, "proteome", "*/")))
    log(f"FASTA {len(seqs)} | phospho {len(phos_dirs)} | proteome {len(prot_dirs)}")
    if not phos_dirs or not prot_dirs:
        log("FAIL: run 目录缺失"); sys.exit(1)

    find_cache, site_info = {}, {}
    num_raw = defaultdict(float)          # (sid, key) -> 原始分子
    run_logs_p = defaultdict(list)        # phospho run key -> log2 行强度 (求 run 中位)

    # ---- 扫 phospho ----
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
                run_logs_p[key].append(math.log2(qty))
                for ppos, res in psites:
                    sid = f"{acc}_{res}{start + ppos}"
                    site_info[sid] = (acc, start + ppos, res)
                    num_raw[(sid, key)] += qty
        log(f"  [phos] {run}")
    if not site_info:
        log("FAIL: 零位点"); sys.exit(1)

    sites_by_acc = defaultdict(list)
    for sid, (acc, pos, _r) in site_info.items():
        sites_by_acc[acc].append((pos, sid))

    # ---- 扫 proteome: 残基级/蛋白级分母 + run 中位 ----
    den_res_raw = defaultdict(float)
    den_prot_raw = defaultdict(float)     # (acc, key)
    run_logs_q = defaultdict(list)
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
                run_logs_q[key].append(math.log2(qty))
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
                    den_prot_raw[(acc, key)] += qty
                    end = start + len(pep)
                    for pos, sid in sites_by_acc[acc]:
                        if start < pos <= end:
                            den_res_raw[(sid, key)] += qty
                    break
        log(f"  [prot] {run}")

    # ---- per-run 中位数中心化因子 (log2), 以各线 run 中位数的中位数为锚 ----
    def factors(run_logs):
        meds = {k: median(v) for k, v in run_logs.items() if v}
        anchor = median(list(meds.values()))
        return {k: 2 ** (anchor - m) for k, m in meds.items()}
    fp = factors(run_logs_p)
    fq = factors(run_logs_q)
    num = {k: v * fp.get(k[1], 1.0) for k, v in num_raw.items()}
    den_res = {k: v * fq.get(k[1], 1.0) for k, v in den_res_raw.items()}
    den_prot = {k: v * fq.get(k[1], 1.0) for k, v in den_prot_raw.items()}

    # ---- 层1: 轨迹 + 分母档位统一 + 解离分类 ----
    traj = defaultdict(list)               # (sid, fr) -> [(tp, rep)]
    for (sid, (tp, fr, rep)) in num:
        traj[(sid, fr)].append((tp, rep))

    def classify(exclude_outliers):
        rows, counts = [], defaultdict(int)
        for (sid, fr), obs in sorted(traj.items()):
            acc = site_info[sid][0]
            obs_use = [(tp, rep) for tp, rep in obs
                       if not (exclude_outliers and (tp, fr, rep) in OUTLIER_OBS)]
            if not obs_use:
                continue
            n_res = sum(1 for tp, rep in obs_use if (sid, (tp, fr, rep)) in den_res)
            tier = "residue" if n_res >= 0.5 * len(obs_use) else "protein"
            dd = den_res if tier == "residue" else {}
            occ, inten = defaultdict(list), defaultdict(list)
            for tp, rep in obs_use:
                key = (tp, fr, rep)
                p = num[(sid, key)]
                inten[tp].append(math.log2(p))
                d = den_res.get((sid, key)) if tier == "residue" else den_prot.get((acc, key))
                if d:
                    occ[tp].append(math.log2(p / d))

            def responder(vals):
                if len(vals.get("CTRL", [])) < MIN_REPS:
                    return False, ""
                base = median(vals["CTRL"])
                for tp in EGF_TPS:
                    v = vals.get(tp, [])
                    if len(v) < MIN_REPS:
                        continue
                    fc = median(v) - base
                    if abs(fc) >= FC_THRESH and all((x - base) * fc > 0 for x in v):
                        return True, tp
                return False, ""
            o_resp, o_tp = responder(occ)
            i_resp, i_tp = responder(inten)
            cls = ("both" if o_resp and i_resp else
                   "occupancy_only" if o_resp else
                   "intensity_only" if i_resp else "none")
            counts[cls] += 1
            counts["tier_" + tier] += 1
            if not exclude_outliers:
                occ_m = [round(median(occ[tp]), 3) if len(occ.get(tp, [])) >= MIN_REPS else "" for tp in TPS]
                int_m = [round(median(inten[tp]), 3) if len(inten.get(tp, [])) >= MIN_REPS else "" for tp in TPS]
                rows.append([sid, acc, fr, tier, cls, o_tp or i_tp] + occ_m + int_m)
        return rows, counts

    rows_full, counts_full = classify(False)
    _, counts_sens = classify(True)

    with open(os.path.join(outdir, "site_traj.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["site_id", "protein", "fraction", "denom_tier", "class", "peak_tp"]
                   + [f"log2occ_{t}" for t in TPS] + [f"log2int_{t}" for t in TPS])
        w.writerows(rows_full)

    # ---- 层2: CLR 轮廓 ----
    # 位点轮廓: (sid, tp) -> {fr: [强度]}
    prof_site = defaultdict(lambda: defaultdict(list))
    for (sid, (tp, fr, rep)), v in num.items():
        prof_site[(sid, tp)][fr].append(v)
    prof_prot = defaultdict(lambda: defaultdict(list))
    for (acc, (tp, fr, rep)), v in den_prot.items():
        prof_prot[(acc, tp)][fr].append(v)

    def clr_profile(frvals, pseudo):
        vec = [median(frvals[fr]) if frvals.get(fr) else 0.0 for fr in FRS]
        vec = [v + pseudo for v in vec]
        s = sum(vec)
        return clr([v / s for v in vec])

    pos_all = [v for m in prof_site.values() for lst in m.values() for v in lst]
    pseudo = 0.5 * min(pos_all)
    pos_q = [v for m in prof_prot.values() for lst in m.values() for v in lst]
    pseudo_q = 0.5 * min(pos_q) if pos_q else pseudo

    def n_obs(frvals):
        return sum(len(v) for v in frvals.values())

    with open(os.path.join(outdir, "site_profiles.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["site_id", "protein", "proteome_tier", "timepoint",
                    "phos_shift_vs_ctrl", "prot_shift_vs_ctrl",
                    "site_vs_prot_divergence", "divergence_delta_vs_ctrl"])
        n_shift = 0
        for sid in sorted(site_info):
            acc = site_info[sid][0]
            has_prot = any((acc, tp) in prof_prot for tp in TPS)
            tier3 = not has_prot
            ctrl_s = prof_site.get((sid, "CTRL"))
            if not ctrl_s or n_obs(ctrl_s) < MIN_REPS:
                continue
            clr_ctrl = clr_profile(ctrl_s, pseudo)
            clr_pc = clr_profile(prof_prot[(acc, "CTRL")], pseudo_q) if (acc, "CTRL") in prof_prot else None
            div_ctrl = dist(clr_ctrl, clr_pc) if clr_pc else None
            for tp in EGF_TPS:
                s = prof_site.get((sid, tp))
                if not s or n_obs(s) < MIN_REPS:
                    continue
                clr_tp = clr_profile(s, pseudo)
                shift = dist(clr_tp, clr_ctrl)
                n_shift += 1
                pshift = div = ddelta = ""
                if clr_pc and (acc, tp) in prof_prot:
                    clr_pt = clr_profile(prof_prot[(acc, tp)], pseudo_q)
                    pshift = round(dist(clr_pt, clr_pc), 4)
                    div = round(dist(clr_tp, clr_pt), 4)
                    ddelta = round(dist(clr_tp, clr_pt) - div_ctrl, 4)
                w.writerow([sid, acc, "T3" if tier3 else "T1/T2", tp,
                            round(shift, 4), pshift, div, ddelta])

    # ---- 汇总 ----
    lines = ["== 层1 解离分类 (位点x组分轨迹) =="]
    total = sum(counts_full[c] for c in ("both", "occupancy_only", "intensity_only", "none"))
    for c in ("both", "occupancy_only", "intensity_only", "none"):
        lines.append(f"  {c}: {counts_full[c]} ({counts_full[c]/total*100:.1f}%)")
    lines.append(f"  分母档位: residue {counts_full['tier_residue']} | protein {counts_full['tier_protein']}")
    lines.append("== 层3 敏感性 (剔除 90min FR1/FR3 Rep3) ==")
    for c in ("both", "occupancy_only", "intensity_only"):
        lines.append(f"  {c}: {counts_full[c]} -> {counts_sens[c]}")
    lines.append("== 层2 轮廓 ==")
    lines.append(f"  可算轮廓位移的 (位点x时间点): {n_shift}")
    txt = "\n".join(lines)
    with open(os.path.join(outdir, "summary07.txt"), "w") as f:
        f.write(txt + "\n")
    log("\n" + txt)
    log(f"输出目录: {outdir}")

if __name__ == "__main__":
    main()
