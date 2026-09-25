#!/usr/bin/env python3
"""15_fpr_null.py — FPR 零试验: CTRL 4 重复劈 2v2 (三种劈法), 跑与 07/08 相同的筛选规则。
真值为零 => 任何通过筛选的位点均为经验假阳性。
镜像 07 的规则: occupancy 响应 = |median log2FC|>=1 且每侧>=2重复且方向一致 (分母档位: 残基优先, 不足半数回退蛋白级);
              轮廓位移 = 6组分 CLR (伪计数=0.5*min正值, 组分内取中位) 的 Aitchison 距离
主判定阈值 = 真实分析的 P90 位移阈值 (默认 8.4018, 可作参数覆盖); 同时报告各劈法内部 P90 供参照。
输出: analysis07/fpr_null_summary.txt; 屏幕同步打印。
运行: login 直接跑, ~2-3 分钟 (48 个 CTRL parquet):
  module load Python/3.13.1-GCCcore-14.2.0 && python3 15_fpr_null.py
"""
import sys, os, re, csv, glob, math
from collections import defaultdict

LOC, QV, FC_THRESH = 0.75, 0.01, 1.0
FRS = [f"FR{i}" for i in range(1, 7)]
SPLITS = [(("Rep1","Rep2"),("Rep3","Rep4")), (("Rep1","Rep3"),("Rep2","Rep4")), (("Rep1","Rep4"),("Rep2","Rep3"))]
REAL_THR = float(sys.argv[1]) if len(sys.argv) > 1 else 8.4018
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

def median(v):
    s = sorted(v); n = len(s)
    return s[n//2] if n % 2 else (s[n//2-1] + s[n//2]) / 2

def clr(props):
    logs = [math.log(p) for p in props]
    m = sum(logs) / len(logs)
    return [x - m for x in logs]

def dist(a, b):
    return math.sqrt(sum((x-y)**2 for x, y in zip(a, b)))

def pctl(sv, p):
    if not sv: return None
    k = (len(sv)-1)*p/100; f = int(k)
    return sv[f] if f == k else sv[f] + (sv[f+1]-sv[f])*(k-f)

def scan(dirs, seqs, need_conf, sites_by_acc=None):
    """返回 num[(id,fr,rep)]->qty 与 run_logs[(fr,rep)]->[log2 qty]; phospho 时 id=sid, proteome 时同时填残基/蛋白分母"""
    import pyarrow.parquet as pq
    find_cache = {}
    num = defaultdict(float); den_res = defaultdict(float); den_prot = defaultdict(float)
    run_logs = defaultdict(list); site_info = {}
    for d in dirs:
        run = os.path.basename(d.rstrip("/")); key = parse_design(run)
        if key is None or key[0] != "CTRL": continue
        _, fr, rep = key
        p = os.path.join(d, "report.parquet")
        if not os.path.isfile(p): continue
        pf = pq.ParquetFile(p); cols = pf.schema_arrow.names
        want = ["Modified.Sequence","Stripped.Sequence","Protein.Ids","Precursor.Quantity","Q.Value"] + (["PTM.Site.Confidence"] if need_conf else [])
        use = [c for c in want if c in cols]
        for b in pf.iter_batches(columns=use, batch_size=100_000):
            cd = {n: b.column(n).to_pylist() for n in use}
            for i in range(b.num_rows):
                if (cd["Q.Value"][i] or 1) > QV: continue
                qty = cd["Precursor.Quantity"][i]
                if not qty or qty <= 0: continue
                pep = cd["Stripped.Sequence"][i] or ""
                accs = [a for a in re.split(r"[;,]", cd["Protein.Ids"][i] or "") if a]
                if need_conf:
                    mod = cd["Modified.Sequence"][i] or ""
                    if "UniMod:21" not in mod: continue
                    conf = cd["PTM.Site.Confidence"][i]
                    if conf is None or conf < LOC: continue
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
                    run_logs[(fr, rep)].append(math.log2(qty))
                    for pp, r in ps:
                        sid = f"{acc}_{r}{st+pp}"
                        site_info[sid] = acc
                        num[(sid, fr, rep)] += qty
                else:
                    run_logs[(fr, rep)].append(math.log2(qty))
                    for acc in accs:
                        if acc not in sites_by_acc: continue
                        s = seqs.get(acc)
                        if s is None: continue
                        ck = (acc, pep)
                        if ck not in find_cache: find_cache[ck] = s.find(pep)
                        st = find_cache[ck]
                        if st < 0: continue
                        den_prot[(acc, fr, rep)] += qty
                        end = st + len(pep)
                        for pos, sid in sites_by_acc[acc]:
                            if st < pos <= end: den_res[(sid, fr, rep)] += qty
                        break
        print(f"  [scan] {run}", flush=True)
    return num, den_res, den_prot, run_logs, site_info

def main():
    base = os.path.join(os.environ["VSC_DATA"], "rerun_2026-08")
    a7 = os.path.join(base, "occupancy/analysis07")
    fasta = sorted(glob.glob(os.path.join(os.environ["VSC_SCRATCH"], "rerun/fasta/human_*_UP*.fasta")))[-1]
    seqs = load_fasta(fasta)
    num, _, _, rl_p, site_info = scan(sorted(glob.glob(os.path.join(base, "phospho", "*CTRL*/"))), seqs, True)
    sites_by_acc = defaultdict(list)
    for sid, acc in site_info.items():
        sites_by_acc[acc].append((int(sid.rsplit("_",1)[1][1:]), sid))
    _, den_res, den_prot, rl_q, _ = scan(sorted(glob.glob(os.path.join(base, "proteome", "*CTRL*/"))), seqs, False, sites_by_acc)

    def factors(rl):
        meds = {k: median(v) for k, v in rl.items() if v}
        anchor = median(list(meds.values()))
        return {k: 2**(anchor - m) for k, m in meds.items()}
    fp, fq = factors(rl_p), factors(rl_q)
    numN = {k: v*fp.get((k[1],k[2]),1.0) for k, v in num.items()}
    denR = {k: v*fq.get((k[1],k[2]),1.0) for k, v in den_res.items()}
    denP = {k: v*fq.get((k[1],k[2]),1.0) for k, v in den_prot.items()}

    pos_all = sorted(numN.values()); pseudo = 0.5*pos_all[0]
    lines = [f"FPR 零试验 (CTRL 2v2, 三种劈法) | 真实分析位移阈值 = {REAL_THR}",
             f"CTRL phospho 位点数: {len(site_info)}"]
    agg = []
    for si, (h1, h2) in enumerate(SPLITS, 1):
        # occupancy 响应 (镜像 07 responder, 半1=基线 半2=伪处理)
        traj = defaultdict(lambda: defaultdict(dict))   # (sid,fr) -> half -> {rep: occ}
        for (sid, fr, rep), p in numN.items():
            half = "H1" if rep in h1 else "H2"
            acc = site_info[sid]
            obs_all = [(r,) for r in (h1+h2) if (sid, fr, r) in numN]
            n_res = sum(1 for (r,) in obs_all if (sid, fr, r) in denR)
            tier_res = n_res >= 0.5*len(obs_all)
            d = denR.get((sid, fr, rep)) if tier_res else denP.get((acc, fr, rep))
            if d: traj[(sid, fr)][half][rep] = math.log2(p/d)
        occ_q = set(); occ_resp = set()
        for (sid, fr), hv in traj.items():
            v1 = list(hv.get("H1", {}).values()); v2 = list(hv.get("H2", {}).values())
            if len(v1) < 2 or len(v2) < 2: continue
            occ_q.add(sid)
            base_m = median(v1); fc = median(v2) - base_m
            if abs(fc) >= FC_THRESH and all((x-base_m)*fc > 0 for x in v2):
                occ_resp.add(sid)
        # CLR 位移 (半1 vs 半2, 各半组分中位)
        prof = defaultdict(lambda: {"H1": defaultdict(list), "H2": defaultdict(list)})
        for (sid, fr, rep), v in numN.items():
            prof[sid]["H1" if rep in h1 else "H2"][fr].append(v)
        shifts = {}
        for sid, hv in prof.items():
            if not hv["H1"] or not hv["H2"]: continue
            if sum(len(x) for x in hv["H1"].values()) < 2 or sum(len(x) for x in hv["H2"].values()) < 2: continue
            def cp(frv):
                vec = [median(frv[fr]) if frv.get(fr) else 0.0 for fr in FRS]
                vec = [v+pseudo for v in vec]; s = sum(vec)
                return clr([v/s for v in vec])
            shifts[sid] = dist(cp(hv["H1"]), cp(hv["H2"]))
        sv = sorted(shifts.values())
        own_p90 = pctl(sv, 90)
        movers = {s for s, v in shifts.items() if v >= REAL_THR}
        joint = occ_resp & movers
        occ_q_mov = occ_q & movers
        agg.append((len(occ_q), len(occ_resp), len(shifts), len(movers), len(occ_q_mov), len(joint), own_p90))
        lines.append(f"劈法{si} {h1}v{h2}: occ可算 {len(occ_q)} | occ假响应 {len(occ_resp)} ({len(occ_resp)/max(1,len(occ_q))*100:.2f}%) | "
                     f"轮廓可算 {len(shifts)} | 位移>=真实阈值 {len(movers)} ({len(movers)/max(1,len(shifts))*100:.2f}%) | "
                     f"occ可算∩位移达标 {len(occ_q_mov)} ({len(occ_q_mov)/max(1,len(occ_q))*100:.2f}%, 真实分析对应值 17.3%) | "
                     f"双重规则同过 {len(joint)} | 本劈法内部P90 {own_p90:.4f}")
    m = [sum(a[i] for a in agg)/3 for i in range(7)]
    lines.append(f"三劈法均值: occ假响应率 {m[1]/max(1,m[0])*100:.2f}% | 位移超真实阈值率 {m[3]/max(1,m[2])*100:.2f}% | occ∩位移 {m[4]/max(1,m[0])*100:.2f}% | 双重规则 {m[5]:.1f} 个")
    lines.append("判读边界: 2v2 每侧重复数少于真实分析的 4v4 结构, 噪声更大, 此估计偏保守(偏高); 比较基准为占比而非绝对数")
    txt = "\n".join(lines)
    for w in FORBIDDEN: assert w not in txt, w
    with open(os.path.join(a7, "fpr_null_summary.txt"), "w") as f: f.write(txt + "\n")
    print("\n" + txt)

if __name__ == "__main__":
    main()
