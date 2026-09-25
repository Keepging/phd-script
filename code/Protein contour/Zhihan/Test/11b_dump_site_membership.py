#!/usr/bin/env python3
"""11b_dump_site_membership.py — Task A2 变体: 导出位点级明细清单 (site membership)。
基于 11_martinez_counts.py 改写: 过滤逻辑/常量/FASTA选择/site_id 构造 100% 不变 (直接复制)。
新增: 把每个检出位点 (口径W, 即 det 中的每个 sid) 展开成一行, 记录该位点的
      accession、overall aligned (口径M, m_pass over all cells)、以及每个
      fraction 内的 aligned_FRi (m_pass 限定该 fraction)。
输入: $VSC_DATA/rerun_2026-08/phospho/*/report.parquet + $VSC_SCRATCH/rerun/fasta/human_*_UP*.fasta
输出: $VSC_DATA/rerun_2026-08/occupancy/analysis07/site_membership.tsv
      (同时保留原 rerun_counts.tsv 计数输出, 不改变原脚本行为)
预计: 与 11_martinez_counts.py 同等运行时间 (~5-10 分钟)
"""
import os, re, csv, glob, sys
from collections import defaultdict

LOC, QV = 0.75, 0.01
FRS = [f"FR{i}" for i in range(1, 7)]

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

def main():
    import pyarrow.parquet as pq
    base = os.path.join(os.environ["VSC_DATA"], "rerun_2026-08")
    fasta = sorted(glob.glob(os.path.join(os.environ["VSC_SCRATCH"], "rerun/fasta/human_*_UP*.fasta")))[-1]
    seqs = load_fasta(fasta)
    find_cache = {}
    det = defaultdict(set)   # (sid) -> set of (tp, fr, rep)
    site_acc = {}
    dirs = sorted(glob.glob(os.path.join(base, "phospho", "*/")))
    for d in dirs:
        run = os.path.basename(d.rstrip("/")); key = parse_design(run)
        if key is None: continue
        p = os.path.join(d, "report.parquet")
        if not os.path.isfile(p): print(f"WARN 缺 {run}", flush=True); continue
        pf = pq.ParquetFile(p); cols = pf.schema_arrow.names
        use = [c for c in ("Modified.Sequence","Stripped.Sequence","Protein.Ids","Q.Value","PTM.Site.Confidence") if c in cols]
        for b in pf.iter_batches(columns=use, batch_size=100_000):
            cd = {n: b.column(n).to_pylist() for n in use}
            for i in range(b.num_rows):
                if (cd["Q.Value"][i] or 1) > QV: continue
                mod = cd["Modified.Sequence"][i] or ""
                if "UniMod:21" not in mod: continue
                conf = cd["PTM.Site.Confidence"][i]
                if conf is None or conf < LOC: continue
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
                    det[sid].add(key)
        print(f"  [scan] {run}", flush=True)

    def m_pass(cells, fr=None):
        cnt = defaultdict(int)
        for tp, f, rep in cells:
            if fr is None or f == fr: cnt[(tp, f)] += 1
        return any(v >= 3 for v in cnt.values())

    rows = [["scope","metric","wide","martinez_aligned"]]
    sW = set(det); sM = {s for s, c in det.items() if m_pass(c)}
    rows.append(["overall","phospho_sites",len(sW),len(sM)])
    rows.append(["overall","phospho_proteins",len({site_acc[s] for s in sW}),len({site_acc[s] for s in sM})])
    for fr in FRS:
        fW = {s for s, c in det.items() if any(f == fr for _, f, _ in c)}
        fM = {s for s, c in det.items() if m_pass(c, fr)}
        rows.append([fr,"phospho_sites",len(fW),len(fM)])
        rows.append([fr,"phospho_proteins",len({site_acc[s] for s in fW}),len({site_acc[s] for s in fM})])
    out = os.path.join(base, "occupancy/analysis07/rerun_counts.tsv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["# 口径M假设: >=3/4 重复在任一 fraction x timepoint 格子; per-fraction 限定该 fraction 内任一时间点"])
        w.writerows(rows)
    for r in rows: print("\t".join(map(str, r)))
    print(f"输出: {out}")

    # ---- 新增: 位点级明细清单 (site membership) ----
    mem_out = os.path.join(base, "occupancy/analysis07/site_membership.tsv")
    mem_header = ["site_id", "accession", "aligned"] + [f"aligned_{fr}" for fr in FRS]
    with open(mem_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(mem_header)
        for sid in sorted(det):
            cells = det[sid]
            row = [sid, site_acc[sid], m_pass(cells)]
            row += [m_pass(cells, fr) for fr in FRS]
            w.writerow(row)
    print(f"输出: {mem_out}")

    # ---- 自校验 ----
    expected = {
        "wide_sites": 18268,
        "wide_proteins": 4996,
        "aligned_sites": 8120,
        "aligned_proteins": 3138,
    }
    actual = {
        "wide_sites": len(sW),
        "wide_proteins": len({site_acc[s] for s in sW}),
        "aligned_sites": len(sM),
        "aligned_proteins": len({site_acc[s] for s in sM}),
    }
    print("\n== 自校验 ==")
    all_ok = True
    for k in expected:
        exp, act = expected[k], actual[k]
        status = "OK" if exp == act else "MISMATCH"
        if status == "MISMATCH": all_ok = False
        print(f"{k}: expected={exp} actual={act} -> {status}")
    print("== 总体: " + ("OK" if all_ok else "MISMATCH") + " ==")

if __name__ == "__main__":
    main()
