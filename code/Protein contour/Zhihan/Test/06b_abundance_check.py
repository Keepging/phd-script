#!/usr/bin/env python3
"""06b_abundance_check.py — 用数据检验"严格子集是否偏向高丰度蛋白"。

读 occupancy_variants_long.tsv, 把蛋白分三档:
  T1 严格档: 至少一个样本里有残基级分母 (v1 算得出)
  T2 回退档: 从无残基级分母, 但有蛋白级分母 (只有 vA 算得出)
  T3 隐形档: 两种分母都从未出现 (proteome 里查无此蛋白)
对 T1/T2 输出蛋白丰度 (蛋白级分母的样本中位数, log10) 的分位数对比。
判读: T1 丰度显著高于 T2 => 偏倚证实; T3 计数 = 任何 occupancy 方案都覆盖不了的蛋白规模。

用法: python3 06b_abundance_check.py [occupancy_variants_long.tsv 路径]
"""
import sys, os, csv, math, statistics as st
from collections import defaultdict

path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.environ["VSC_DATA"], "rerun_2026-08/occupancy/occupancy_variants_long.tsv")

has_res = set()
prot_int = defaultdict(list)   # protein -> 蛋白级分母观测
site_prot = {}
with open(path) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        prot = r["site_id"].rsplit("_", 1)[0]
        site_prot[r["site_id"]] = prot
        if r["denom_residue"]:
            has_res.add(prot)
        if r["denom_protein"]:
            prot_int[prot].append(float(r["denom_protein"]))

all_prots = set(site_prot.values())
t1 = sorted(has_res)
t2 = sorted(p for p in prot_int if p not in has_res)
t3 = sorted(p for p in all_prots if p not in prot_int)

def q(vals, p):
    s = sorted(vals); k = (len(s) - 1) * p; f = int(k)
    return s[f] if f == k else s[f] + (s[f + 1] - s[f]) * (k - f)

def describe(tier, prots):
    med = [math.log10(st.median(prot_int[p])) for p in prots if prot_int.get(p)]
    if not med:
        print(f"{tier}: n={len(prots)} (无丰度可算)"); return None
    print(f"{tier}: 蛋白数 {len(prots)} | log10 丰度 中位 {st.median(med):.2f} | 四分位 [{q(med,0.25):.2f}, {q(med,0.75):.2f}]")
    return med

print(f"输入: {path}")
print(f"位点数 {len(site_prot)} | 涉及蛋白 {len(all_prots)}\n")
m1 = describe("T1 严格档(v1可算)", t1)
m2 = describe("T2 回退档(仅vA可算)", t2)
print(f"T3 隐形档(occupancy不可算): 蛋白数 {len(t3)} | 位点数 {sum(1 for s,p in site_prot.items() if p in set(t3))}")
if m1 and m2:
    diff = st.median(m1) - st.median(m2)
    print(f"\nT1 与 T2 的中位丰度差: {diff:.2f} 个 log10 单位 (= {10**diff:.1f} 倍)")
    print("判读: 差值 > 0.5 (约3倍) 即可认为丰度偏倚成立, 高丰度假说获数据支持")
