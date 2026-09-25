#!/usr/bin/env python3
"""14_breadth_regression.py — C 补充: 强度校正回归 (敏感性第三口径)。
每特征三个模型: M0 feature~breadth | M1 feature~breadth+log2int | M2 feature~G2,G3哑变量+log2int
蛋白级聚类稳健标准误 (同蛋白多位点不独立); M2 对两哑变量做联合 Wald (df=2)。
BH 校正 M1 的 breadth p (16 特征)。
输入: analysis07/breadth_sites.tsv | 输出: analysis07/breadth_regression.tsv + 摘要追加到 breadth_summary.txt
运行: login 节点直接跑, 秒级 —
  module load Python/3.13.1-GCCcore-14.2.0 && python3 14_breadth_regression.py
判读边界(印在输出): breadth 检出本身由强度驱动, 二者部分同源; 显著仅表示"给定强度下的边际关联", 不表示独立空间效应。
"""
import os, csv, math
import numpy as np

FORBIDDEN = ["证明","验证了","breakthrough","突破","隐藏的信号","hidden signal"]

def cluster_ols(X, y, clusters):
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta
    meat = np.zeros((X.shape[1], X.shape[1]))
    for g in np.unique(clusters):
        m = clusters == g
        Xg, eg = X[m], resid[m]
        v = Xg.T @ eg
        meat += np.outer(v, v)
    G = len(np.unique(clusters))
    cov = XtX_inv @ meat @ XtX_inv * (G / (G - 1))
    return beta, cov

def p_z(z): return math.erfc(abs(z) / math.sqrt(2))

def bh(ps):
    n = len(ps); order = sorted(range(n), key=lambda i: ps[i]); out = [0.0]*n; prev = 1.0
    for k, i in enumerate(reversed(order)):
        prev = min(prev, ps[i]*n/(n-k)); out[i] = prev
    return out

def main():
    a7 = os.path.join(os.environ["VSC_DATA"], "rerun_2026-08/occupancy/analysis07")
    rows = list(csv.DictReader(open(os.path.join(a7, "breadth_sites.tsv")), delimiter="\t"))
    feats = [c for c in rows[0] if c.startswith(("single_","win5_"))]
    keep = [r for r in rows if r["log2_intensity_ctrl"] not in ("", "nan")]
    print(f"位点: 总 {len(rows)} | 有强度 {len(keep)}")
    br = np.array([float(r["breadth"]) for r in keep])
    li = np.array([float(r["log2_intensity_ctrl"]) for r in keep])
    prot = np.array([r["protein"] for r in keep])
    g2 = np.array([1.0 if r["group"] == "G2_moderate" else 0.0 for r in keep])
    g3 = np.array([1.0 if r["group"] == "G3_ubiquitous" else 0.0 for r in keep])
    one = np.ones(len(keep))
    zli = (li - li.mean()) / li.std()
    zbr = (br - br.mean()) / br.std()

    out, p1s = [], []
    for ft in feats:
        y = np.array([float(r[ft]) for r in keep])
        zy = (y - y.mean()) / (y.std() if y.std() > 0 else 1)
        # M0: y ~ breadth
        b0, c0 = cluster_ols(np.column_stack([one, zbr]), zy, prot)
        p0 = p_z(b0[1] / math.sqrt(c0[1, 1]))
        # M1: y ~ breadth + intensity
        b1, c1 = cluster_ols(np.column_stack([one, zbr, zli]), zy, prot)
        z1 = b1[1] / math.sqrt(c1[1, 1]); p1 = p_z(z1)
        # M2: y ~ G2 + G3 + intensity, 联合 Wald df=2
        X2 = np.column_stack([one, g2, g3, zli])
        b2, c2 = cluster_ols(X2, zy, prot)
        sub, csub = b2[1:3], c2[1:3, 1:3]
        W = float(sub @ np.linalg.pinv(csub) @ sub)
        p2 = math.exp(-W / 2)
        out.append([ft, round(b0[1], 4), f"{p0:.3e}", round(b1[1], 4), f"{p1:.3e}",
                    round(b1[2], 4), round(b2[1], 4), round(b2[2], 4), f"{p2:.3e}"])
        p1s.append(p1)
    pbh = bh(p1s)
    with open(os.path.join(a7, "breadth_regression.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["# 判读边界: breadth 的检出由强度驱动, 二者部分同源; M1/M2 显著仅为给定强度下的边际关联, 非独立空间效应"])
        w.writerow(["feature","beta_breadth_M0","p_M0","beta_breadth_M1_adj","p_M1_adj",
                    "beta_intensity_M1","beta_G2_M2","beta_G3_M2","p_joint_M2","p_M1_BH"])
        for r, pb in zip(out, pbh): w.writerow(r + [f"{pb:.3e}"])
    n0 = sum(1 for r in out if float(r[2]) < 0.05)
    n1 = sum(1 for pb in pbh if pb < 0.05)
    n2 = sum(1 for r in out if float(r[8]) < 0.05)
    att = [abs(float(r[3])) / max(abs(float(r[1])), 1e-9) for r in out]
    med_att = sorted(att)[len(att)//2]
    S = (f"回归口径: M0 未校正显著 {n0}/16 | M1 强度校正后 BH<0.05: {n1}/16 | M2 因子联合 Wald<0.05: {n2}/16 | "
         f"校正后 |beta| 中位保留比例 {med_att:.2f} | 最大校正后 |beta_breadth| = "
         f"{max(abs(float(r[3])) for r in out):.4f}")
    for wbad in FORBIDDEN: assert wbad not in S, wbad
    with open(os.path.join(a7, "breadth_summary.txt"), "a") as f: f.write(S + "\n")
    print(S)
    print("输出: breadth_regression.tsv")

if __name__ == "__main__":
    main()
