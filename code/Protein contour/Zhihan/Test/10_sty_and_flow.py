#!/usr/bin/env python3
"""10_sty_and_flow.py — Task2: 656 候选的 S/T/Y 拆解; Task3: FR_from→FR_to 流向矩阵。
输入: analysis07/{top_candidates.tsv, site_traj.tsv}
输出: analysis07/{sty_breakdown.txt, sty_class_contingency.tsv,
       flow_matrix_intensity.tsv, flow_matrix_occupancy.tsv,
       flow_heatmap_intensity.png, flow_heatmap_occupancy.png (matplotlib 可用时)}
标注: candidate set = 位移度量 v1。卡方 df=2 的 p 用精确式 exp(-x/2)。
用法: python3 10_sty_and_flow.py [analysis07目录]
"""
import sys, os, csv, math
from collections import defaultdict, Counter

TPS = ["CTRL", "2min", "8min", "20min", "90min"]
EGF = ["90min", "20min", "8min", "2min"]   # 从晚到早找"最晚可用"
FRS = [f"FR{i}" for i in range(1, 7)]

def main():
    a7 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ["VSC_DATA"], "rerun_2026-08/occupancy/analysis07")
    cand = {}
    with open(os.path.join(a7, "top_candidates.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cand[r["site_id"]] = r
    traj = defaultdict(dict)   # sid -> fr -> row
    all_sites = set()
    for r in csv.DictReader(open(os.path.join(a7, "site_traj.tsv")), delimiter="\t"):
        all_sites.add(r["site_id"])
        traj[r["site_id"]][r["fraction"]] = r

    out = []
    def emit(s):
        print(s, flush=True); out.append(s)

    # ================= TASK 2 =================
    res = lambda sid: sid.rsplit("_", 1)[1][0]
    c_cnt = Counter(res(s) for s in cand)
    a_cnt = Counter(res(s) for s in all_sites)
    n_c, n_a = len(cand), len(all_sites)
    emit("=== TASK 2: S/T/Y 拆解 (candidate set = 位移度量 v1) ===")
    emit(f"候选 656 集合: " + " | ".join(f"{k}={c_cnt.get(k,0)} ({c_cnt.get(k,0)/n_c*100:.1f}%)" for k in "STY"))
    emit(f"全体 {n_a} 位点: " + " | ".join(f"{k}={a_cnt.get(k,0)} ({a_cnt.get(k,0)/n_a*100:.1f}%)" for k in "STY"))
    # 卡方: 候选 vs 非候选 × S/T/Y
    chi = 0.0
    for k in "STY":
        for grp, tot, cnt in (("cand", n_c, c_cnt.get(k, 0)), ("rest", n_a - n_c, a_cnt.get(k, 0) - c_cnt.get(k, 0))):
            exp = a_cnt.get(k, 0) * tot / n_a
            if exp > 0:
                chi += (cnt - exp) ** 2 / exp
    p_chi = math.exp(-chi / 2)   # df=2 精确
    emit(f"卡方 (候选 vs 非候选 × S/T/Y): chi2={chi:.3f}, df=2, p={p_chi:.3e}")
    if c_cnt.get("Y", 0) < 20:
        emit(f"注: Y 仅 {c_cnt.get('Y',0)} 个 — underpowered for Y, 结果仍报告")
    rows_ct = []
    emit("残基 × 解离分类 | 中位|log2occFC| | 峰值组分 top3:")
    for k in "STY":
        sids = [s for s in cand if res(s) == k]
        cl = Counter(cand[s]["site_class"] for s in sids)
        fcs = sorted(float(cand[s]["max_abs_log2occ_fc"]) for s in sids if cand[s]["max_abs_log2occ_fc"])
        med = fcs[len(fcs)//2] if fcs else float("nan")
        frc = Counter(cand[s]["occ_fc_fraction"] for s in sids)
        emit(f"  {k} (n={len(sids)}): " + " ".join(f"{c}={cl.get(c,0)}" for c in ("occ_and_int","occ_not_int"))
             + f" | 中位FC={med:.2f} | " + " ".join(f"{f}:{n}" for f, n in frc.most_common(3)))
        rows_ct.append([k, len(sids), cl.get("occ_and_int",0), cl.get("occ_not_int",0), round(med,3)])
    with open(os.path.join(a7, "sty_class_contingency.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["residue","n","occ_and_int","occ_not_int","median_abs_log2occ_fc"]); w.writerows(rows_ct)
    emit("=== TASK 2 COMPLETE ===\n")

    # ================= TASK 3 =================
    def argmax_at(sid, tp, cols):
        vals = {}
        for fr, r in traj[sid].items():
            v = r.get(f"{cols}_{tp}", "")
            if v != "":
                vals[fr] = float(v)
        if not vals:
            return None, 0
        mx = max(vals.values())
        winners = [f for f, v in vals.items() if v == mx]
        return winners[0], len(winners)

    for metric, col in (("intensity", "log2int"), ("occupancy", "log2occ")):
        M = [[0]*6 for _ in range(6)]
        excl_ctrl = excl_egf = ties = 0
        t_to_used = Counter()
        for sid in cand:
            fr_from, tie1 = argmax_at(sid, "CTRL", col)
            if fr_from is None:
                excl_ctrl += 1; continue
            fr_to, t_used = None, None
            for tp in EGF:
                fr_to, tie2 = argmax_at(sid, tp, col)
                if fr_to is not None:
                    t_used = tp; break
            if fr_to is None:
                excl_egf += 1; continue
            ties += (tie1 > 1) + (tie2 > 1)
            t_to_used[t_used] += 1
            M[FRS.index(fr_from)][FRS.index(fr_to)] += 1
        n_ok = sum(map(sum, M))
        emit(f"=== TASK 3 ({metric} 版): FR_from(行) × FR_to(列), n={n_ok} ===")
        emit("      " + "  ".join(f"{f:>5s}" for f in FRS))
        for i, f in enumerate(FRS):
            emit(f"{f:>5s} " + "  ".join(f"{M[i][j]:5d}" for j in range(6)))
        diag = sum(M[i][i] for i in range(6))
        flows = sorted(((M[i][j], FRS[i], FRS[j]) for i in range(6) for j in range(6) if i != j), reverse=True)[:5]
        emit(f"对角线(峰值组分未变): {diag}/{n_ok} ({diag/n_ok*100:.1f}%)" if n_ok else "无可算位点")
        emit("非对角 top5 流向: " + "; ".join(f"{a}→{b}:{n}" for n, a, b in flows if n > 0))
        emit(f"排除: 无 CTRL 值 {excl_ctrl} | 无 EGF 值 {excl_egf} | 并列取首 {ties} 次")
        emit(f"t_to 实际使用: " + " ".join(f"{t}:{n}" for t, n in t_to_used.most_common()))
        with open(os.path.join(a7, f"flow_matrix_{metric}.tsv"), "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["from\\to"] + FRS)
            for i, fr in enumerate(FRS): w.writerow([fr] + M[i])
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(5.5, 4.8))
            im = ax.imshow(M, cmap="viridis")
            ax.set_xticks(range(6), FRS); ax.set_yticks(range(6), FRS)
            ax.set_xlabel("FR_to (latest EGF tp)"); ax.set_ylabel("FR_from (CTRL)")
            ax.set_title(f"Candidate phospho-form peak-fraction flow ({metric}, metric v1)")
            for i in range(6):
                for j in range(6):
                    ax.text(j, i, M[i][j], ha="center", va="center",
                            color="white" if M[i][j] < max(map(max, M))*0.6 else "black", fontsize=8)
            fig.colorbar(im); fig.tight_layout()
            png = os.path.join(a7, f"flow_heatmap_{metric}.png")
            fig.savefig(png, dpi=150); plt.close(fig)
            emit(f"热力图: {png}")
        except ImportError:
            emit("matplotlib 不可用 — 热力图降级为上方文本矩阵 (数字产出不受影响)")
        emit(f"=== TASK 3 ({metric}) COMPLETE ===\n")
    emit("注记: intensity 版跨组分取 argmax 受组分深度差影响 (FR3 天然浅); occupancy 版部分免疫但缺失更多。两版并列, 不做取舍。")
    with open(os.path.join(a7, "sty_breakdown.txt"), "w") as f:
        f.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
