"""
Task A: Null distribution (random gene pairs) vs real kinase-substrate correlations
Task B: Top pairs after filtering self-correlations
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, mannwhitneyu
import os
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 150

outdir = r'D:\博士\Protein contour\Phospho\kinase_correlation'

timepoints = ['0min', '2min', '8min', '20min', '90min']
compartments = ['Cyt', 'Mem', 'Nuc']
cond_order = [f'{c}_{t}' for c in compartments for t in timepoints]

# ============================================================
# Load data
# ============================================================
abund_matrix = pd.read_csv(os.path.join(outdir, 'protein_abundance_matrix.csv'), index_col=0)
corr_df = pd.read_csv(os.path.join(outdir, 'kinase_substrate_correlation.csv'))

print(f"Abundance matrix: {abund_matrix.shape}")
print(f"Correlation table: {corr_df.shape}")

# Filter to non-self pairs
nonself = corr_df[corr_df['Kinase'] != corr_df['Substrate']].copy()
nonself_valid = nonself.dropna(subset=['Overall_Pearson_r']).copy()
N = len(nonself_valid)
print(f"Non-self pairs with valid correlation: {N}")

# ============================================================
# TASK A: Null distribution
# ============================================================
print("\n" + "=" * 60)
print("TASK A: Null distribution from random gene pairs")
print("=" * 60)

# Precompute values matrix (N_genes x 15), numpy for speed
abund_vals = abund_matrix[cond_order].values.astype(float)
all_genes = abund_matrix.index.tolist()
n_genes = len(all_genes)

def safe_pearson(x, y, min_pairs=5):
    mask = (~np.isnan(x)) & (~np.isnan(y))
    if mask.sum() < min_pairs:
        return np.nan
    xm = x[mask]
    ym = y[mask]
    if np.std(xm) == 0 or np.std(ym) == 0:
        return np.nan
    try:
        r, _ = pearsonr(xm, ym)
        return r
    except Exception:
        return np.nan

# 3 rounds of random sampling
round_stats = []
all_random_r = []

rng = np.random.default_rng(seed=42)

for round_idx in range(3):
    rs = []
    # Sample N pairs of distinct gene indices
    # Use rejection-free sampling: sample 2*N indices and take pairs
    attempts = 0
    while len(rs) < N and attempts < N * 10:
        i1 = rng.integers(0, n_genes, size=N)
        i2 = rng.integers(0, n_genes, size=N)
        for a, b in zip(i1, i2):
            if a == b:
                continue
            r = safe_pearson(abund_vals[a], abund_vals[b])
            if not np.isnan(r):
                rs.append(r)
            if len(rs) >= N:
                break
        attempts += 1

    rs = np.array(rs[:N])
    all_random_r.append(rs)
    round_stats.append({
        'round': round_idx + 1,
        'n': len(rs),
        'mean': np.mean(rs),
        'median': np.median(rs),
        'std': np.std(rs),
    })
    print(f"Round {round_idx+1}: n={len(rs)}, mean={np.mean(rs):.3f}, "
          f"median={np.median(rs):.3f}, std={np.std(rs):.3f}")

# Aggregate across rounds
random_all = np.concatenate(all_random_r)
random_mean = np.mean([s['mean'] for s in round_stats])
random_median = np.mean([s['median'] for s in round_stats])
random_std = np.mean([s['std'] for s in round_stats])

real_r = nonself_valid['Overall_Pearson_r'].values
real_mean = np.mean(real_r)
real_median = np.median(real_r)
real_std = np.std(real_r)

print("\n--- Summary (averaged over 3 rounds for random) ---")
print(f"Random pairs:           mean={random_mean:.3f}, median={random_median:.3f}, std={random_std:.3f}")
print(f"Kinase-substrate pairs: mean={real_mean:.3f}, median={real_median:.3f}, std={real_std:.3f}")
print(f"Δmean   = {real_mean - random_mean:+.3f}")
print(f"Δmedian = {real_median - random_median:+.3f}")

# Mann-Whitney U test (use round 1 random vs real)
u_stat, u_p = mannwhitneyu(real_r, all_random_r[0], alternative='two-sided')
print(f"\nMann-Whitney U test (real vs random round 1):")
print(f"  U = {u_stat:.0f}, p = {u_p:.3e}")

# Conclusion
print("\n--- Conclusion ---")
if random_median > 0.4:
    print(f"⚠️ Random median r = {random_median:.3f} > 0.4 — SEVERE compartment confound.")
    print("   Positive correlations in kinase-substrate pairs are NOT trustworthy as evidence of regulation.")
elif random_median < 0.2:
    print(f"✓ Random median r = {random_median:.3f} < 0.2 — Kinase-substrate positive correlations")
    print("  represent a REAL signal above the background.")
else:
    print(f"⚠️ Random median r = {random_median:.3f} is in grey zone (0.2-0.4).")
    print("  Some confound exists; interpret kinase-substrate correlations with caution.")
    print(f"  Δmedian = {real_median - random_median:+.3f} over null — {'still meaningful' if (real_median - random_median) > 0.1 else 'marginal'}")

# ============================================================
# TASK B: Non-self top pairs analysis
# ============================================================
print("\n" + "=" * 60)
print("TASK B: Non-self top pairs (filtered)")
print("=" * 60)

sig = nonself_valid[nonself_valid['Overall_Pearson_p'] < 0.05].copy()
print(f"Non-self + significant (p<0.05): {len(sig)}")
print(f"  Positive: {(sig['Overall_Pearson_r'] > 0).sum()}")
print(f"  Negative: {(sig['Overall_Pearson_r'] < 0).sum()}")

top20_pos = sig.nlargest(20, 'Overall_Pearson_r').reset_index(drop=True)
top20_neg = sig.nsmallest(20, 'Overall_Pearson_r').reset_index(drop=True)

print("\nTop 20 POSITIVELY correlated (non-self, p<0.05):")
print(top20_pos[['Kinase', 'Substrate', 'Sites', 'Overall_Pearson_r', 'Overall_Pearson_p']].to_string(index=False))

print("\nTop 20 NEGATIVELY correlated (non-self, p<0.05):")
print(top20_neg[['Kinase', 'Substrate', 'Sites', 'Overall_Pearson_r', 'Overall_Pearson_p']].to_string(index=False))

# --- Helper: z-score a vector ignoring NaN ---
def zscore(v):
    v = v.astype(float)
    mask = ~np.isnan(v)
    if mask.sum() < 2:
        return np.zeros_like(v)
    m = np.nanmean(v)
    s = np.nanstd(v)
    if s == 0 or np.isnan(s):
        return np.zeros_like(v)
    out = (v - m) / s
    out[~mask] = 0  # for plotting fill NaN with 0
    return out

# --- Figure 1: top 20 positive heatmap ---
def make_pair_heatmap(top_df, title, filename):
    rows = []
    labels = []
    for _, row in top_df.iterrows():
        k = row['Kinase']
        s = row['Substrate']
        if k not in abund_matrix.index or s not in abund_matrix.index:
            continue
        kv = abund_matrix.loc[k, cond_order].values.astype(float)
        sv = abund_matrix.loc[s, cond_order].values.astype(float)
        rows.append(zscore(kv))
        labels.append(f"[K] {k} → {s} (r={row['Overall_Pearson_r']:.2f}, p={row['Overall_Pearson_p']:.1e})")
        rows.append(zscore(sv))
        labels.append(f"[S] {k} → {s}")

    mat = pd.DataFrame(rows, index=labels, columns=cond_order)
    fig, ax = plt.subplots(figsize=(14, 16))
    sns.heatmap(mat, cmap='RdBu_r', center=0, ax=ax,
                cbar_kws={'label': 'Z-score (per row)'},
                linewidths=0.3, linecolor='#eeeeee',
                vmin=-2, vmax=2)
    ax.set_title(title, fontsize=13)
    ax.set_xlabel('Condition', fontsize=11)
    # Compartment separators
    for i in [5, 10]:
        ax.axvline(x=i, color='black', linewidth=2)
    # Pair separators
    for i in range(2, len(labels), 2):
        ax.axhline(y=i, color='black', linewidth=0.8)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8, rotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, filename), dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename}")

make_pair_heatmap(top20_pos,
                  'Top 20 Positively Correlated Kinase→Substrate Pairs (non-self, p<0.05)',
                  'taskB_top20_positive.png')
make_pair_heatmap(top20_neg,
                  'Top 20 Negatively Correlated Kinase→Substrate Pairs (non-self, p<0.05)',
                  'taskB_top20_negative.png')

# --- Figure 3: scatter examples ---
top3_pos = sig.nlargest(3, 'Overall_Pearson_r').reset_index(drop=True)
top3_neg = sig.nsmallest(3, 'Overall_Pearson_r').reset_index(drop=True)
examples = pd.concat([top3_pos, top3_neg]).reset_index(drop=True)

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
comp_colors = {'Cyt': '#1f77b4', 'Mem': '#ff7f0e', 'Nuc': '#2ca02c'}

for idx, (_, row) in enumerate(examples.iterrows()):
    ax = axes[idx // 3][idx % 3]
    k = row['Kinase']
    s = row['Substrate']
    kv = abund_matrix.loc[k, cond_order].values.astype(float)
    sv = abund_matrix.loc[s, cond_order].values.astype(float)

    for i, cond in enumerate(cond_order):
        comp, tp = cond.split('_')
        if np.isnan(kv[i]) or np.isnan(sv[i]):
            continue
        ax.scatter(kv[i], sv[i], color=comp_colors[comp], s=140, alpha=0.85,
                   edgecolors='black', linewidth=0.6)
        ax.annotate(tp.replace('min',''), (kv[i], sv[i]), fontsize=9,
                    xytext=(6, 6), textcoords='offset points')

    for comp, c in comp_colors.items():
        ax.scatter([], [], color=c, s=80, label=comp)
    ax.legend(loc='best', fontsize=9)

    sign = '(+)' if row['Overall_Pearson_r'] > 0 else '(-)'
    ax.set_title(f"{sign} {k} → {s}\nr={row['Overall_Pearson_r']:.3f}, p={row['Overall_Pearson_p']:.1e}",
                 fontsize=11)
    ax.set_xlabel(f'{k} abundance', fontsize=10)
    ax.set_ylabel(f'{s} abundance', fontsize=10)
    ax.grid(True, alpha=0.3)

plt.suptitle('Non-self Kinase-Substrate Correlation Examples (Top 3 +/-)',
             fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'taskB_scatter_nonself.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: taskB_scatter_nonself.png")

# Save non-self correlation table too
nonself_sig = sig.sort_values('Overall_Pearson_p')
nonself_sig.to_csv(os.path.join(outdir, 'kinase_substrate_correlation_nonself_sig.csv'), index=False)
print(f"\nSaved: kinase_substrate_correlation_nonself_sig.csv ({len(nonself_sig)} rows)")

# ============================================================
# Final summary
# ============================================================
print("\n" + "=" * 60)
print("FINAL SUMMARY")
print("=" * 60)
print(f"N (non-self kinase-substrate pairs): {N}")
print(f"Random null distribution (3 rounds averaged):")
print(f"  mean={random_mean:.3f}, median={random_median:.3f}, std={random_std:.3f}")
print(f"Real kinase-substrate distribution:")
print(f"  mean={real_mean:.3f}, median={real_median:.3f}, std={real_std:.3f}")
print(f"Δmean={real_mean-random_mean:+.3f}, Δmedian={real_median-random_median:+.3f}")
print(f"Mann-Whitney U p = {u_p:.3e}")
print(f"\nFiles saved:")
for f in ['taskB_top20_positive.png', 'taskB_top20_negative.png',
          'taskB_scatter_nonself.png', 'kinase_substrate_correlation_nonself_sig.csv']:
    fp = os.path.join(outdir, f)
    if os.path.exists(fp):
        print(f"  {f} ({os.path.getsize(fp)/1024:.0f} KB)")
