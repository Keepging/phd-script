"""
Analysis 1: Length-matched Hub vs Low-only — 7 biophys features
Analysis 2: Per-compartment temporal correlation vs random baseline
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, pearsonr
import json
import os
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 150

OUTDIR = r'D:\博士\Protein contour\Phospho\length_matched_temporal'
os.makedirs(OUTDIR, exist_ok=True)

BIOPHYS_DIR = r'D:\博士\Protein contour\Phospho\biophys_json'
FEATURES = ['backbone', 'sidechain', 'disoMine', 'helix', 'sheet', 'coil', 'earlyFolding']
FEATURE_LABELS = {
    'backbone': 'Backbone dynamics',
    'sidechain': 'Sidechain dynamics',
    'disoMine': 'Disorder propensity',
    'helix': 'Helix propensity',
    'sheet': 'Sheet propensity',
    'coil': 'Coil propensity',
    'earlyFolding': 'Early-folding propensity',
}

# ============================================================
# Reconstruct hub / low-only lists
# ============================================================
high = pd.read_csv(r'D:\博士\Protein contour\Phospho\high_mobility_sites.csv')
low = pd.read_csv(r'D:\博士\Protein contour\Phospho\low_mobility_sites.csv')
high_uniprots = set(high['uniprot_id'].dropna().unique())
low_uniprots = set(low['uniprot_id'].dropna().unique())
hub_proteins = high_uniprots & low_uniprots
low_only_proteins = low_uniprots - high_uniprots
print(f"Hub: {len(hub_proteins)}, Low-only: {len(low_only_proteins)}")

# ============================================================
# Load biophys + compute per-protein S/T/Y feature means + lengths
# ============================================================
def load_biophys(uniprot_id):
    fname = f'{uniprot_id}.json'
    fp = os.path.join(BIOPHYS_DIR, fname)
    if not os.path.exists(fp):
        return None
    try:
        with open(fp) as f:
            return json.load(f)
    except Exception:
        return None

def compute_protein_record(up):
    data = load_biophys(up)
    if data is None:
        return None
    sty = [r for r in data['residues'] if r['aa'] in ('S', 'T', 'Y')]
    if not sty:
        return None
    rec = {'uniprot_id': up, 'length': len(data['residues']), 'n_sty': len(sty)}
    for f in FEATURES:
        vals = [r[f] for r in sty if f in r]
        rec[f] = np.mean(vals) if vals else np.nan
    return rec

print("Computing per-protein S/T/Y feature means + lengths...")
hub_records = [compute_protein_record(up) for up in hub_proteins]
hub_records = [r for r in hub_records if r is not None]
low_records = [compute_protein_record(up) for up in low_only_proteins]
low_records = [r for r in low_records if r is not None]

hub_df = pd.DataFrame(hub_records)
low_df = pd.DataFrame(low_records)
print(f"Hub  records: {len(hub_df)} (length median={hub_df['length'].median():.0f})")
print(f"Low-only records: {len(low_df)} (length median={low_df['length'].median():.0f})")

# ============================================================
# ANALYSIS 1
# ============================================================
print("\n" + "=" * 60)
print("ANALYSIS 1: Length-matched Hub vs Low-only")
print("=" * 60)

# --- Original (unmatched) MW U test ---
print("\n--- Original (unmatched) ---")
original_results = {}
for f in FEATURES:
    a = hub_df[f].dropna().values
    b = low_df[f].dropna().values
    u, p = mannwhitneyu(a, b, alternative='two-sided')
    delta = np.median(a) - np.median(b)
    original_results[f] = {
        'orig_p': p,
        'orig_delta_median': delta,
        'hub_median': np.median(a),
        'low_median': np.median(b),
    }
    print(f"  {f:14s}: hub={np.median(a):.3f}, low={np.median(b):.3f}, "
          f"Δ={delta:+.4f}, p={p:.2e}")

# --- Length-matched (1000 iterations, bin = 100 aa) ---
print("\n--- Length-matched (bin=100 aa, 1000 iterations) ---")
BIN_SIZE = 100
N_ITER = 1000

hub_df['length_bin'] = (hub_df['length'] // BIN_SIZE).astype(int)
low_df['length_bin'] = (low_df['length'] // BIN_SIZE).astype(int)

# Hub bin counts
hub_bin_counts = hub_df['length_bin'].value_counts().to_dict()

# Pre-index low-only by bin for fast sampling
low_by_bin = {b: low_df[low_df['length_bin'] == b].reset_index(drop=True)
              for b in hub_bin_counts.keys()}

# Check coverage
shortfall_bins = []
for b, n in hub_bin_counts.items():
    avail = len(low_by_bin.get(b, []))
    if avail < n:
        shortfall_bins.append((b, n, avail))
if shortfall_bins:
    print(f"  WARNING: {len(shortfall_bins)} bins have fewer low-only than hub:")
    for b, n, avail in shortfall_bins[:5]:
        print(f"    bin {b*BIN_SIZE}-{(b+1)*BIN_SIZE} aa: hub={n}, low={avail}")

rng = np.random.default_rng(seed=42)

# Storage: per feature, list of p-values and deltas
matched_pvals = {f: [] for f in FEATURES}
matched_deltas = {f: [] for f in FEATURES}

# Pre-extract feature values per protein (fast access)
hub_feat = {f: hub_df[f].values for f in FEATURES}
low_feat = {f: low_df[f].values for f in FEATURES}

# Track length-matched low subset across iterations to verify
matched_length_medians = []

for it in range(N_ITER):
    # Sample matched low-only indices
    matched_idx = []
    for b, n in hub_bin_counts.items():
        pool = low_by_bin.get(b)
        if pool is None or len(pool) == 0:
            continue
        if len(pool) >= n:
            sampled = rng.choice(len(pool), size=n, replace=False)
        else:
            # If pool smaller, sample with replacement
            sampled = rng.choice(len(pool), size=n, replace=True)
        # Get original df index
        orig_idx = pool.index[sampled].tolist() if hasattr(pool.index, 'tolist') else list(sampled)
        # We stored low_by_bin with reset_index, so 'sampled' are positions inside the original row order
        # Need to map back via stored uniprot_ids
        matched_idx.extend(pool.iloc[sampled]['uniprot_id'].tolist())

    # Get matched low subset rows
    matched_low = low_df[low_df['uniprot_id'].isin(matched_idx)]
    # In case of duplicates (sampled with replacement), we lose them via .isin; rebuild via merge
    # For correctness use a proper merge by uniprot_id with multiplicity
    # but for MW test, multiplicity within the matched set matters; use list-based assembly
    matched_uniprots_with_mult = matched_idx  # list with possible duplicates (rare)
    # Build feature vectors with multiplicity
    up_to_row = low_df.set_index('uniprot_id')

    if it == 0:
        # First iteration: track length distribution to verify matching
        try:
            matched_lengths = up_to_row.loc[matched_uniprots_with_mult, 'length'].values
            matched_length_medians.append(np.median(matched_lengths))
            print(f"  Iter 1: matched low length median = {np.median(matched_lengths):.0f} "
                  f"(hub median = {hub_df['length'].median():.0f})")
        except Exception as e:
            print(f"  Length check failed: {e}")

    for f in FEATURES:
        try:
            b_vals = up_to_row.loc[matched_uniprots_with_mult, f].dropna().values
        except KeyError:
            continue
        a_vals = hub_feat[f]
        a_vals = a_vals[~np.isnan(a_vals)]
        if len(b_vals) < 5 or len(a_vals) < 5:
            continue
        u, p = mannwhitneyu(a_vals, b_vals, alternative='two-sided')
        matched_pvals[f].append(p)
        matched_deltas[f].append(np.median(a_vals) - np.median(b_vals))

    if (it + 1) % 200 == 0:
        print(f"  Completed {it+1}/{N_ITER} iterations")

# Aggregate matched results
print("\n--- Length-matched results ---")
matched_results = {}
for f in FEATURES:
    pvals = np.array(matched_pvals[f])
    deltas = np.array(matched_deltas[f])
    matched_results[f] = {
        'matched_median_p': np.median(pvals),
        'matched_pct_sig': (pvals < 0.05).mean() * 100,
        'matched_median_delta': np.median(deltas),
    }
    print(f"  {f:14s}: median p={np.median(pvals):.3e}, "
          f"%sig(<0.05)={(pvals<0.05).mean()*100:.0f}%, "
          f"Δ_median={np.median(deltas):+.4f}")

# Combine into summary table
summary_rows = []
for f in FEATURES:
    o = original_results[f]
    m = matched_results[f]
    summary_rows.append({
        'feature': f,
        'hub_median': o['hub_median'],
        'low_median': o['low_median'],
        'orig_delta_median': o['orig_delta_median'],
        'orig_p': o['orig_p'],
        'matched_median_delta': m['matched_median_delta'],
        'matched_median_p': m['matched_median_p'],
        'matched_pct_iter_sig': m['matched_pct_sig'],
        'effect_attenuation': abs(m['matched_median_delta']) / abs(o['orig_delta_median']) if o['orig_delta_median'] != 0 else np.nan,
    })
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUTDIR, 'analysis1_length_matched_summary.csv'), index=False)
print(f"\nSaved: analysis1_length_matched_summary.csv")
print(summary_df.to_string(index=False))

# Conclusion
print("\n--- Conclusion: which features survive length matching? ---")
for _, row in summary_df.iterrows():
    f = row['feature']
    survives = row['matched_median_p'] < 0.05
    pct_iter_sig = row['matched_pct_iter_sig']
    attenuation = (1 - row['effect_attenuation']) * 100
    flag = '✓ SURVIVES' if survives and pct_iter_sig > 50 else '✗ LOST'
    print(f"  {f:14s}: {flag}  (median p={row['matched_median_p']:.2e}, "
          f"effect attenuated by {attenuation:.0f}%)")

# Violin plot: original vs matched (use last iteration's matched set for visualization)
# Use median-p iteration for fairness — pick iteration with median behavior
final_matched_uniprots = matched_uniprots_with_mult  # last iteration
final_matched_low = up_to_row.loc[final_matched_uniprots].reset_index()

fig, axes = plt.subplots(2, 7, figsize=(28, 10))
for i, f in enumerate(FEATURES):
    # Top row: original
    ax_top = axes[0][i]
    plot_df_orig = pd.concat([
        pd.DataFrame({'value': hub_df[f].dropna(), 'group': 'Hub'}),
        pd.DataFrame({'value': low_df[f].dropna(), 'group': 'Low-only'}),
    ])
    sns.violinplot(data=plot_df_orig, x='group', y='value', ax=ax_top,
                   palette=['#9b59b6', '#3498db'], inner='quartile', cut=0)
    p_orig = original_results[f]['orig_p']
    ax_top.set_title(f"{FEATURE_LABELS[f]}\nORIGINAL: p={p_orig:.1e}", fontsize=10)
    ax_top.set_xlabel('')
    ax_top.set_ylabel('Per-protein S/T/Y mean' if i == 0 else '', fontsize=9)

    # Bottom row: length-matched
    ax_bot = axes[1][i]
    plot_df_matched = pd.concat([
        pd.DataFrame({'value': hub_df[f].dropna(), 'group': 'Hub'}),
        pd.DataFrame({'value': final_matched_low[f].dropna(), 'group': 'Low-only\n(matched)'}),
    ])
    sns.violinplot(data=plot_df_matched, x='group', y='value', ax=ax_bot,
                   palette=['#9b59b6', '#85c1e9'], inner='quartile', cut=0)
    p_med = matched_results[f]['matched_median_p']
    pct_sig = matched_results[f]['matched_pct_sig']
    ax_bot.set_title(f"MATCHED: median p={p_med:.1e}\n%iter sig={pct_sig:.0f}%", fontsize=10)
    ax_bot.set_xlabel('')
    ax_bot.set_ylabel('Per-protein S/T/Y mean' if i == 0 else '', fontsize=9)

plt.suptitle(f'Hub (n={len(hub_df)}) vs Low-only — Original vs Length-matched (1000 iter, bin=100 aa)\n'
             f'Top: original | Bottom: length-matched (last iteration shown)',
             fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'analysis1_violin_orig_vs_matched.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: analysis1_violin_orig_vs_matched.png")

# ============================================================
# ANALYSIS 2: Per-compartment temporal correlation
# ============================================================
print("\n" + "=" * 60)
print("ANALYSIS 2: Per-compartment temporal correlation")
print("=" * 60)

corr_df = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\kinase_substrate_correlation.csv')
abund_matrix = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\protein_abundance_matrix.csv',
                            index_col=0)

# Filter auto-phos
nonself = corr_df[corr_df['Kinase'] != corr_df['Substrate']].copy()
print(f"Total pairs: {len(corr_df)}, non-self: {len(nonself)}")

# Distributions
ovr_r = nonself['Overall_Pearson_r'].dropna().values
cyt_r = nonself['Cyt_r'].dropna().values
mem_r = nonself['Mem_r'].dropna().values
nuc_r = nonself['Nuc_r'].dropna().values

print(f"\n--- Real kinase-substrate (non-self) ---")
print(f"  Overall  : n={len(ovr_r):4d}, median r={np.median(ovr_r):.3f}, mean={np.mean(ovr_r):.3f}")
print(f"  Cyt only : n={len(cyt_r):4d}, median r={np.median(cyt_r):.3f}, mean={np.mean(cyt_r):.3f}")
print(f"  Mem only : n={len(mem_r):4d}, median r={np.median(mem_r):.3f}, mean={np.mean(mem_r):.3f}")
print(f"  Nuc only : n={len(nuc_r):4d}, median r={np.median(nuc_r):.3f}, mean={np.mean(nuc_r):.3f}")

# --- Random baseline: per-compartment + overall ---
print("\nComputing random baseline (per-compartment)...")
timepoints = ['0min', '2min', '8min', '20min', '90min']
compartments = ['Cyt', 'Mem', 'Nuc']
cond_order = [f'{c}_{t}' for c in compartments for t in timepoints]
cyt_idx = [cond_order.index(f'Cyt_{t}') for t in timepoints]
mem_idx = [cond_order.index(f'Mem_{t}') for t in timepoints]
nuc_idx = [cond_order.index(f'Nuc_{t}') for t in timepoints]

abund_vals = abund_matrix[cond_order].values.astype(float)
n_genes = len(abund_matrix)
N_RAND = len(nonself)  # match size of real

def safe_pearson(x, y, min_pairs=3):
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

rng2 = np.random.default_rng(seed=2024)
rand_overall = []
rand_cyt = []
rand_mem = []
rand_nuc = []

attempts = 0
while len(rand_overall) < N_RAND and attempts < N_RAND * 8:
    a = rng2.integers(0, n_genes)
    b = rng2.integers(0, n_genes)
    if a == b:
        attempts += 1
        continue
    va = abund_vals[a]
    vb = abund_vals[b]
    r_all = safe_pearson(va, vb, min_pairs=5)
    if not np.isnan(r_all):
        rand_overall.append(r_all)
    rand_cyt.append(safe_pearson(va[cyt_idx], vb[cyt_idx], min_pairs=3))
    rand_mem.append(safe_pearson(va[mem_idx], vb[mem_idx], min_pairs=3))
    rand_nuc.append(safe_pearson(va[nuc_idx], vb[nuc_idx], min_pairs=3))
    attempts += 1

rand_overall = np.array(rand_overall)
rand_cyt = np.array([x for x in rand_cyt if not np.isnan(x)])
rand_mem = np.array([x for x in rand_mem if not np.isnan(x)])
rand_nuc = np.array([x for x in rand_nuc if not np.isnan(x)])

print(f"\n--- Random baseline ---")
print(f"  Overall  : n={len(rand_overall):4d}, median r={np.median(rand_overall):.3f}")
print(f"  Cyt only : n={len(rand_cyt):4d}, median r={np.median(rand_cyt):.3f}")
print(f"  Mem only : n={len(rand_mem):4d}, median r={np.median(rand_mem):.3f}")
print(f"  Nuc only : n={len(rand_nuc):4d}, median r={np.median(rand_nuc):.3f}")

# Mann-Whitney for each
results2 = []
for label, real, rand in [('Overall', ovr_r, rand_overall),
                           ('Cyt', cyt_r, rand_cyt),
                           ('Mem', mem_r, rand_mem),
                           ('Nuc', nuc_r, rand_nuc)]:
    u, p = mannwhitneyu(real, rand, alternative='two-sided')
    results2.append({
        'level': label,
        'real_n': len(real),
        'real_median': np.median(real),
        'rand_n': len(rand),
        'rand_median': np.median(rand),
        'delta_median': np.median(real) - np.median(rand),
        'mw_u': u,
        'mw_p': p,
    })

results2_df = pd.DataFrame(results2)
results2_df.to_csv(os.path.join(OUTDIR, 'analysis2_per_compartment_summary.csv'), index=False)
print(f"\nSaved: analysis2_per_compartment_summary.csv")
print(results2_df.to_string(index=False))

# Conclusion
print("\n--- Conclusion: temporal signal? ---")
for _, row in results2_df.iterrows():
    if row['level'] == 'Overall':
        continue
    delta = row['delta_median']
    if abs(delta) < 0.05:
        flag = '≈ NO temporal signal'
    elif delta > 0.1:
        flag = '✓ TEMPORAL signal present'
    else:
        flag = '~ weak temporal signal'
    print(f"  {row['level']:7s}: real median={row['real_median']:.3f}, "
          f"random={row['rand_median']:.3f}, Δ={delta:+.3f} → {flag}")

# 4-panel violin plot
fig, axes = plt.subplots(1, 4, figsize=(20, 6))
panels = [
    ('Overall (15 conditions)', ovr_r, rand_overall),
    ('Cyt only (5 timepoints)', cyt_r, rand_cyt),
    ('Mem only (5 timepoints)', mem_r, rand_mem),
    ('Nuc only (5 timepoints)', nuc_r, rand_nuc),
]
for i, (title, real, rand) in enumerate(panels):
    ax = axes[i]
    plot_df = pd.concat([
        pd.DataFrame({'r': real, 'group': f'Kinase-substrate\n(n={len(real)})'}),
        pd.DataFrame({'r': rand, 'group': f'Random\n(n={len(rand)})'}),
    ])
    sns.violinplot(data=plot_df, x='group', y='r', ax=ax,
                   palette=['#e74c3c', '#95a5a6'], inner='quartile', cut=0)
    ax.axhline(0, color='black', linewidth=0.8, linestyle='-', alpha=0.5)
    p_val = results2[i]['mw_p']
    delta = results2[i]['delta_median']
    ax.set_title(f"{title}\nΔmedian={delta:+.3f}, p={p_val:.1e}", fontsize=11)
    ax.set_xlabel('')
    if i == 0:
        ax.set_ylabel('Pearson r', fontsize=12)
    else:
        ax.set_ylabel('')
    ax.set_ylim(-1.05, 1.05)

plt.suptitle('Spatial vs Temporal Decomposition of Kinase-Substrate Correlation\n'
             '(Non-self pairs; per-compartment correlations use only 5 timepoints within compartment)',
             fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'analysis2_4panel_violin.png'), dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved: analysis2_4panel_violin.png")

# ============================================================
# Final summary
# ============================================================
print("\n" + "=" * 60)
print("FINAL SUMMARY")
print("=" * 60)

print(f"\n[Analysis 1: Length-matched Hub vs Low-only]")
n_survive = sum(1 for _, row in summary_df.iterrows()
                if row['matched_median_p'] < 0.05 and row['matched_pct_iter_sig'] > 50)
print(f"  {n_survive}/7 features survive length matching (p<0.05 in >50% of iterations)")
for _, row in summary_df.iterrows():
    f = row['feature']
    survives = row['matched_median_p'] < 0.05 and row['matched_pct_iter_sig'] > 50
    print(f"    {f:14s}: orig p={row['orig_p']:.1e} → matched p={row['matched_median_p']:.2e} "
          f"({'✓' if survives else '✗'})")

print(f"\n[Analysis 2: Per-compartment temporal correlation]")
for _, row in results2_df.iterrows():
    print(f"  {row['level']:7s}: real={row['real_median']:.3f}, "
          f"random={row['rand_median']:.3f}, Δ={row['delta_median']:+.3f}, p={row['mw_p']:.1e}")

print(f"\nFiles saved to {OUTDIR}:")
for f in sorted(os.listdir(OUTDIR)):
    fsize = os.path.getsize(os.path.join(OUTDIR, f))
    print(f"  {f} ({fsize/1024:.0f} KB)")
