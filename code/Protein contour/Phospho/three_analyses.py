"""
Three independent analyses:
1. Biophysical features of S/T/Y in phospho-bearing vs non-phospho proteins
2. Auto- vs trans-kinase-substrate correlations
3. Hub (both High+Low) vs Low-only protein length comparison
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import json
import os
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 150

OUTDIR = r'D:\博士\Protein contour\Phospho\three_analyses'
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
# Load data
# ============================================================
success_df = pd.read_csv(r'D:\博士\Protein contour\Phospho\success_df_M1M2_cleaned.csv')
high = pd.read_csv(r'D:\博士\Protein contour\Phospho\high_mobility_sites.csv')
low = pd.read_csv(r'D:\博士\Protein contour\Phospho\low_mobility_sites.csv')

high_uniprots = set(high['uniprot_id'].dropna().unique())
low_uniprots = set(low['uniprot_id'].dropna().unique())
hub_proteins = high_uniprots & low_uniprots          # 364
low_only_proteins = low_uniprots - high_uniprots     # 2467

print(f"Hub proteins (High+Low both):  {len(hub_proteins)}")
print(f"Low-only proteins:             {len(low_only_proteins)}")

# Phospho-bearing proteins (from success_df)
phospho_uniprots = set(success_df['uniprot_id'].dropna().unique())
print(f"Proteins with >=1 phosphosite: {len(phospho_uniprots)}")

# Phospho positions per protein
phos_positions_per_protein = success_df.groupby('uniprot_id')['Position'].apply(set).to_dict()

# All available biophys JSONs
all_jsons = set(os.listdir(BIOPHYS_DIR))
print(f"Total biophys JSON files: {len(all_jsons)}")

def load_biophys(uniprot_id):
    fname = f'{uniprot_id}.json'
    if fname not in all_jsons:
        return None
    try:
        with open(os.path.join(BIOPHYS_DIR, fname)) as f:
            return json.load(f)
    except Exception:
        return None

# ============================================================
# ANALYSIS 1: Biophysical features of S/T/Y residues
# ============================================================
print("\n" + "=" * 60)
print("ANALYSIS 1: S/T/Y residue features in phospho vs non-phospho proteins")
print("=" * 60)

phospho_protein_means = []   # per-protein mean of non-phosphorylated S/T/Y residues
nonphospho_protein_means = [] # per-protein mean of all S/T/Y residues

# Proteins with phosphosites
n_processed_phos = 0
n_no_sty_phos = 0
for up in phospho_uniprots:
    data = load_biophys(up)
    if data is None:
        continue
    phos_pos = phos_positions_per_protein.get(up, set())
    sty_residues = []
    for r in data['residues']:
        if r['aa'] in ('S', 'T', 'Y') and r['seqpos'] not in phos_pos:
            sty_residues.append(r)
    if not sty_residues:
        n_no_sty_phos += 1
        continue
    means = {f: np.mean([r[f] for r in sty_residues if f in r]) for f in FEATURES}
    means['uniprot_id'] = up
    means['n_sty'] = len(sty_residues)
    phospho_protein_means.append(means)
    n_processed_phos += 1

# Non-phospho proteins (all proteins with biophys data but NOT in phospho set)
all_uniprots_with_biophys = {f.replace('.json', '') for f in all_jsons}
nonphospho_uniprots = all_uniprots_with_biophys - phospho_uniprots
print(f"Non-phospho proteins available: {len(nonphospho_uniprots)}")

n_processed_nonphos = 0
n_no_sty_nonphos = 0
for up in nonphospho_uniprots:
    data = load_biophys(up)
    if data is None:
        continue
    sty_residues = [r for r in data['residues'] if r['aa'] in ('S', 'T', 'Y')]
    if not sty_residues:
        n_no_sty_nonphos += 1
        continue
    means = {f: np.mean([r[f] for r in sty_residues if f in r]) for f in FEATURES}
    means['uniprot_id'] = up
    means['n_sty'] = len(sty_residues)
    nonphospho_protein_means.append(means)
    n_processed_nonphos += 1

phos_df = pd.DataFrame(phospho_protein_means)
nonphos_df = pd.DataFrame(nonphospho_protein_means)
print(f"Processed phospho-bearing proteins: {n_processed_phos} (skipped {n_no_sty_phos} with no available S/T/Y)")
print(f"Processed non-phospho proteins: {n_processed_nonphos} (skipped {n_no_sty_nonphos} with no S/T/Y)")

# Stats per feature
print("\n--- Mann-Whitney U tests per feature ---")
stats_rows = []
for f in FEATURES:
    a = phos_df[f].dropna().values
    b = nonphos_df[f].dropna().values
    u, p = mannwhitneyu(a, b, alternative='two-sided')
    stats_rows.append({
        'feature': f,
        'phos_median': np.median(a),
        'nonphos_median': np.median(b),
        'delta_median': np.median(a) - np.median(b),
        'U': u,
        'p_value': p,
    })
    print(f"  {f:14s}: phos median={np.median(a):.3f}, "
          f"non-phos median={np.median(b):.3f}, "
          f"Δ={np.median(a)-np.median(b):+.3f}, p={p:.2e}")

stats_df = pd.DataFrame(stats_rows)
stats_df.to_csv(os.path.join(OUTDIR, 'analysis1_feature_stats.csv'), index=False)

# Violin plot
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
axes = axes.flatten()
for i, f in enumerate(FEATURES):
    ax = axes[i]
    plot_df = pd.concat([
        pd.DataFrame({'value': phos_df[f].dropna(), 'group': 'Phospho-bearing\nproteins'}),
        pd.DataFrame({'value': nonphos_df[f].dropna(), 'group': 'Non-phospho\nproteins'}),
    ])
    sns.violinplot(data=plot_df, x='group', y='value', ax=ax,
                   palette=['#e74c3c', '#3498db'], inner='quartile', cut=0)
    p_val = stats_df[stats_df['feature'] == f]['p_value'].values[0]
    ax.set_title(f"{FEATURE_LABELS[f]}\np = {p_val:.2e}", fontsize=11)
    ax.set_xlabel('')
    ax.set_ylabel('Per-protein mean (S/T/Y)', fontsize=9)

# Hide last subplot
axes[7].set_visible(False)
plt.suptitle(f'Per-Protein S/T/Y Biophysical Features\n'
             f'Phospho-bearing (n={len(phos_df)}, non-phos S/T/Y residues only) vs Non-phospho proteins (n={len(nonphos_df)})',
             fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'analysis1_violin.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: analysis1_violin.png")

# ============================================================
# ANALYSIS 2: Auto- vs trans-kinase correlations
# ============================================================
print("\n" + "=" * 60)
print("ANALYSIS 2: Auto- vs trans-kinase-substrate correlations")
print("=" * 60)

corr_df = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\kinase_substrate_correlation.csv')
corr_df['type'] = np.where(corr_df['Kinase'] == corr_df['Substrate'], 'auto', 'trans')
valid = corr_df.dropna(subset=['Overall_Pearson_r']).copy()

auto_r = valid[valid['type'] == 'auto']['Overall_Pearson_r'].values
trans_r = valid[valid['type'] == 'trans']['Overall_Pearson_r'].values

print(f"Auto pairs (Kinase==Substrate): n={len(auto_r)}")
print(f"  median r = {np.median(auto_r):.3f}, mean = {np.mean(auto_r):.3f}")
print(f"Trans pairs: n={len(trans_r)}")
print(f"  median r = {np.median(trans_r):.3f}, mean = {np.mean(trans_r):.3f}")

u_at, p_at = mannwhitneyu(auto_r, trans_r, alternative='two-sided')
print(f"\nMann-Whitney U: U={u_at:.0f}, p={p_at:.3e}")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(15, 6))
# Histogram overlay
ax = axes[0]
ax.hist(trans_r, bins=50, color='#3498db', alpha=0.6,
        label=f'Trans (n={len(trans_r)})', edgecolor='white', density=True)
ax.hist(auto_r, bins=50, color='#e74c3c', alpha=0.7,
        label=f'Auto (n={len(auto_r)})', edgecolor='white', density=True)
ax.axvline(np.median(trans_r), color='#1f5582', linestyle='--', linewidth=2,
           label=f'Trans median = {np.median(trans_r):.3f}')
ax.axvline(np.median(auto_r), color='#922b21', linestyle='--', linewidth=2,
           label=f'Auto median = {np.median(auto_r):.3f}')
ax.set_xlabel('Overall Pearson r', fontsize=12)
ax.set_ylabel('Density', fontsize=12)
ax.set_title(f'Distribution: Auto- vs Trans-correlation\n(Mann-Whitney p = {p_at:.2e})', fontsize=12)
ax.legend(fontsize=10)

# Violin
ax = axes[1]
plot_df = pd.concat([
    pd.DataFrame({'r': auto_r, 'type': 'Auto\n(Kinase=Substrate)'}),
    pd.DataFrame({'r': trans_r, 'type': 'Trans\n(Kinase≠Substrate)'}),
])
sns.violinplot(data=plot_df, x='type', y='r', ax=ax,
               palette=['#e74c3c', '#3498db'], inner='quartile', cut=0)
ax.axhline(0, color='black', linewidth=0.8, linestyle='-', alpha=0.5)
ax.set_xlabel('')
ax.set_ylabel('Overall Pearson r', fontsize=12)
ax.set_title(f'Auto vs Trans (p = {p_at:.2e})', fontsize=12)

plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'analysis2_auto_vs_trans.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: analysis2_auto_vs_trans.png")

# ============================================================
# ANALYSIS 3: Hub vs Low-only protein length
# ============================================================
print("\n" + "=" * 60)
print("ANALYSIS 3: Hub vs Low-only protein length")
print("=" * 60)

def get_length(uniprot_id):
    data = load_biophys(uniprot_id)
    if data is None:
        return None
    return len(data['residues'])

hub_lengths = []
hub_missing = 0
for up in hub_proteins:
    L = get_length(up)
    if L is None:
        hub_missing += 1
    else:
        hub_lengths.append(L)

low_lengths = []
low_missing = 0
for up in low_only_proteins:
    L = get_length(up)
    if L is None:
        low_missing += 1
    else:
        low_lengths.append(L)

hub_lengths = np.array(hub_lengths)
low_lengths = np.array(low_lengths)

print(f"Hub proteins:      n={len(hub_lengths)} (missing {hub_missing})")
print(f"  median length = {np.median(hub_lengths):.0f}, mean = {np.mean(hub_lengths):.0f}")
print(f"Low-only proteins: n={len(low_lengths)} (missing {low_missing})")
print(f"  median length = {np.median(low_lengths):.0f}, mean = {np.mean(low_lengths):.0f}")

u_hl, p_hl = mannwhitneyu(hub_lengths, low_lengths, alternative='two-sided')
print(f"\nMann-Whitney U: U={u_hl:.0f}, p={p_hl:.3e}")

# Violin plot — log-scale length
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
plot_df = pd.concat([
    pd.DataFrame({'length': hub_lengths, 'group': f'Hub\n(High+Low, n={len(hub_lengths)})'}),
    pd.DataFrame({'length': low_lengths, 'group': f'Low-only\n(n={len(low_lengths)})'}),
])

ax = axes[0]
sns.violinplot(data=plot_df, x='group', y='length', ax=ax,
               palette=['#9b59b6', '#3498db'], inner='quartile', cut=0)
ax.set_yscale('log')
ax.set_xlabel('')
ax.set_ylabel('Protein length (residues, log scale)', fontsize=12)
ax.set_title(f'Protein Length Distribution\n(Mann-Whitney p = {p_hl:.2e})', fontsize=12)

ax = axes[1]
sns.violinplot(data=plot_df, x='group', y='length', ax=ax,
               palette=['#9b59b6', '#3498db'], inner='quartile', cut=0)
ax.set_xlabel('')
ax.set_ylabel('Protein length (residues, linear)', fontsize=12)
ax.set_title(f'Linear scale\nHub median={np.median(hub_lengths):.0f}, Low-only median={np.median(low_lengths):.0f}',
             fontsize=12)

plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'analysis3_protein_length.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: analysis3_protein_length.png")

# ============================================================
# Final summary
# ============================================================
print("\n" + "=" * 60)
print("FINAL SUMMARY")
print("=" * 60)

print(f"\n[Analysis 1] Biophysical features (S/T/Y, per-protein means):")
print(f"  Phospho-bearing proteins: n={len(phos_df)}")
print(f"  Non-phospho proteins:     n={len(nonphos_df)}")
print(f"  Significant features (p<0.05): {(stats_df['p_value']<0.05).sum()}/7")
top_feat = stats_df.iloc[stats_df['p_value'].argmin()]
print(f"  Strongest difference: {top_feat['feature']} (Δmedian={top_feat['delta_median']:+.3f}, p={top_feat['p_value']:.2e})")

print(f"\n[Analysis 2] Kinase-substrate correlations:")
print(f"  Auto  (n={len(auto_r)}): median r = {np.median(auto_r):.3f}")
print(f"  Trans (n={len(trans_r)}): median r = {np.median(trans_r):.3f}")
print(f"  Mann-Whitney U p = {p_at:.3e}")

print(f"\n[Analysis 3] Protein length:")
print(f"  Hub (n={len(hub_lengths)}): median={np.median(hub_lengths):.0f} aa")
print(f"  Low-only (n={len(low_lengths)}): median={np.median(low_lengths):.0f} aa")
print(f"  Mann-Whitney U p = {p_hl:.3e}")

print(f"\nFiles saved to: {OUTDIR}")
for f in sorted(os.listdir(OUTDIR)):
    print(f"  {f} ({os.path.getsize(os.path.join(OUTDIR, f))/1024:.0f} KB)")
