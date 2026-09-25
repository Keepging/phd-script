"""
Data-driven top-20 kinase heatmap — replaces the manual 18-kinase panel.

Ranking: top 20 by SD of log10(abundance) across 15 spatiotemporal conditions
(equivalently: log fold-change spread). Picks the kinases with the most
informative spatiotemporal patterns.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 150

OUTDIR = r'D:\博士\Protein contour\Phospho\kinase_correlation'

timepoints = ['0min', '2min', '8min', '20min', '90min']
compartments = ['Cyt', 'Mem', 'Nuc']
cond_order = [f'{c}_{t}' for c in compartments for t in timepoints]

# ============================================================
# Load data
# ============================================================
abund = pd.read_csv(os.path.join(OUTDIR, 'protein_abundance_matrix.csv'), index_col=0)
print(f"Abundance matrix: {abund.shape}")

# All KinAce kinase names
kinase_df = pd.read_csv(r'C:\Users\cc\Zotero\storage\LX2UQSAS\phospho_data\kinace_ksi_source_full.csv')
kinase_df = kinase_df.dropna(subset=['Kinase Name'])
all_kinases = kinase_df['Kinase Name'].str.strip().unique()
kinases_in_matrix = [k for k in all_kinases if k in abund.index]
print(f"Kinases with abundance data: {len(kinases_in_matrix)}")

# Restrict to kinase rows
kinase_abund = abund.loc[kinases_in_matrix].copy()

# Require >= 8 non-NaN values across 15 conditions (at least half coverage)
n_nonan = kinase_abund.notna().sum(axis=1)
kinase_abund = kinase_abund[n_nonan >= 8]
print(f"Kinases with >=8/15 non-NaN values: {len(kinase_abund)}")

# ============================================================
# Compute multiple ranking metrics
# ============================================================
log_abund = np.log10(kinase_abund.replace(0, np.nan))

metrics = pd.DataFrame(index=log_abund.index)
metrics['log_SD'] = log_abund.std(axis=1)                # variability metric
metrics['log_range'] = log_abund.max(axis=1) - log_abund.min(axis=1)  # max log-fold-change
metrics['mean_log_abund'] = log_abund.mean(axis=1)

# Compartment shift: max compartment mean - min compartment mean (in log space)
comp_means = pd.DataFrame({
    c: log_abund[[f'{c}_{t}' for t in timepoints]].mean(axis=1)
    for c in compartments
})
metrics['compartment_shift'] = comp_means.max(axis=1) - comp_means.min(axis=1)

# Temporal range within best compartment (max within-compartment temporal range)
temp_ranges = []
for c in compartments:
    cols = [f'{c}_{t}' for t in timepoints]
    temp_ranges.append(log_abund[cols].max(axis=1) - log_abund[cols].min(axis=1))
metrics['max_temporal_range'] = pd.concat(temp_ranges, axis=1).max(axis=1)

print("\nMetric correlations:")
print(metrics[['log_SD', 'log_range', 'compartment_shift', 'max_temporal_range']].corr().round(2))

# ============================================================
# Top 20 by log_SD (primary metric)
# ============================================================
top20 = metrics.nlargest(20, 'log_SD').index.tolist()
print(f"\nTop 20 kinases by log10-SD across 15 conditions:")
top_table = metrics.loc[top20, ['log_SD', 'log_range', 'compartment_shift', 'max_temporal_range', 'mean_log_abund']]
top_table['rank'] = range(1, 21)
print(top_table.round(3).to_string())

# Compare with manual 18-list
manual_18 = ['EGFR', 'MAPK1', 'MAPK3', 'AKT1', 'AKT2', 'MTOR', 'CDK1', 'CDK2', 'CDK9',
             'CSNK2A1', 'GSK3B', 'SRC', 'PRKCA', 'PRKCD', 'RPS6KB1', 'PLK1', 'CHEK1', 'AURKB']

overlap = set(top20) & set(manual_18)
new_in_top20 = set(top20) - set(manual_18)
manual_not_in_top20 = set(manual_18) - set(top20)

print(f"\nOverlap with manual 18-list: {len(overlap)}/18")
print(f"  Shared: {sorted(overlap)}")
print(f"  NEW in data-driven top20 (not in manual): {sorted(new_in_top20)}")
print(f"  In manual but NOT in top20: {sorted(manual_not_in_top20)}")

# Show where the manual 18 ranked
print(f"\nManual 18-list ranks by log_SD:")
manual_ranks = metrics.sort_values('log_SD', ascending=False).reset_index()
manual_ranks['rank'] = range(1, len(manual_ranks) + 1)
manual_ranks_subset = manual_ranks[manual_ranks['Gene'].isin(manual_18)][['Gene', 'rank', 'log_SD', 'compartment_shift']]
print(manual_ranks_subset.to_string(index=False))

# ============================================================
# Z-score and plot heatmap
# ============================================================
def zscore_row(row):
    vals = row.values.astype(float)
    mask = ~np.isnan(vals)
    if mask.sum() < 2:
        return row
    m = np.nanmean(vals)
    s = np.nanstd(vals)
    if s == 0 or np.isnan(s):
        return pd.Series(np.zeros_like(vals), index=row.index)
    return pd.Series((vals - m) / s, index=row.index)

top20_z = kinase_abund.loc[top20].apply(zscore_row, axis=1)

# Cluster rows by similarity for prettier heatmap
from scipy.cluster.hierarchy import linkage, leaves_list
filled = top20_z.fillna(0).values
link = linkage(filled, method='average', metric='euclidean')
order = leaves_list(link)
top20_z_sorted = top20_z.iloc[order]

# Build a side annotation with the metric values (rank label)
row_labels = []
for k in top20_z_sorted.index:
    rank = top20.index(k) + 1
    sd = metrics.loc[k, 'log_SD']
    cs = metrics.loc[k, 'compartment_shift']
    label = f"{k} (#{rank}, SD={sd:.2f}, Δcomp={cs:.2f})"
    row_labels.append(label)
top20_z_sorted_labeled = top20_z_sorted.copy()
top20_z_sorted_labeled.index = row_labels

# --- Plot ---
fig, ax = plt.subplots(figsize=(13, 9))
sns.heatmap(top20_z_sorted_labeled, cmap='RdBu_r', center=0,
            ax=ax, annot=True, fmt='.2f',
            cbar_kws={'label': 'Z-score (per kinase row)'},
            linewidths=0.5, linecolor='white',
            annot_kws={'size': 7})
ax.set_title('Top 20 Most Spatially-Temporally Dynamic Kinases\n'
             '(Data-driven: ranked by log10-SD across 15 conditions, hierarchically clustered)',
             fontsize=13)
ax.set_xlabel('Condition', fontsize=11)
ax.set_ylabel('')
# Compartment separators
for i in [5, 10]:
    ax.axvline(x=i, color='black', linewidth=2)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=9, rotation=0)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'step2_top20_data_driven_heatmap.png'), dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved: step2_top20_data_driven_heatmap.png")

# --- Bonus: top-20 by compartment_shift only ---
top20_cs = metrics.nlargest(20, 'compartment_shift').index.tolist()
print(f"\nTop 20 by compartment_shift (most spatially redistributed):")
print(metrics.loc[top20_cs, ['compartment_shift', 'log_SD', 'mean_log_abund']].round(3).to_string())

top20_cs_z = kinase_abund.loc[top20_cs].apply(zscore_row, axis=1)
filled_cs = top20_cs_z.fillna(0).values
link_cs = linkage(filled_cs, method='average', metric='euclidean')
order_cs = leaves_list(link_cs)
top20_cs_z_sorted = top20_cs_z.iloc[order_cs]
row_labels_cs = [f"{k} (Δcomp={metrics.loc[k, 'compartment_shift']:.2f})"
                  for k in top20_cs_z_sorted.index]
top20_cs_z_sorted.index = row_labels_cs

fig, ax = plt.subplots(figsize=(13, 9))
sns.heatmap(top20_cs_z_sorted, cmap='RdBu_r', center=0,
            ax=ax, annot=True, fmt='.2f',
            cbar_kws={'label': 'Z-score'},
            linewidths=0.5, linecolor='white',
            annot_kws={'size': 7})
ax.set_title('Top 20 Kinases by Compartment Shift\n'
             '(Largest difference between max and min compartment mean abundance)',
             fontsize=13)
ax.set_xlabel('Condition', fontsize=11)
for i in [5, 10]:
    ax.axvline(x=i, color='black', linewidth=2)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=9, rotation=0)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'step2_top20_compartment_shift_heatmap.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: step2_top20_compartment_shift_heatmap.png")

# --- Save metric table for reference ---
metrics_full = metrics.sort_values('log_SD', ascending=False)
metrics_full.to_csv(os.path.join(OUTDIR, 'kinase_dynamics_metrics.csv'))
print(f"Saved: kinase_dynamics_metrics.csv ({len(metrics_full)} kinases ranked)")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Total kinases analyzed: {len(kinase_abund)}")
print(f"Ranking metric: SD of log10(abundance) across 15 spatiotemporal conditions")
print(f"\nTop 20 by log_SD:")
for i, k in enumerate(top20, 1):
    in_manual = '★' if k in manual_18 else ' '
    print(f"  {in_manual} #{i:2d} {k:10s} (log_SD={metrics.loc[k,'log_SD']:.2f}, "
          f"Δcomp={metrics.loc[k,'compartment_shift']:.2f})")

print(f"\n★ = also in manual 18-list ({len(overlap)}/18 overlap)")
print(f"\nManual kinases NOT making top 20: {sorted(manual_not_in_top20)}")
print(f"   (These are present but more uniformly distributed across conditions.)")
