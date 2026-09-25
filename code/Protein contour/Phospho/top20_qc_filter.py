"""
QC the data-driven top 20 kinases:
- Detection frequency (non-NaN conditions / 15)
- Absolute intensity level
- Per-compartment detection (does it have data in all 3 compartments?)
- Coefficient of variation on log scale (relative dynamics)
- Flag "pseudo-dynamic" kinases driven by sparse detection
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
metrics = pd.read_csv(os.path.join(OUTDIR, 'kinase_dynamics_metrics.csv'), index_col='Gene')

# Top 20 by log_SD (the original list)
top20 = metrics.sort_values('log_SD', ascending=False).head(20).index.tolist()

# ============================================================
# QC each top-20 kinase
# ============================================================
print("=" * 100)
print("QC of Top 20 Data-Driven Dynamic Kinases")
print("=" * 100)

qc_rows = []
for k in top20:
    row = abund.loc[k, cond_order]
    n_detected = int(row.notna().sum())
    detection_pct = n_detected / 15 * 100

    # Per-compartment detection
    cyt_det = int(row[[f'Cyt_{t}' for t in timepoints]].notna().sum())
    mem_det = int(row[[f'Mem_{t}' for t in timepoints]].notna().sum())
    nuc_det = int(row[[f'Nuc_{t}' for t in timepoints]].notna().sum())
    n_comp_with_data = sum([1 if x > 0 else 0 for x in [cyt_det, mem_det, nuc_det]])

    vals = row.dropna().values.astype(float)
    mean_int = np.mean(vals)
    median_int = np.median(vals)
    log_vals = np.log10(vals[vals > 0])
    log_sd = np.std(log_vals)
    log_range = np.max(log_vals) - np.min(log_vals)

    qc_rows.append({
        'Kinase': k,
        'rank_log_SD': metrics.loc[k, 'log_SD'],
        'n_detected': int(n_detected),
        'detection_pct': detection_pct,
        'Cyt_n': int(cyt_det),
        'Mem_n': int(mem_det),
        'Nuc_n': int(nuc_det),
        'n_compartments': int(n_comp_with_data),
        'mean_intensity': mean_int,
        'median_intensity': median_int,
        'log_SD': log_sd,
        'log_range': log_range,
    })

qc_df = pd.DataFrame(qc_rows)

# Determine median intensity of all kinases for context
all_kinases = metrics.index
all_intensities = []
for k in all_kinases:
    if k in abund.index:
        vals = abund.loc[k, cond_order].dropna().values
        if len(vals) > 0:
            all_intensities.append(np.mean(vals))
all_int_arr = np.array(all_intensities)
median_all = np.median(all_int_arr)
q25_all = np.percentile(all_int_arr, 25)
print(f"\nReference: All-kinase mean-intensity median = {median_all:,.0f}")
print(f"           25th percentile = {q25_all:,.0f}")

# ============================================================
# Print QC table
# ============================================================
print("\n" + "=" * 100)
print(f"{'Kinase':<10s} {'#det':>5s} {'%det':>5s} {'Cyt':>4s} {'Mem':>4s} {'Nuc':>4s} "
      f"{'#comp':>5s} {'mean_int':>12s} {'log_SD':>7s} {'log_rng':>7s} {'flag':<25s}")
print("-" * 100)

for _, r in qc_df.iterrows():
    flags = []
    # Hard filter: < 8/15 → reject (sparse detection)
    if r['n_detected'] < 8:
        flags.append('SPARSE')
    elif r['n_detected'] < 12:
        flags.append('partial')
    # Single-compartment dominance
    if r['n_compartments'] < 3:
        flags.append('NOT-3COMP')
    # Low intensity (close to detection limit)
    if r['mean_intensity'] < q25_all:
        flags.append('LOW-INT')
    # Possibly fake dynamics: very few detections
    if r['n_detected'] <= 5:
        flags.append('FAKE-DYN!')

    flag_str = '/'.join(flags) if flags else 'OK'
    print(f"{r['Kinase']:<10s} {r['n_detected']:>5d} {r['detection_pct']:>4.0f}% "
          f"{r['Cyt_n']:>4d} {r['Mem_n']:>4d} {r['Nuc_n']:>4d} "
          f"{r['n_compartments']:>5d} {r['mean_intensity']:>12,.0f} "
          f"{r['log_SD']:>7.3f} {r['log_range']:>7.3f} {flag_str:<25s}")

# Save QC table
qc_df.to_csv(os.path.join(OUTDIR, 'top20_qc_table.csv'), index=False)
print(f"\nSaved: top20_qc_table.csv")

# ============================================================
# Apply filters and produce CLEAN top 20
# ============================================================
print("\n" + "=" * 100)
print("FILTER LOGIC")
print("=" * 100)
print("Reject if any of:")
print("  - n_detected < 10  (must be measured in >=2/3 of conditions; filters 3-4 cond pseudo-dynamics)")
print("  - n_compartments < 2  (must have data in at least 2 compartments)")
print("  - mean_intensity < 25th percentile of all kinases (close to detection limit)")
print()

rejected = qc_df[(qc_df['n_detected'] < 10) |
                 (qc_df['n_compartments'] < 2) |
                 (qc_df['mean_intensity'] < q25_all)]
print(f"Rejected from top 20: {len(rejected)}")
if len(rejected) > 0:
    print(rejected[['Kinase', 'n_detected', 'n_compartments', 'mean_intensity']].to_string(index=False))

passed_top20 = qc_df[~qc_df['Kinase'].isin(rejected['Kinase'])]
print(f"\nPassed top 20: {len(passed_top20)}")

# ============================================================
# Find replacements: extend list with kinases ranked >20 that pass QC
# ============================================================
extended_pool = metrics.sort_values('log_SD', ascending=False).index.tolist()

# Recompute QC for the entire pool
def compute_qc(k):
    if k not in abund.index:
        return None
    row = abund.loc[k, cond_order]
    n_det = int(row.notna().sum())
    cyt_n = int(row[[f'Cyt_{t}' for t in timepoints]].notna().sum())
    mem_n = int(row[[f'Mem_{t}' for t in timepoints]].notna().sum())
    nuc_n = int(row[[f'Nuc_{t}' for t in timepoints]].notna().sum())
    n_comp = sum([1 if x > 0 else 0 for x in [cyt_n, mem_n, nuc_n]])
    vals = row.dropna().values.astype(float)
    mean_int = np.mean(vals) if len(vals) > 0 else np.nan
    return n_det, n_comp, mean_int

# Build clean top 20 by rank order, filter as we go
clean_top20 = []
for k in extended_pool:
    if len(clean_top20) >= 20:
        break
    qc = compute_qc(k)
    if qc is None:
        continue
    n_det, n_comp, mean_int = qc
    if n_det >= 10 and n_comp >= 2 and mean_int >= q25_all:
        clean_top20.append(k)

new_in_clean = set(clean_top20) - set(top20)
removed = set(top20) - set(clean_top20)
print(f"\n--- CLEAN top 20 (filtered) ---")
for i, k in enumerate(clean_top20, 1):
    is_new = '★' if k in new_in_clean else ' '
    n_det, n_comp, mean_int = compute_qc(k)
    log_sd = metrics.loc[k, 'log_SD']
    print(f"  {is_new} #{i:2d} {k:10s} (log_SD={log_sd:.2f}, det={n_det}/15, mean_int={mean_int:,.0f})")
print(f"\n★ = newly added (replaced sparse-detection kinase)")
print(f"Removed from original top 20: {sorted(removed)}")
print(f"Newly added: {sorted(new_in_clean)}")

# ============================================================
# Plot CLEAN top 20 heatmap
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

clean_z = abund.loc[clean_top20, cond_order].apply(zscore_row, axis=1)
filled = clean_z.fillna(0).values
from scipy.cluster.hierarchy import linkage, leaves_list
link = linkage(filled, method='average', metric='euclidean')
order = leaves_list(link)
clean_z_sorted = clean_z.iloc[order]

row_labels = []
for k in clean_z_sorted.index:
    rank = clean_top20.index(k) + 1
    sd = metrics.loc[k, 'log_SD']
    n_det = compute_qc(k)[0]
    label = f"{k} (#{rank}, SD={sd:.2f}, det={n_det}/15)"
    row_labels.append(label)
clean_z_sorted_labeled = clean_z_sorted.copy()
clean_z_sorted_labeled.index = row_labels

fig, ax = plt.subplots(figsize=(13, 9))
sns.heatmap(clean_z_sorted_labeled, cmap='RdBu_r', center=0,
            ax=ax, annot=True, fmt='.2f',
            cbar_kws={'label': 'Z-score (per kinase row)'},
            linewidths=0.5, linecolor='white',
            annot_kws={'size': 7})
ax.set_title('Top 20 Most Dynamic Kinases — QC-FILTERED\n'
             '(Detection >=10/15, present in >=2 compartments, intensity >25th percentile)',
             fontsize=13)
ax.set_xlabel('Condition', fontsize=11)
ax.set_ylabel('')
for i in [5, 10]:
    ax.axvline(x=i, color='black', linewidth=2)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=9, rotation=0)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'step2_top20_clean_heatmap.png'), dpi=200, bbox_inches='tight')
plt.close()
print(f"\nSaved: step2_top20_clean_heatmap.png")

# ============================================================
# Visualization: detection vs log_SD scatter (showing the trap)
# ============================================================
all_qc = []
for k in metrics.index:
    qc = compute_qc(k)
    if qc is None:
        continue
    n_det, n_comp, mean_int = qc
    all_qc.append({
        'Kinase': k,
        'n_detected': n_det,
        'mean_intensity': mean_int,
        'log_SD': metrics.loc[k, 'log_SD'],
        'in_orig_top20': k in top20,
        'in_clean_top20': k in clean_top20,
    })
all_qc_df = pd.DataFrame(all_qc)

fig, ax = plt.subplots(figsize=(12, 7))
# Background: all kinases
not_top = all_qc_df[~all_qc_df['in_orig_top20'] & ~all_qc_df['in_clean_top20']]
ax.scatter(not_top['n_detected'], not_top['log_SD'],
           c='#cccccc', s=25, alpha=0.5, label=f'Other ({len(not_top)})')

# Original top 20 (will be split into kept vs rejected)
orig_kept = all_qc_df[all_qc_df['in_orig_top20'] & all_qc_df['in_clean_top20']]
orig_rejected = all_qc_df[all_qc_df['in_orig_top20'] & ~all_qc_df['in_clean_top20']]
new_added = all_qc_df[~all_qc_df['in_orig_top20'] & all_qc_df['in_clean_top20']]

ax.scatter(orig_kept['n_detected'], orig_kept['log_SD'],
           c='#27ae60', s=120, alpha=0.85, edgecolor='black',
           label=f'Original top 20 - PASS QC ({len(orig_kept)})')
ax.scatter(orig_rejected['n_detected'], orig_rejected['log_SD'],
           c='#e74c3c', s=140, alpha=0.85, edgecolor='black', marker='X',
           label=f'Original top 20 - REJECTED ({len(orig_rejected)})')
ax.scatter(new_added['n_detected'], new_added['log_SD'],
           c='#3498db', s=120, alpha=0.85, edgecolor='black', marker='s',
           label=f'Newly promoted into top 20 ({len(new_added)})')

# Annotate rejected and newly added
for _, r in orig_rejected.iterrows():
    ax.annotate(r['Kinase'], (r['n_detected'], r['log_SD']),
                fontsize=9, fontweight='bold', color='#922b21',
                xytext=(5, 5), textcoords='offset points')
for _, r in new_added.iterrows():
    ax.annotate(r['Kinase'], (r['n_detected'], r['log_SD']),
                fontsize=9, color='#1f5582',
                xytext=(5, -10), textcoords='offset points')

# Threshold line
ax.axvline(x=9.5, color='red', linestyle='--', linewidth=1.5, alpha=0.7,
           label='Detection threshold (10/15)')

ax.set_xlabel('Detection frequency (conditions with non-NaN, out of 15)', fontsize=12)
ax.set_ylabel('log10-SD of abundance (dynamics)', fontsize=12)
ax.set_title('QC of Top-20 Dynamic Kinases — Detection vs Apparent Dynamics\n'
             '(Sparsely-detected kinases can show inflated SD)',
             fontsize=13)
ax.legend(fontsize=10, loc='lower left')
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'step2_top20_qc_scatter.png'), dpi=200, bbox_inches='tight')
plt.close()
print(f"Saved: step2_top20_qc_scatter.png")

# ============================================================
# Final summary
# ============================================================
print("\n" + "=" * 100)
print("FINAL SUMMARY")
print("=" * 100)
print(f"Original top 20 by log_SD:")
print(f"  Passed QC: {len(orig_kept)}")
print(f"  Rejected:  {len(orig_rejected)} ({sorted(removed)})")
print(f"\nClean top 20 (after filter + replacement):")
for i, k in enumerate(clean_top20, 1):
    is_new = '★' if k in new_in_clean else ' '
    n_det, _, mi = compute_qc(k)
    print(f"  {is_new} #{i:2d} {k:10s} det={n_det}/15, mean_int={mi:,.0f}")
