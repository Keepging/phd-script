"""
Check whether top multi-kinase hub substrates have high spatial variability
in protein abundance across compartments.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

OUTDIR = r'D:\博士\Protein contour\Phospho\multi_kinase_analysis'

def save_tiff(fig, filename, dpi=300):
    path = filename if filename.endswith('.tiff') else filename.rsplit('.', 1)[0] + '.tiff'
    fig.savefig(path, dpi=dpi, bbox_inches='tight', format='tiff',
                pil_kwargs={'compression': 'tiff_lzw'})

# ============================================================
# Load data
# ============================================================
pam = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\protein_abundance_matrix.csv')
print(f'protein_abundance_matrix: {pam.shape}')

t3 = pd.read_csv(os.path.join(OUTDIR, 'task3_multi_kinase_substrates.csv'))
print(f'task3_multi_kinase_substrates: {t3.shape}')

# ============================================================
# Step 1: Compute per-protein spatial and temporal SD
# ============================================================
cyt_cols = [c for c in pam.columns if c.startswith('Cyt_')]
mem_cols = [c for c in pam.columns if c.startswith('Mem_')]
nuc_cols = [c for c in pam.columns if c.startswith('Nuc_')]
all_cond = cyt_cols + mem_cols + nuc_cols
print(f'Condition columns: {len(cyt_cols)} Cyt + {len(mem_cols)} Mem + {len(nuc_cols)} Nuc = {len(all_cond)}')

# Convert to numeric
for c in all_cond:
    pam[c] = pd.to_numeric(pam[c], errors='coerce')

pam['mean_Cyt'] = pam[cyt_cols].mean(axis=1)
pam['mean_Mem'] = pam[mem_cols].mean(axis=1)
pam['mean_Nuc'] = pam[nuc_cols].mean(axis=1)

# cross_compartment_SD: SD of the 3 compartment means
pam['cross_compartment_SD'] = pam[['mean_Cyt', 'mean_Mem', 'mean_Nuc']].std(axis=1, ddof=1)

# temporal_SD: SD across all 15 conditions
pam['temporal_SD'] = pam[all_cond].std(axis=1, ddof=1)

# max_compartment
def get_max_comp(row):
    vals = {'Cyt': row['mean_Cyt'], 'Mem': row['mean_Mem'], 'Nuc': row['mean_Nuc']}
    valid = {k: v for k, v in vals.items() if not np.isnan(v)}
    if not valid:
        return np.nan
    return max(valid, key=valid.get)

pam['max_compartment'] = pam.apply(get_max_comp, axis=1)

# Drop rows where cross_compartment_SD is NaN (need at least 2 compartments)
valid_pam = pam.dropna(subset=['cross_compartment_SD']).copy()
print(f'Proteins with valid cross_compartment_SD: {len(valid_pam)}')

# ============================================================
# Step 2: Top 20 multi-kinase hubs
# ============================================================
top20_hubs = t3.nlargest(20, 'N_competing_kinases')
print(f'\nTop 20 multi-kinase hubs:')
for _, row in top20_hubs.iterrows():
    print(f'  {row["Substrate"]:12s}  N_kinases={row["N_competing_kinases"]}  '
          f'Top_r={row["Top_r"]:.4f}  Spread={row["Spread"]:.4f}')

# ============================================================
# Step 3: Look up hub proteins in abundance matrix
# ============================================================
hub_genes = top20_hubs['Substrate'].tolist()

# Compute percentile rank for all proteins
valid_pam['SD_percentile_rank'] = valid_pam['cross_compartment_SD'].rank(pct=True) * 100

hub_data = []
for gene in hub_genes:
    match = valid_pam[valid_pam['Gene'] == gene]
    if len(match) == 0:
        hub_data.append({
            'Substrate': gene,
            'N_competing_kinases': int(top20_hubs[top20_hubs['Substrate'] == gene]['N_competing_kinases'].values[0]),
            'cross_compartment_SD': np.nan,
            'temporal_SD': np.nan,
            'max_compartment': np.nan,
            'SD_percentile_rank': np.nan,
            'mean_Cyt': np.nan,
            'mean_Mem': np.nan,
            'mean_Nuc': np.nan,
        })
    else:
        row = match.iloc[0]
        hub_data.append({
            'Substrate': gene,
            'N_competing_kinases': int(top20_hubs[top20_hubs['Substrate'] == gene]['N_competing_kinases'].values[0]),
            'cross_compartment_SD': round(row['cross_compartment_SD'], 2),
            'temporal_SD': round(row['temporal_SD'], 2),
            'max_compartment': row['max_compartment'],
            'SD_percentile_rank': round(row['SD_percentile_rank'], 1),
            'mean_Cyt': round(row['mean_Cyt'], 1),
            'mean_Mem': round(row['mean_Mem'], 1),
            'mean_Nuc': round(row['mean_Nuc'], 1),
        })

hub_df = pd.DataFrame(hub_data)
hub_path = os.path.join(OUTDIR, 'hub_spatial_variability.csv')
hub_df.to_csv(hub_path, index=False)
print(f'\nSaved: {hub_path}')

# ============================================================
# Step 4: Compare with proteome distribution
# ============================================================
proteome_sd = valid_pam['cross_compartment_SD'].dropna()
hub_sd = hub_df['cross_compartment_SD'].dropna()

print(f'\n=== CROSS-COMPARTMENT SD COMPARISON ===')
print(f'Proteome (n={len(proteome_sd)}):')
print(f'  Median: {proteome_sd.median():.2f}')
print(f'  Mean:   {proteome_sd.mean():.2f}')
print(f'  Q25:    {proteome_sd.quantile(0.25):.2f}')
print(f'  Q75:    {proteome_sd.quantile(0.75):.2f}')
print(f'  Q90:    {proteome_sd.quantile(0.90):.2f}')

print(f'\nTop 20 hubs (n={len(hub_sd)}):')
print(f'  Median: {hub_sd.median():.2f}')
print(f'  Mean:   {hub_sd.mean():.2f}')

n_top25 = int((hub_df['SD_percentile_rank'] >= 75).sum())
n_top10 = int((hub_df['SD_percentile_rank'] >= 90).sum())
n_found = int(hub_df['cross_compartment_SD'].notna().sum())
print(f'\nHub percentile distribution (out of {n_found} found):')
print(f'  In top 25% (pctl >= 75): {n_top25}/{n_found}')
print(f'  In top 10% (pctl >= 90): {n_top10}/{n_found}')

print(f'\nPer-hub details:')
for _, row in hub_df.iterrows():
    if np.isnan(row['cross_compartment_SD']):
        print(f'  {row["Substrate"]:12s}  N_kin={int(row["N_competing_kinases"])}  NOT FOUND')
    else:
        print(f'  {row["Substrate"]:12s}  N_kin={int(row["N_competing_kinases"])}  '
              f'SD={row["cross_compartment_SD"]:10.2f}  pctl={row["SD_percentile_rank"]:5.1f}%  '
              f'max={row["max_compartment"]:3s}  '
              f'Cyt={row["mean_Cyt"]:10.1f} Mem={row["mean_Mem"]:10.1f} Nuc={row["mean_Nuc"]:10.1f}')

# ============================================================
# Figure: Proteome SD distribution + hub overlay
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [2, 1]})

# Left panel: violin + strip of hubs
ax = axes[0]
parts = ax.violinplot(proteome_sd.values, positions=[1], showmedians=True,
                      showextrema=False, widths=0.6)
for pc in parts['bodies']:
    pc.set_facecolor('#3498db')
    pc.set_alpha(0.4)
parts['cmedians'].set_color('#2c3e50')
parts['cmedians'].set_linewidth(2)

# Overlay hub points
hub_valid = hub_df.dropna(subset=['cross_compartment_SD'])
jitter = np.random.default_rng(42).uniform(-0.08, 0.08, size=len(hub_valid))
ax.scatter(np.ones(len(hub_valid)) + jitter, hub_valid['cross_compartment_SD'],
           c='#e74c3c', s=60, zorder=5, edgecolors='white', linewidth=0.5,
           label=f'Top 20 hubs (n={len(hub_valid)})')

# Annotate top 5 hubs
for _, row in hub_valid.nlargest(5, 'cross_compartment_SD').iterrows():
    ax.annotate(row['Substrate'],
                xy=(1, row['cross_compartment_SD']),
                xytext=(1.35, row['cross_compartment_SD']),
                fontsize=7, color='#c0392b',
                arrowprops=dict(arrowstyle='-', color='#c0392b', lw=0.5))

ax.axhline(proteome_sd.median(), color='#2c3e50', ls='--', lw=1, alpha=0.5,
           label=f'Proteome median ({proteome_sd.median():.0f})')
ax.axhline(proteome_sd.quantile(0.75), color='#7f8c8d', ls=':', lw=1, alpha=0.5,
           label=f'Proteome Q75 ({proteome_sd.quantile(0.75):.0f})')

ax.set_ylabel('Cross-compartment SD (protein abundance)', fontsize=11)
ax.set_xticks([1])
ax.set_xticklabels(['Full proteome\n+ hub overlay'])
ax.legend(fontsize=9, loc='upper left')
ax.set_title('Spatial Variability: Multi-kinase Hubs vs Proteome', fontsize=12)

# Right panel: ranked hub SD with percentile coloring
ax2 = axes[1]
hub_sorted = hub_valid.sort_values('cross_compartment_SD', ascending=True).reset_index(drop=True)
colors = ['#e74c3c' if p >= 75 else '#f39c12' if p >= 50 else '#3498db'
          for p in hub_sorted['SD_percentile_rank']]
ax2.barh(range(len(hub_sorted)), hub_sorted['cross_compartment_SD'],
         color=colors, edgecolor='white', linewidth=0.4)
ax2.set_yticks(range(len(hub_sorted)))
ax2.set_yticklabels([f"{r['Substrate']} ({int(r['N_competing_kinases'])}K)"
                     for _, r in hub_sorted.iterrows()], fontsize=8)
ax2.set_xlabel('Cross-compartment SD', fontsize=10)
ax2.set_title('Hub Proteins Ranked\n(red=top25%, orange=25-50%)', fontsize=10)
ax2.axvline(proteome_sd.median(), color='#2c3e50', ls='--', lw=1, alpha=0.5)
ax2.axvline(proteome_sd.quantile(0.75), color='#7f8c8d', ls=':', lw=1, alpha=0.5)

plt.tight_layout()
fig_path = os.path.join(OUTDIR, 'hub_vs_proteome_sd_boxplot.png')
fig.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
save_tiff(fig, fig_path)
plt.close(fig)
print(f'\nSaved figure: {fig_path} + .tiff')
print(f'All outputs in: {OUTDIR}')
