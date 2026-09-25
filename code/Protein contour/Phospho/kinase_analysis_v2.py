"""
Kinase-substrate mapping analysis v2 — fixed classification + raw detection matrix
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency
import os
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.family'] = 'sans-serif'

outdir = r'D:\博士\Protein contour\Phospho\kinase_results\v2'
os.makedirs(outdir, exist_ok=True)

# ============================================================
# Load matched data from v1
# ============================================================
merged = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_results\kinase_matched_sites.csv')
print(f"Loaded kinase_matched_sites: {merged.shape}")
print(f"Unique sites: {merged['PTM_collapse_key'].nunique()}")
print(f"Unique kinases: {merged['Kinase Name'].nunique()}")

# Save old classification for comparison
old_other_pct = (merged['kinase_family'] == 'Other').mean() * 100

# ============================================================
# Fix 1: Reclassify kinase families using Kinase Name
# ============================================================
print("\n" + "=" * 60)
print("FIX 1: Reclassify kinase families")
print("=" * 60)

family_rules = {
    'MAPK cascade': ['MAPK1', 'MAPK3', 'MAPK8', 'MAPK9', 'MAPK14',
                     'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K7',
                     'MAPKAPK2', 'MAPKAPK3', 'MAPKAPK5',
                     'RPS6KA1', 'RPS6KA2', 'RPS6KA3', 'RPS6KA6'],
    'CDK': ['CDK1', 'CDK2', 'CDK3', 'CDK4', 'CDK5', 'CDK6', 'CDK7',
            'CDK9', 'CDK12', 'CDK13', 'CDK18'],
    'CK1/CK2': ['CSNK1A1', 'CSNK1D', 'CSNK1E', 'CSNK1G1', 'CSNK1G2', 'CSNK1G3',
                 'CSNK2A1', 'CSNK2A2', 'CSNK2B'],
    'PKA/PKC': ['PRKACA', 'PRKACB', 'PRKCA', 'PRKCB', 'PRKCD', 'PRKCE',
                'PRKCI', 'PRKCZ', 'PRKD1', 'PRKD2', 'SGK1', 'SGK3'],
    'AKT/mTOR': ['AKT1', 'AKT2', 'AKT3', 'MTOR', 'RPS6KB1', 'RPS6KB2', 'PDPK1'],
    'RTK/Src': ['EGFR', 'ERBB2', 'ERBB3', 'MET', 'IGF1R', 'INSR', 'FGFR1',
                'SRC', 'FYN', 'YES1', 'LYN', 'ABL1', 'ABL2'],
    'Cell cycle': ['PLK1', 'PLK2', 'PLK3', 'AURKA', 'AURKB', 'AURKC',
                   'CHEK1', 'CHEK2', 'ATM', 'ATR', 'NEK2', 'NEK6', 'NEK7',
                   'LATS1', 'LATS2'],
    'GSK3': ['GSK3A', 'GSK3B'],
    'CAMK': ['CAMK1', 'CAMK1D', 'CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G',
             'CAMKK1', 'CAMKK2', 'DAPK1', 'DAPK3'],
    'AMPK': ['PRKAA1', 'PRKAA2'],
    'CLK/DYRK': ['CLK1', 'CLK2', 'CLK3', 'CLK4', 'DYRK1A', 'DYRK1B', 'DYRK2'],
}

# Build reverse map (kinase_name -> family), case-insensitive
kinase_to_family = {}
for fam, members in family_rules.items():
    for m in members:
        kinase_to_family[m.upper()] = fam

# Apply classification
merged['kinase_family'] = merged['Kinase Name'].str.upper().map(kinase_to_family).fillna('Other')

new_other_pct = (merged['kinase_family'] == 'Other').mean() * 100

print(f"\nOther % before fix: {old_other_pct:.1f}%")
print(f"Other % after fix:  {new_other_pct:.1f}%")

print(f"\nKinase family counts (by row, descending):")
fam_counts = merged['kinase_family'].value_counts()
for fam, cnt in fam_counts.items():
    n_sites = merged[merged['kinase_family'] == fam]['PTM_collapse_key'].nunique()
    print(f"  {fam:15s}: {cnt:5d} rows, {n_sites:4d} unique sites")

print(f"\nTop 10 kinases still in 'Other':")
other_kinases = merged[merged['kinase_family'] == 'Other']['Kinase Name'].value_counts().head(10)
print(other_kinases)

# ============================================================
# Fix 2: Task 4 — Temporal kinase family composition (High mobility)
# ============================================================
print("\n" + "=" * 60)
print("FIX 2: Temporal kinase family heatmap (High mobility)")
print("=" * 60)

high_m = merged[merged['MobilityLevel'] == 'High'].copy()
print(f"High-mobility matched rows: {len(high_m)}")

# Deduplicate: one (site, kinase_family) pair per row
high_site_fam = high_m.drop_duplicates(subset=['PTM_collapse_key', 'kinase_family'])

timepoints = ['2min', '8min', '20min', '90min']
high_valid = high_site_fam[high_site_fam['movement_best_time'].isin(timepoints)]

print(f"\nPer-timepoint breakdown:")
for tp in timepoints:
    subset = high_valid[high_valid['movement_best_time'] == tp]
    n = subset['PTM_collapse_key'].nunique()
    top = subset['kinase_family'].value_counts()
    top_fam = top.index[0] if len(top) > 0 else 'N/A'
    print(f"  {tp:5s}: n={n:4d} sites, top = {top_fam} ({top.iloc[0] if len(top)>0 else 0})")

# Cross-tabulation
ct = pd.crosstab(high_valid['movement_best_time'], high_valid['kinase_family'])
ct = ct.reindex(timepoints).fillna(0)
ct_prop = ct.div(ct.sum(axis=1), axis=0)

# Order families: by total count descending, Other last
fam_order = ct.sum().sort_values(ascending=False).index.tolist()
if 'Other' in fam_order:
    fam_order.remove('Other')
    fam_order.append('Other')

ct_prop = ct_prop[fam_order]
ct_ordered = ct[fam_order]

print("\nProportion table:")
print(ct_prop.round(3).to_string())

# --- Color palette ---
n_fam = len(fam_order)
palette = sns.color_palette("tab20", n_fam)
color_map = dict(zip(fam_order, palette))

# --- Fig A: Stacked bar chart ---
fig, ax = plt.subplots(figsize=(10, 6))
bottom = np.zeros(len(timepoints))
for fam in fam_order:
    vals = ct_prop[fam].values
    bars = ax.bar(timepoints, vals, bottom=bottom, label=fam, color=color_map[fam], edgecolor='white', linewidth=0.5)
    bottom += vals

# Annotate n= on top of each bar
for i, tp in enumerate(timepoints):
    n = int(ct_ordered.loc[tp].sum())
    ax.text(i, 1.02, f'n={n}', ha='center', va='bottom', fontsize=10, fontweight='bold')

ax.set_xlabel('Movement Best Time', fontsize=12)
ax.set_ylabel('Proportion', fontsize=12)
ax.set_title('Kinase Family Composition by Movement Timepoint\n(High-mobility matched sites, v2)', fontsize=13)
ax.set_ylim(0, 1.12)
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9, framealpha=0.9)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'fix2_stacked_bar.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: fix2_stacked_bar.png")

# --- Fig B: Heatmap ---
fig, ax = plt.subplots(figsize=(10, 8))
hm_data = ct_prop[fam_order].T
sns.heatmap(hm_data, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax,
            linewidths=0.5, vmin=0, vmax=hm_data.max().max() * 1.1)
ax.set_title('Kinase Family Proportion Heatmap\n(High-mobility matched sites, v2)', fontsize=13)
ax.set_xlabel('Movement Best Time', fontsize=12)
ax.set_ylabel('Kinase Family', fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'fix2_heatmap.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: fix2_heatmap.png")

# ============================================================
# Fix 3: Task 5 — High vs Low mobility comparison
# ============================================================
print("\n" + "=" * 60)
print("FIX 3: High vs Low mobility kinase comparison")
print("=" * 60)

site_fam = merged.drop_duplicates(subset=['PTM_collapse_key', 'kinase_family'])
site_fam_hl = site_fam[site_fam['MobilityLevel'].isin(['High', 'Low'])]

ct_hl = pd.crosstab(site_fam_hl['MobilityLevel'], site_fam_hl['kinase_family'])
# Reorder columns
cols_present = [f for f in fam_order if f in ct_hl.columns]
ct_hl = ct_hl[cols_present]
ct_hl_prop = ct_hl.div(ct_hl.sum(axis=1), axis=0)

n_high = ct_hl.loc['High'].sum() if 'High' in ct_hl.index else 0
n_low = ct_hl.loc['Low'].sum() if 'Low' in ct_hl.index else 0

print(f"High mobility: {n_high} (site, family) pairs")
print(f"Low mobility:  {n_low} (site, family) pairs")
print("\nCounts:")
print(ct_hl)
print("\nProportions:")
print(ct_hl_prop.round(3).to_string())

chi2, p_val, dof, expected = chi2_contingency(ct_hl)
print(f"\nChi-square: chi2={chi2:.2f}, p={p_val:.2e}, dof={dof}")

# Grouped bar chart
fig, ax = plt.subplots(figsize=(12, 6))
x = np.arange(len(cols_present))
width = 0.35

h_vals = ct_hl_prop.loc['High'].values if 'High' in ct_hl_prop.index else np.zeros(len(x))
l_vals = ct_hl_prop.loc['Low'].values if 'Low' in ct_hl_prop.index else np.zeros(len(x))

bars_h = ax.bar(x - width/2, h_vals, width, label=f'High Mobility (n={int(n_high)})',
                color='#e74c3c', alpha=0.85, edgecolor='white')
bars_l = ax.bar(x + width/2, l_vals, width, label=f'Low Mobility (n={int(n_low)})',
                color='#3498db', alpha=0.85, edgecolor='white')

ax.set_xlabel('Kinase Family', fontsize=12)
ax.set_ylabel('Proportion', fontsize=12)
ax.set_title(f'Kinase Family Distribution: High vs Low Mobility\n(Chi-square p = {p_val:.2e})', fontsize=13)
ax.set_xticks(x)
ax.set_xticklabels(cols_present, rotation=45, ha='right', fontsize=10)
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'fix3_high_vs_low.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: fix3_high_vs_low.png")

# ============================================================
# Fix 4: Task 6 — Detection matrix using binary columns
# ============================================================
print("\n" + "=" * 60)
print("FIX 4: Detection matrix (binary columns from success_df)")
print("=" * 60)

# success_df has binary 0/1 columns: 0min_Cyt, 0min_Mem, ..., 90min_Nuc
# These are TRUE detection flags (0 = not detected, 1 = detected)
success_df = pd.read_csv(r'D:\博士\Protein contour\Phospho\success_df_M1M2_cleaned.csv')

# For each timepoint, detected = max across compartments
all_tps = ['0min', '2min', '8min', '20min', '90min']
tp_bin_cols = {
    '0min': ['0min_Cyt', '0min_Mem', '0min_Nuc'],
    '2min': ['2min_Cyt', '2min_Mem', '2min_Nuc'],
    '8min': ['8min_Cyt', '8min_Mem', '8min_Nuc'],
    '20min': ['20min_Cyt', '20min_Mem', '20min_Nuc'],
    '90min': ['90min_Cyt', '90min_Mem', '90min_Nuc'],
}
for tp, cols in tp_bin_cols.items():
    success_df[f'detected_{tp}'] = success_df[cols].max(axis=1).astype(int)

# Check: are there any 0s in detected columns?
for tp in all_tps:
    n_zero = (success_df[f'detected_{tp}'] == 0).sum()
    print(f"  detected_{tp}: {n_zero} zeros out of {len(success_df)} ({n_zero/len(success_df)*100:.1f}%)")

any_zeros = any((success_df[f'detected_{tp}'] == 0).sum() > 0 for tp in all_tps)

if not any_zeros:
    print("\n WARNING: All detected columns are 1 (fully imputed). Trying per-compartment level instead.")
    # Use per-compartment binary columns directly for more granular view
    use_compartment_level = True
else:
    use_compartment_level = False
    print("\n Binary detection columns have real 0s — using max-across-compartments.")

# High-mobility matched sites, proteins with >= 2 matched sites
high_m_sites = high_m['PTM_collapse_key'].unique()
# Get detection info from success_df for these sites
det_cols = [f'detected_{tp}' for tp in all_tps]
success_high = success_df[success_df['PTM_collapse_key'].isin(high_m_sites)].copy()

# Count sites per protein
protein_counts = success_high.groupby('Gene')['PTM_collapse_key'].nunique().sort_values(ascending=False)
proteins_2plus = protein_counts[protein_counts >= 2]
print(f"\nProteins with >= 2 high-mob matched sites: {len(proteins_2plus)}")

# Look for proteins with temporal variation
protein_variation = {}
for gene in proteins_2plus.index:
    gene_data = success_high[success_high['Gene'] == gene]
    det_matrix = gene_data[det_cols].values
    # variation = not all 1s and not all 0s in the matrix
    has_variation = det_matrix.min() != det_matrix.max()
    protein_variation[gene] = has_variation

proteins_with_var = [g for g, v in protein_variation.items() if v]
proteins_without_var = [g for g, v in protein_variation.items() if not v]
print(f"Proteins with temporal variation: {len(proteins_with_var)}")
print(f"Proteins without variation (all 1s): {len(proteins_without_var)}")

if len(proteins_with_var) == 0:
    print("\n WARNING: No proteins with temporal variation found at site-level (max across compartments).")
    print("This confirms detection is binary by site design. Switching to per-compartment detection matrix.")
    use_compartment_level = True

# Select top 6 proteins to display
if use_compartment_level:
    # Use per-compartment binary columns for richer variation
    comp_cols = []
    comp_labels = []
    for tp in all_tps:
        for comp in ['Cyt', 'Mem', 'Nuc']:
            col = f'{tp}_{comp}'
            comp_cols.append(col)
            comp_labels.append(f'{tp}\n{comp}')

    # Find proteins with per-compartment variation
    protein_comp_var = {}
    for gene in proteins_2plus.index:
        gene_data = success_high[success_high['Gene'] == gene]
        mat = gene_data[comp_cols].values
        protein_comp_var[gene] = mat.min() != mat.max()

    candidates = [g for g, v in protein_comp_var.items() if v]
    print(f"Proteins with per-compartment variation: {len(candidates)}")

    # Prefer proteins with more sites
    candidates_sorted = sorted(candidates, key=lambda g: proteins_2plus.get(g, 0), reverse=True)
    top_proteins = candidates_sorted[:6]

    if len(top_proteins) == 0:
        print("\n ERROR: No proteins with any variation. Data is fully imputed. Skipping Task 6.")
    else:
        print(f"Selected proteins: {top_proteins}")

        # Get kinase labels
        site_kinase_map = high_m.groupby('PTM_collapse_key')['Kinase Name'].apply(
            lambda x: ','.join(sorted(set(x)))
        ).to_dict()

        n_prots = len(top_proteins)
        ncols_fig = 3 if n_prots > 3 else n_prots
        nrows_fig = (n_prots + ncols_fig - 1) // ncols_fig
        fig, axes = plt.subplots(nrows_fig, ncols_fig, figsize=(16, 4 * nrows_fig))
        if n_prots == 1:
            axes = np.array([[axes]])
        axes = np.atleast_2d(axes)

        for idx, gene in enumerate(top_proteins):
            r, c = divmod(idx, ncols_fig)
            ax = axes[r][c]
            gene_data = success_high[success_high['Gene'] == gene].drop_duplicates('PTM_collapse_key')
            gene_data = gene_data.sort_values('Position')

            row_labels = []
            mat_data = []
            for _, row in gene_data.iterrows():
                kinases = site_kinase_map.get(row['PTM_collapse_key'], '?')
                # Shorten kinase list
                kin_short = kinases if len(kinases) < 20 else kinases[:17] + '...'
                label = f"{row['AA']}{row['Position']} ({kin_short})"
                row_labels.append(label)
                mat_data.append([int(row[col]) for col in comp_cols])

            mat = pd.DataFrame(mat_data, index=row_labels, columns=comp_labels)

            # Custom colormap: 0=light gray, 1=dark green
            from matplotlib.colors import ListedColormap
            cmap_01 = ListedColormap(['#f0f0f0', '#2d6a4f'])

            sns.heatmap(mat, annot=True, fmt='d', cmap=cmap_01, ax=ax,
                        cbar=False, linewidths=0.5, vmin=0, vmax=1,
                        annot_kws={'size': 7})
            ax.set_title(f'{gene} ({len(mat)} sites)', fontsize=11, fontweight='bold')
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=7, rotation=0)
            ax.set_xticklabels(ax.get_xticklabels(), fontsize=7, rotation=0)

        # Hide empty subplots
        for idx in range(n_prots, nrows_fig * ncols_fig):
            r, c = divmod(idx, ncols_fig)
            axes[r][c].set_visible(False)

        plt.suptitle('Per-Compartment Detection Matrix (High-Mobility Matched Proteins, v2)\n'
                     '0 = not detected (gray), 1 = detected (green)',
                     y=1.02, fontsize=13)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, 'fix4_detection_matrix.png'), dpi=200, bbox_inches='tight')
        plt.close()
        print("Saved: fix4_detection_matrix.png")

else:
    # Use site-level detection (rare case where there's actual variation)
    candidates_sorted = sorted(proteins_with_var, key=lambda g: proteins_2plus.get(g, 0), reverse=True)
    top_proteins = candidates_sorted[:6]
    print(f"Selected proteins (with temporal variation): {top_proteins}")

    site_kinase_map = high_m.groupby('PTM_collapse_key')['Kinase Name'].apply(
        lambda x: ','.join(sorted(set(x)))
    ).to_dict()

    n_prots = len(top_proteins)
    ncols_fig = 3 if n_prots > 3 else n_prots
    nrows_fig = (n_prots + ncols_fig - 1) // ncols_fig
    fig, axes = plt.subplots(nrows_fig, ncols_fig, figsize=(16, 4 * nrows_fig))
    if n_prots == 1:
        axes = np.array([[axes]])
    axes = np.atleast_2d(axes)

    from matplotlib.colors import ListedColormap
    cmap_01 = ListedColormap(['#f0f0f0', '#2d6a4f'])

    for idx, gene in enumerate(top_proteins):
        r, c = divmod(idx, ncols_fig)
        ax = axes[r][c]
        gene_data = success_high[success_high['Gene'] == gene].drop_duplicates('PTM_collapse_key')
        gene_data = gene_data.sort_values('Position')

        row_labels = []
        mat_data = []
        for _, row in gene_data.iterrows():
            kinases = site_kinase_map.get(row['PTM_collapse_key'], '?')
            kin_short = kinases if len(kinases) < 20 else kinases[:17] + '...'
            label = f"{row['AA']}{row['Position']} ({kin_short})"
            row_labels.append(label)
            mat_data.append([int(row[f'detected_{tp}']) for tp in all_tps])

        mat = pd.DataFrame(mat_data, index=row_labels, columns=all_tps)

        sns.heatmap(mat, annot=True, fmt='d', cmap=cmap_01, ax=ax,
                    cbar=False, linewidths=0.5, vmin=0, vmax=1)
        ax.set_title(f'{gene} ({len(mat)} sites)', fontsize=11, fontweight='bold')
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=8, rotation=0)

    for idx in range(n_prots, nrows_fig * ncols_fig):
        r, c = divmod(idx, ncols_fig)
        axes[r][c].set_visible(False)

    plt.suptitle('Site-level Detection Matrix (High-Mobility Matched Proteins, v2)\n'
                 '0 = not detected (gray), 1 = detected (green)',
                 y=1.02, fontsize=13)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'fix4_detection_matrix.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: fix4_detection_matrix.png")

# ============================================================
# Task 7: Per-compartment (LocalizationGroup) kinase analysis
# ============================================================
print("\n" + "=" * 60)
print("TASK 7: Per-compartment kinase family analysis")
print("=" * 60)

site_fam_all = merged.drop_duplicates(subset=['PTM_collapse_key', 'kinase_family'])
print(f"Total (site, family) pairs: {len(site_fam_all)}")

ct_loc = pd.crosstab(site_fam_all['kinase_family'], site_fam_all['LocalizationGroup'])

# Order rows by total count, Other last
row_totals = ct_loc.sum(axis=1).sort_values(ascending=False)
row_order = row_totals.index.tolist()
if 'Other' in row_order:
    row_order.remove('Other')
    row_order.append('Other')
ct_loc = ct_loc.reindex(row_order)

# Normalize by column (per localization group)
ct_loc_prop = ct_loc.div(ct_loc.sum(axis=0), axis=1)

# Sort columns by total count
col_order = ct_loc.sum(axis=0).sort_values(ascending=False).index.tolist()
ct_loc_prop = ct_loc_prop[col_order]

print("\nProportion table:")
print(ct_loc_prop.round(2).to_string())

fig, ax = plt.subplots(figsize=(12, 8))
sns.heatmap(ct_loc_prop, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax,
            linewidths=0.5, vmin=0)
ax.set_title('Kinase Family Proportion by Localization Group\n(All matched sites)', fontsize=13)
ax.set_xlabel('Localization Group', fontsize=12)
ax.set_ylabel('Kinase Family', fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'task7_compartment_heatmap.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: task7_compartment_heatmap.png")

# ============================================================
# Save results
# ============================================================
print("\n" + "=" * 60)
print("SAVING RESULTS")
print("=" * 60)

merged.to_csv(os.path.join(outdir, 'kinase_matched_sites_v2.csv'), index=False)
print(f"Saved: kinase_matched_sites_v2.csv ({merged.shape[0]} rows, {merged['PTM_collapse_key'].nunique()} unique sites)")

# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("FINAL SUMMARY")
print("=" * 60)

n_total_sites = pd.read_csv(r'D:\博士\Protein contour\Phospho\success_df_M1M2_cleaned.csv',
                             usecols=['PTM_collapse_key'])['PTM_collapse_key'].nunique()
n_matched = merged['PTM_collapse_key'].nunique()
match_rate = n_matched / n_total_sites * 100

print(f"Match rate: {n_matched}/{n_total_sites} = {match_rate:.1f}%")
print(f"Other % before fix: {old_other_pct:.1f}%")
print(f"Other % after fix:  {new_other_pct:.1f}%")

high_matched_n = merged[merged['MobilityLevel'] == 'High']['PTM_collapse_key'].nunique()
print(f"High-mobility matched unique sites: {high_matched_n}")

print(f"\nPer-timepoint sample size (high-mobility matched):")
for tp in timepoints:
    n = high_valid[high_valid['movement_best_time'] == tp]['PTM_collapse_key'].nunique()
    print(f"  {tp}: {n}")

print(f"\nTask 5 chi-square p-value: {p_val:.2e}")

# Main findings
print("\n--- Key temporal patterns ---")
for tp in timepoints:
    tp_data = ct_prop.loc[tp] if tp in ct_prop.index else None
    if tp_data is not None:
        top3 = tp_data.sort_values(ascending=False).head(3)
        top_str = ', '.join([f"{f}({v:.0%})" for f, v in top3.items()])
        print(f"  {tp}: {top_str}")

print(f"\nAll figures saved to: {outdir}")
for f in sorted(os.listdir(outdir)):
    fsize = os.path.getsize(os.path.join(outdir, f))
    print(f"  {f} ({fsize/1024:.0f} KB)")

if match_rate < 10:
    print(f"\n WARNING: Match rate only {match_rate:.1f}% — supplement with OmniPath!")
