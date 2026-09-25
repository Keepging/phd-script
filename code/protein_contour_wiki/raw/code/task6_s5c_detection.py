"""
Task 6 v3: Detection matrix from S5C raw data vs imputed success_df
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
import re
import os
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 150

outdir = r'D:\博士\Protein contour\Phospho\kinase_results\v2'
os.makedirs(outdir, exist_ok=True)

# ============================================================
# Step 1: Load S5C, parse columns
# ============================================================
print("=" * 60)
print("STEP 1: Load S5C raw data")
print("=" * 60)

s5c = pd.read_excel(r'D:\博士\Protein contour\Phospho\Full_data.xlsx',
                     sheet_name='S5C-Collapse PHOS Hela+EGF')
print(f"S5C shape: {s5c.shape}")
print(f"NaN rate: {s5c.iloc[:, 1:].isna().mean().mean()*100:.1f}%")

# Parse column names
pattern = r'EGF_(CTRL|\d+min)_(FR\d+)_(Rep\d+)'
col_meta = {}  # col_name -> (timepoint, fraction, replicate)
for c in s5c.columns[1:]:
    m = re.search(pattern, str(c))
    if m:
        tp = m.group(1)
        if tp == 'CTRL':
            tp = '0min'
        col_meta[c] = (tp, m.group(2), m.group(3))

print(f"Parsed {len(col_meta)} data columns")

# Map fractions to compartments
fr_to_comp = {
    'FR1': 'Cyt', 'FR2': 'Cyt',
    'FR3': 'Mem', 'FR4': 'Mem',
    'FR5': 'Nuc', 'FR6': 'Nuc',
}

# Group columns by (timepoint, compartment)
tp_comp_cols = {}  # (tp, comp) -> list of column names
for col, (tp, fr, rep) in col_meta.items():
    comp = fr_to_comp[fr]
    key = (tp, comp)
    tp_comp_cols.setdefault(key, []).append(col)

timepoints = ['0min', '2min', '8min', '20min', '90min']
compartments = ['Cyt', 'Mem', 'Nuc']

print("\nColumns per (timepoint, compartment):")
for tp in timepoints:
    for comp in compartments:
        key = (tp, comp)
        n = len(tp_comp_cols.get(key, []))
        print(f"  {tp:5s} {comp}: {n} columns (2 fractions x 4 reps = 8 expected)")

# ============================================================
# Step 2: Build binary detection matrix
# ============================================================
print("\n" + "=" * 60)
print("STEP 2: Build binary detection matrix")
print("=" * 60)

# For each site, each (tp, comp): count non-NaN values. If >= 2 -> detected (1), else 0
binary_results = {}
for tp in timepoints:
    for comp in compartments:
        col_name = f'{tp}_{comp}'
        cols = tp_comp_cols.get((tp, comp), [])
        if cols:
            non_nan_count = s5c[cols].notna().sum(axis=1)
            binary_results[col_name] = (non_nan_count >= 2).astype(int)
        else:
            binary_results[col_name] = pd.Series(0, index=s5c.index)

binary_df = pd.DataFrame(binary_results)
binary_df.insert(0, 'PTM_collapse_key', s5c['PTM_collapse_key'])

print(f"Binary detection matrix shape: {binary_df.shape}")
det_rate = binary_df.iloc[:, 1:].mean().mean()
print(f"Overall detection rate (fraction of 1s): {det_rate*100:.1f}%")
print(f"Per condition detection rate:")
for col in binary_df.columns[1:]:
    rate = binary_df[col].mean()
    print(f"  {col:10s}: {rate*100:.1f}%")

# Compare with success_df binary columns
success_df = pd.read_csv(r'D:\博士\Protein contour\Phospho\success_df_M1M2_cleaned.csv')
# success_df has 0min_Cyt, 0min_Mem, etc. as binary columns
compare_cols = [f'{tp}_{comp}' for tp in timepoints for comp in compartments]

# Merge on PTM_collapse_key
compare = binary_df.merge(success_df[['PTM_collapse_key'] + compare_cols],
                          on='PTM_collapse_key', suffixes=('_s5c', '_imp'))

n_cells = 0
n_diff = 0
for col in compare_cols:
    s5c_col = f'{col}_s5c'
    imp_col = f'{col}_imp'
    if s5c_col in compare.columns and imp_col in compare.columns:
        diff = (compare[s5c_col] != compare[imp_col]).sum()
        n_diff += diff
        n_cells += len(compare)

print(f"\nS5C vs success_df comparison (sites in both: {len(compare)}):")
print(f"  Total cells compared: {n_cells}")
print(f"  Discrepant cells: {n_diff}")
print(f"  Discrepancy rate: {n_diff/n_cells*100:.2f}%")

# Direction of discrepancies
n_imp_extra = 0  # imputed=1 but s5c=0
n_s5c_extra = 0  # s5c=1 but imputed=0
for col in compare_cols:
    s5c_col = f'{col}_s5c'
    imp_col = f'{col}_imp'
    if s5c_col in compare.columns and imp_col in compare.columns:
        n_imp_extra += ((compare[imp_col] == 1) & (compare[s5c_col] == 0)).sum()
        n_s5c_extra += ((compare[s5c_col] == 1) & (compare[imp_col] == 0)).sum()

print(f"  Imputed=1 but S5C=0 (false positives from imputation): {n_imp_extra}")
print(f"  S5C=1 but Imputed=0 (false negatives in imputed):      {n_s5c_extra}")

# ============================================================
# Step 3: Merge with kinase match results
# ============================================================
print("\n" + "=" * 60)
print("STEP 3: Merge with kinase match data")
print("=" * 60)

kinase_matched = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_results\v2\kinase_matched_sites_v2.csv')
high_matched = kinase_matched[kinase_matched['MobilityLevel'] == 'High'].copy()
print(f"High-mobility matched rows: {len(high_matched)}")
print(f"High-mobility matched unique sites: {high_matched['PTM_collapse_key'].nunique()}")

# Get kinase labels per site
site_kinase = high_matched.groupby('PTM_collapse_key')['Kinase Name'].apply(
    lambda x: ','.join(sorted(set(x)))
).to_dict()

# Get site info (AA, Position, Gene)
site_info = high_matched.drop_duplicates('PTM_collapse_key')[['PTM_collapse_key', 'Gene', 'AA', 'Position']].copy()

# Merge binary detection with site info
high_sites = site_info['PTM_collapse_key'].unique()
binary_high = binary_df[binary_df['PTM_collapse_key'].isin(high_sites)].copy()
binary_high = binary_high.merge(site_info, on='PTM_collapse_key', how='left')

print(f"High-mobility matched sites with S5C data: {len(binary_high)}")

# ============================================================
# Step 4: Select proteins and draw S5C detection matrix
# ============================================================
print("\n" + "=" * 60)
print("STEP 4: Select proteins and draw S5C detection matrix")
print("=" * 60)

# Count sites per protein
protein_counts = binary_high.groupby('Gene')['PTM_collapse_key'].nunique().sort_values(ascending=False)
proteins_2plus = protein_counts[protein_counts >= 2]
print(f"Proteins with >= 2 matched high-mob sites: {len(proteins_2plus)}")
print(proteins_2plus)

# Check temporal variation
det_cols_ordered = [f'{tp}_{comp}' for tp in timepoints for comp in compartments]
protein_variation = {}
for gene in proteins_2plus.index:
    gene_data = binary_high[binary_high['Gene'] == gene]
    mat = gene_data[det_cols_ordered].values
    has_zero = (mat == 0).any()
    has_one = (mat == 1).any()
    n_zeros = (mat == 0).sum()
    protein_variation[gene] = {'has_variation': has_zero and has_one, 'n_zeros': n_zeros}

proteins_with_var = [g for g, v in protein_variation.items() if v['has_variation']]
print(f"\nProteins with temporal/compartment variation: {len(proteins_with_var)}")
for g in proteins_with_var:
    print(f"  {g}: {protein_variation[g]['n_zeros']} zero cells")

# Priority list
priority = ['MARCKS', 'RB1', 'SRRM2']
selected = []
for g in priority:
    if g in proteins_with_var:
        selected.append(g)

# Fill remaining from variation list sorted by n_zeros (more variation = more interesting)
remaining = [g for g in proteins_with_var if g not in selected]
remaining.sort(key=lambda g: protein_variation[g]['n_zeros'], reverse=True)
for g in remaining:
    if len(selected) >= 6:
        break
    selected.append(g)

# If still < 6, add from non-variation proteins
if len(selected) < 6:
    for g in proteins_2plus.index:
        if g not in selected and len(selected) < 6:
            selected.append(g)

print(f"\nSelected proteins: {selected}")

# --- Draw S5C detection matrix ---
x_labels = [f'{comp[0]}-{tp.replace("min","")}'
            for tp in timepoints for comp in compartments]

cmap_01 = ListedColormap(['#ffffff', '#2d6a4f'])

n_prots = len(selected)
nrows_fig = 2
ncols_fig = 3
fig, axes = plt.subplots(nrows_fig, ncols_fig, figsize=(18, 14))

for idx, gene in enumerate(selected):
    r, c = divmod(idx, ncols_fig)
    ax = axes[r][c]

    gene_data = binary_high[binary_high['Gene'] == gene].sort_values('Position')

    row_labels = []
    mat_data = []
    for _, row in gene_data.iterrows():
        kinases = site_kinase.get(row['PTM_collapse_key'], '?')
        kin_short = kinases if len(kinases) < 25 else kinases[:22] + '...'
        label = f"{row['AA']}{int(row['Position'])} ({kin_short})"
        row_labels.append(label)
        mat_data.append([int(row[col]) for col in det_cols_ordered])

    mat = pd.DataFrame(mat_data, index=row_labels, columns=x_labels)

    sns.heatmap(mat, annot=True, fmt='d', cmap=cmap_01, ax=ax,
                cbar=False, linewidths=0.8, linecolor='#cccccc',
                vmin=0, vmax=1, annot_kws={'size': 8, 'color': '#333333'})
    ax.set_title(f'{gene} ({len(mat)} sites)', fontsize=12, fontweight='bold')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8, rotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=7, rotation=45, ha='right')

    # Add compartment group separators
    for i in range(1, 5):
        ax.axvline(x=i*3, color='black', linewidth=1.5)

# Hide empty subplots
for idx in range(n_prots, nrows_fig * ncols_fig):
    r, c = divmod(idx, ncols_fig)
    axes[r][c].set_visible(False)

fig.suptitle('Phosphosite Detection Matrix — S5C Raw Data\n'
             '(High-mobility, kinase-matched | white=not detected, green=detected)\n'
             'Detection threshold: >=2 non-NaN values per compartment×timepoint',
             fontsize=13, y=1.03)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'task6_s5c_detection_matrix.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: task6_s5c_detection_matrix.png")

# ============================================================
# Step 5: Draw imputed (success_df) version for comparison
# ============================================================
print("\n" + "=" * 60)
print("STEP 5: Draw imputed (success_df) version for comparison")
print("=" * 60)

# Get success_df binary data for same proteins
imp_cols = [f'{tp}_{comp}' for tp in timepoints for comp in compartments]
success_high = success_df[success_df['PTM_collapse_key'].isin(high_sites)].copy()

fig2, axes2 = plt.subplots(nrows_fig, ncols_fig, figsize=(18, 14))

for idx, gene in enumerate(selected):
    r, c = divmod(idx, ncols_fig)
    ax = axes2[r][c]

    gene_s5c = binary_high[binary_high['Gene'] == gene].sort_values('Position')
    ptm_keys = gene_s5c['PTM_collapse_key'].tolist()
    gene_imp = success_df[success_df['PTM_collapse_key'].isin(ptm_keys)].copy()
    # Reorder to match s5c
    gene_imp = gene_imp.set_index('PTM_collapse_key').reindex(ptm_keys).reset_index()

    row_labels = []
    mat_data = []
    for _, row in gene_imp.iterrows():
        kinases = site_kinase.get(row['PTM_collapse_key'], '?')
        kin_short = kinases if len(kinases) < 25 else kinases[:22] + '...'
        label = f"{row['AA']}{int(row['Position'])} ({kin_short})"
        row_labels.append(label)
        mat_data.append([int(row[col]) for col in imp_cols])

    mat = pd.DataFrame(mat_data, index=row_labels, columns=x_labels)

    sns.heatmap(mat, annot=True, fmt='d', cmap=cmap_01, ax=ax,
                cbar=False, linewidths=0.8, linecolor='#cccccc',
                vmin=0, vmax=1, annot_kws={'size': 8, 'color': '#333333'})
    ax.set_title(f'{gene} ({len(mat)} sites)', fontsize=12, fontweight='bold')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8, rotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=7, rotation=45, ha='right')

    for i in range(1, 5):
        ax.axvline(x=i*3, color='black', linewidth=1.5)

for idx in range(n_prots, nrows_fig * ncols_fig):
    r, c = divmod(idx, ncols_fig)
    axes2[r][c].set_visible(False)

fig2.suptitle('Phosphosite Detection Matrix — Imputed (success_df)\n'
              '(Same proteins as S5C version for comparison)\n'
              'Note: imputed data fills in missing values → more 1s',
              fontsize=13, y=1.03)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'task6_imputed_detection_matrix.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: task6_imputed_detection_matrix.png")

# ============================================================
# Step 6: Save results
# ============================================================
print("\n" + "=" * 60)
print("STEP 6: Save results and summary")
print("=" * 60)

binary_df.to_csv(os.path.join(outdir, 's5c_binary_detection.csv'), index=False)
print(f"Saved: s5c_binary_detection.csv ({binary_df.shape})")

# Per-protein discrepancy for selected proteins
print("\n--- Per-protein S5C vs Imputed discrepancy ---")
for gene in selected:
    gene_s5c = binary_high[binary_high['Gene'] == gene].sort_values('Position')
    ptm_keys = gene_s5c['PTM_collapse_key'].tolist()
    gene_imp = success_df[success_df['PTM_collapse_key'].isin(ptm_keys)]

    n_cells_gene = 0
    n_diff_gene = 0
    for ptm_key in ptm_keys:
        s5c_row = gene_s5c[gene_s5c['PTM_collapse_key'] == ptm_key]
        imp_row = gene_imp[gene_imp['PTM_collapse_key'] == ptm_key]
        if len(s5c_row) == 0 or len(imp_row) == 0:
            continue
        for col in det_cols_ordered:
            v_s5c = int(s5c_row[col].iloc[0])
            v_imp = int(imp_row[col].iloc[0])
            n_cells_gene += 1
            if v_s5c != v_imp:
                n_diff_gene += 1

    rate = n_diff_gene / n_cells_gene * 100 if n_cells_gene > 0 else 0
    n_sites = len(ptm_keys)
    print(f"  {gene:10s}: {n_sites} sites, {n_diff_gene}/{n_cells_gene} cells differ ({rate:.1f}%)")

# Final summary
print("\n" + "=" * 60)
print("FINAL SUMMARY")
print("=" * 60)
print(f"S5C total sites: {len(binary_df)}")
print(f"S5C overall detection rate: {det_rate*100:.1f}%")
print(f"S5C vs imputed total discrepancy: {n_diff}/{n_cells} = {n_diff/n_cells*100:.2f}%")
print(f"  Imputed=1, S5C=0 (imputation artifacts): {n_imp_extra}")
print(f"  S5C=1, Imputed=0 (lost in processing):   {n_s5c_extra}")
print(f"Selected proteins: {selected}")
print(f"\nFiles saved:")
for f in sorted(os.listdir(outdir)):
    if 'task6' in f or 's5c' in f:
        fsize = os.path.getsize(os.path.join(outdir, f))
        print(f"  {f} ({fsize/1024:.0f} KB)")
