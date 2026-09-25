"""
Kinase-substrate mapping analysis for EGF-stimulated HeLa phosphoproteomics
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
outdir = r'D:\博士\Protein contour\Phospho\kinase_results'
os.makedirs(outdir, exist_ok=True)

# ============================================================
# Task 1: Load data
# ============================================================
print("=" * 60)
print("TASK 1: Loading data")
print("=" * 60)

# Load success_df
success_df = pd.read_csv(r'D:\博士\Protein contour\Phospho\success_df_M1M2_cleaned.csv')
print(f"success_df shape: {success_df.shape}")
print(f"success_df columns: {list(success_df.columns)}")

# Load high/low mobility
high_mob = pd.read_csv(r'D:\博士\Protein contour\Phospho\high_mobility_sites.csv')
low_mob = pd.read_csv(r'D:\博士\Protein contour\Phospho\low_mobility_sites.csv')
print(f"\nhigh_mobility shape: {high_mob.shape}")
print(f"low_mobility shape: {low_mob.shape}")

# Check if success_df already has MobilityLevel
if 'MobilityLevel' not in success_df.columns:
    # Merge mobility info
    high_mob_slim = high_mob[['PTM_collapse_key', 'MobilityLevel', 'movement_best_time', 'movement_score']].copy()
    low_mob_slim = low_mob[['PTM_collapse_key', 'MobilityLevel', 'movement_best_time', 'movement_score']].copy()
    mob_info = pd.concat([high_mob_slim, low_mob_slim], ignore_index=True)
    success_df = success_df.merge(mob_info, on='PTM_collapse_key', how='left')
    print(f"\nAfter merging mobility info: {success_df.shape}")
else:
    print("\nMobilityLevel already present in success_df")

print(f"\nMobilityLevel distribution:\n{success_df['MobilityLevel'].value_counts(dropna=False)}")
print(f"\nmovement_best_time distribution:\n{success_df['movement_best_time'].value_counts(dropna=False)}")
print(f"\nsuccess_df first 5 rows:")
print(success_df[['PTM_collapse_key', 'Gene', 'AA', 'Position', 'MobilityLevel', 'movement_best_time']].head())

# Load kinase table
kinase_raw = pd.read_csv(r'C:\Users\cc\Zotero\storage\LX2UQSAS\phospho_data\kinace_ksi_source_full.csv')
print(f"\nkinase_raw shape: {kinase_raw.shape}")
print(f"kinase_raw columns: {list(kinase_raw.columns)}")

# Remove empty rows (all NaN or empty strings)
kinase_df = kinase_raw.dropna(subset=['Kinase Name', 'Substrate Name', 'Site'], how='any')
kinase_df = kinase_df[kinase_df['Kinase Name'].str.strip() != '']
kinase_df = kinase_df[kinase_df['Substrate Name'].str.strip() != '']
kinase_df = kinase_df[kinase_df['Site'].str.strip() != '']
print(f"kinase_df after cleaning: {kinase_df.shape}")
print(f"\nkinase_df first 5 rows:")
print(kinase_df[['Kinase Name', 'Substrate Name', 'Site']].head())

# ============================================================
# Task 2: Format alignment + Match
# ============================================================
print("\n" + "=" * 60)
print("TASK 2: Format alignment + Merge")
print("=" * 60)

# Extract residue and position from kinase Site column (e.g., S422 -> S, 422)
kinase_df = kinase_df.copy()
kinase_df['Site_clean'] = kinase_df['Site'].str.strip()
kinase_df['Residue_kinase'] = kinase_df['Site_clean'].str.extract(r'^([STY])', expand=False)
kinase_df['Position_kinase'] = kinase_df['Site_clean'].str.extract(r'(\d+)$', expand=False)
kinase_df['Position_kinase'] = pd.to_numeric(kinase_df['Position_kinase'], errors='coerce')

# Drop rows where extraction failed
kinase_df = kinase_df.dropna(subset=['Residue_kinase', 'Position_kinase'])
kinase_df['Position_kinase'] = kinase_df['Position_kinase'].astype(int)
print(f"kinase_df after extracting residue/position: {kinase_df.shape}")
print(f"Residue distribution in kinase_df: {kinase_df['Residue_kinase'].value_counts().to_dict()}")

# Prepare success_df merge keys
success_df['Position'] = success_df['Position'].astype(int)
# AA column is the residue in success_df
success_df_merge = success_df.rename(columns={'AA': 'Residue'}) if 'Residue' not in success_df.columns else success_df

# Merge on Gene = Substrate Name, Residue = Residue_kinase, Position = Position_kinase
merged = success_df.merge(
    kinase_df[['Kinase', 'Kinase Name', 'Substrate Name', 'Residue_kinase', 'Position_kinase', 'Site', 'Source Database', 'Evidence']],
    left_on=['Gene', 'AA', 'Position'],
    right_on=['Substrate Name', 'Residue_kinase', 'Position_kinase'],
    how='inner'
)

n_matched_sites = merged['PTM_collapse_key'].nunique()
n_total = success_df['PTM_collapse_key'].nunique()
match_rate = n_matched_sites / n_total * 100

print(f"\nTotal unique sites in success_df: {n_total}")
print(f"Matched sites (with kinase info): {n_matched_sites}")
print(f"Match rate: {match_rate:.1f}%")

# High/Low mobility match counts
high_matched = merged[merged['MobilityLevel'] == 'High']['PTM_collapse_key'].nunique()
low_matched = merged[merged['MobilityLevel'] == 'Low']['PTM_collapse_key'].nunique()
print(f"High-mobility matched sites: {high_matched}")
print(f"Low-mobility matched sites: {low_matched}")

if match_rate < 10:
    print("\n⚠️  WARNING: Match rate < 10%! Consider supplementing with OmniPath data.")

# ============================================================
# Task 3: Kinase family classification
# ============================================================
print("\n" + "=" * 60)
print("TASK 3: Kinase family classification")
print("=" * 60)

kinase_family_map = {}
families = {
    'MAPK cascade': ['MAPK1', 'MAPK3', 'MAPK8', 'MAPK14', 'MAP2K1', 'MAP2K2'],
    'CDK': ['CDK1', 'CDK2', 'CDK4', 'CDK5', 'CDK6', 'CDK7', 'CDK9', 'CDK12', 'CDK13'],
    'CK1/CK2': ['CSNK1A1', 'CSNK1D', 'CSNK1E', 'CSNK2A1', 'CSNK2A2'],
    'PKA/PKC': ['PRKACA', 'PRKACB', 'PRKCA', 'PRKCD', 'PRKCI'],
    'AKT/mTOR': ['AKT1', 'AKT2', 'MTOR', 'RPS6KB1'],
    'RTK': ['EGFR', 'ERBB2', 'MET', 'IGF1R', 'SRC', 'FYN'],
    'GSK3': ['GSK3A', 'GSK3B'],
}
for fam, members in families.items():
    for m in members:
        kinase_family_map[m] = fam

merged['kinase_family'] = merged['Kinase Name'].map(kinase_family_map).fillna('Other')

print("\nKinase family site counts (note: one site can have multiple kinases):")
fam_counts = merged.groupby('kinase_family')['PTM_collapse_key'].nunique().sort_values(ascending=False)
print(fam_counts)

# ============================================================
# Task 4: Kinase family by movement_best_time (High mobility only)
# ============================================================
print("\n" + "=" * 60)
print("TASK 4: Kinase family by movement_best_time")
print("=" * 60)

high_matched_df = merged[merged['MobilityLevel'] == 'High'].copy()
print(f"High-mobility matched rows (with possible multiple kinases per site): {len(high_matched_df)}")

# For proportion analysis, use unique (site, kinase_family) pairs to avoid double-counting within same family
high_site_fam = high_matched_df.drop_duplicates(subset=['PTM_collapse_key', 'kinase_family'])

timepoints = ['0min', '2min', '8min', '20min', '90min']
# Filter to valid timepoints
high_site_fam_valid = high_site_fam[high_site_fam['movement_best_time'].isin(timepoints)]

print(f"\nSites per timepoint:")
for tp in timepoints:
    n = high_site_fam_valid[high_site_fam_valid['movement_best_time'] == tp]['PTM_collapse_key'].nunique()
    top_fam = high_site_fam_valid[high_site_fam_valid['movement_best_time'] == tp]['kinase_family'].value_counts()
    top = top_fam.index[0] if len(top_fam) > 0 else 'N/A'
    print(f"  {tp}: {n} sites, top family = {top}")

# Build proportion table
ct = pd.crosstab(high_site_fam_valid['movement_best_time'], high_site_fam_valid['kinase_family'])
ct = ct.reindex(timepoints).fillna(0)
ct_prop = ct.div(ct.sum(axis=1), axis=0)

print("\nProportion table:")
print(ct_prop.round(3))

# --- Stacked bar chart ---
fig, ax = plt.subplots(figsize=(10, 6))
# Use a nice color palette
all_families = sorted(ct_prop.columns.tolist())
# Put "Other" last
if 'Other' in all_families:
    all_families.remove('Other')
    all_families.append('Other')

colors = sns.color_palette("Set2", len(all_families))
color_map = dict(zip(all_families, colors))

bottom = np.zeros(len(ct_prop))
for fam in all_families:
    if fam in ct_prop.columns:
        vals = ct_prop[fam].values
        ax.bar(ct_prop.index, vals, bottom=bottom, label=fam, color=color_map[fam])
        bottom += vals

ax.set_xlabel('Movement Best Time')
ax.set_ylabel('Proportion')
ax.set_title('Kinase Family Composition by Movement Timepoint\n(High-mobility matched sites)')
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'task4_stacked_bar.png'), dpi=200, bbox_inches='tight')
plt.show()
print("Saved: task4_stacked_bar.png")

# --- Heatmap ---
fig, ax = plt.subplots(figsize=(10, 6))
hm_data = ct_prop[all_families].T
hm_data = hm_data.reindex(all_families)
sns.heatmap(hm_data, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax, linewidths=0.5)
ax.set_title('Kinase Family Proportion Heatmap\n(High-mobility matched sites)')
ax.set_xlabel('Movement Best Time')
ax.set_ylabel('Kinase Family')
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'task4_heatmap.png'), dpi=200, bbox_inches='tight')
plt.show()
print("Saved: task4_heatmap.png")

# ============================================================
# Task 5: High vs Low mobility kinase composition
# ============================================================
print("\n" + "=" * 60)
print("TASK 5: High vs Low mobility kinase comparison")
print("=" * 60)

# Use unique (site, kinase_family) pairs
site_fam = merged.drop_duplicates(subset=['PTM_collapse_key', 'kinase_family'])
site_fam_hl = site_fam[site_fam['MobilityLevel'].isin(['High', 'Low'])]

ct_hl = pd.crosstab(site_fam_hl['MobilityLevel'], site_fam_hl['kinase_family'])
ct_hl_prop = ct_hl.div(ct_hl.sum(axis=1), axis=0)

print("Counts:")
print(ct_hl)
print("\nProportions:")
print(ct_hl_prop.round(3))

# Chi-square test
chi2, p_val, dof, expected = chi2_contingency(ct_hl)
print(f"\nChi-square test: chi2={chi2:.2f}, p={p_val:.2e}, dof={dof}")

# Grouped bar chart
fig, ax = plt.subplots(figsize=(12, 6))
ct_hl_prop_plot = ct_hl_prop[all_families] if all(f in ct_hl_prop.columns for f in all_families) else ct_hl_prop
x = np.arange(len(ct_hl_prop_plot.columns))
width = 0.35

bars_high = ax.bar(x - width/2, ct_hl_prop_plot.loc['High'] if 'High' in ct_hl_prop_plot.index else [0]*len(x),
                   width, label='High Mobility', color='#e74c3c', alpha=0.8)
bars_low = ax.bar(x + width/2, ct_hl_prop_plot.loc['Low'] if 'Low' in ct_hl_prop_plot.index else [0]*len(x),
                  width, label='Low Mobility', color='#3498db', alpha=0.8)

ax.set_xlabel('Kinase Family')
ax.set_ylabel('Proportion')
ax.set_title(f'Kinase Family Distribution: High vs Low Mobility\n(Chi-square p = {p_val:.2e})')
ax.set_xticks(x)
ax.set_xticklabels(ct_hl_prop_plot.columns, rotation=45, ha='right')
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'task5_high_vs_low.png'), dpi=200, bbox_inches='tight')
plt.show()
print("Saved: task5_high_vs_low.png")

# ============================================================
# Task 6: 0/1 detection matrix
# ============================================================
print("\n" + "=" * 60)
print("TASK 6: 0/1 detection matrix for top proteins")
print("=" * 60)

# High-mobility matched sites
high_m = merged[merged['MobilityLevel'] == 'High'].copy()

# Count phosphosites per protein (unique sites)
protein_site_count = high_m.groupby('Gene')['PTM_collapse_key'].nunique().sort_values(ascending=False)
proteins_3plus = protein_site_count[protein_site_count >= 3]
print(f"Proteins with >= 3 high-mobility matched sites: {len(proteins_3plus)}")
print(f"Top 10:\n{protein_site_count.head(10)}")

# Select top 5 proteins
top5 = protein_site_count.head(5).index.tolist()

# Timepoint columns in success_df (binary detection)
tp_cols = {
    '0min': ['0min_Cyt', '0min_Mem', '0min_Nuc'],
    '2min': ['2min_Cyt', '2min_Mem', '2min_Nuc'],
    '8min': ['8min_Cyt', '8min_Mem', '8min_Nuc'],
    '20min': ['20min_Cyt', '20min_Mem', '20min_Nuc'],
    '90min': ['90min_Cyt', '90min_Mem', '90min_Nuc'],
}

# For each timepoint, a site is "detected" if any compartment has signal (binary 1)
for tp, cols in tp_cols.items():
    success_df[f'detected_{tp}'] = success_df[cols].max(axis=1).astype(int)

# Build per-protein matrix
fig, axes = plt.subplots(len(top5), 1, figsize=(10, 3 * len(top5)))
if len(top5) == 1:
    axes = [axes]

for idx, gene in enumerate(top5):
    # Get unique sites for this gene from merged data
    gene_sites = high_m[high_m['Gene'] == gene].drop_duplicates('PTM_collapse_key')
    # Get kinase info for labeling
    site_kinase = high_m[high_m['Gene'] == gene].groupby('PTM_collapse_key')['Kinase Name'].apply(
        lambda x: ','.join(sorted(set(x)))
    ).to_dict()

    # Build matrix
    site_keys = sorted(gene_sites['PTM_collapse_key'].unique())
    matrix_data = []
    row_labels = []

    for sk in site_keys:
        row = success_df[success_df['PTM_collapse_key'] == sk]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        kinases = site_kinase.get(sk, '?')
        label = f"{row['AA']}{row['Position']} ({kinases})"
        row_labels.append(label)
        matrix_data.append([row[f'detected_{tp}'] for tp in timepoints])

    mat = pd.DataFrame(matrix_data, index=row_labels, columns=timepoints)

    sns.heatmap(mat, annot=True, fmt='d', cmap='YlGn', ax=axes[idx],
                cbar_kws={'label': 'Detected'}, linewidths=0.5,
                vmin=0, vmax=1)
    axes[idx].set_title(f'{gene} ({len(mat)} phosphosites)', fontsize=12, fontweight='bold')
    axes[idx].set_yticklabels(axes[idx].get_yticklabels(), fontsize=8)

plt.suptitle('Phosphosite Detection Matrix (Top 5 Proteins, High-Mobility Matched)', y=1.01, fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'task6_detection_matrix.png'), dpi=200, bbox_inches='tight')
plt.show()
print("Saved: task6_detection_matrix.png")

# ============================================================
# Task 7: Save results
# ============================================================
print("\n" + "=" * 60)
print("TASK 7: Saving results")
print("=" * 60)

merged.to_csv(os.path.join(outdir, 'kinase_matched_sites.csv'), index=False)
print(f"Saved: kinase_matched_sites.csv ({merged.shape[0]} rows, {merged['PTM_collapse_key'].nunique()} unique sites)")

# Summary
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Total sites in success_df: {n_total}")
print(f"Matched sites with kinase info: {n_matched_sites} ({match_rate:.1f}%)")
print(f"  - High mobility matched: {high_matched}")
print(f"  - Low mobility matched: {low_matched}")
print(f"Kinase families found: {merged['kinase_family'].nunique()}")
print(f"Unique kinases in matched data: {merged['Kinase Name'].nunique()}")
print(f"\nFiles saved to {outdir}:")
for f in sorted(os.listdir(outdir)):
    print(f"  {f}")

if match_rate < 10:
    print(f"\n⚠️  Match rate is only {match_rate:.1f}% — consider supplementing with OmniPath or PhosphoSitePlus data!")
