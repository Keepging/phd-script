"""
Kinase-Substrate Protein-Level Abundance Correlation Analysis
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from matplotlib.colors import ListedColormap
import re
import os
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 150

outdir = r'D:\博士\Protein contour\Phospho\kinase_correlation'
os.makedirs(outdir, exist_ok=True)

timepoints = ['0min', '2min', '8min', '20min', '90min']
compartments = ['Cyt', 'Mem', 'Nuc']
cond_order = [f'{c}_{t}' for c in compartments for t in timepoints]  # Cyt_0min, ..., Nuc_90min

# ============================================================
# STEP 1: Build protein-level abundance matrix from S5A
# ============================================================
print("=" * 60)
print("STEP 1: Build protein abundance matrix from S5A")
print("=" * 60)

s5a = pd.read_excel(r'D:\博士\Protein contour\Phospho\Full_data.xlsx',
                     sheet_name='S5A-Output SN Hela+EGF PROT')
print(f"S5A shape: {s5a.shape}")
print(f"Metadata columns: {list(s5a.columns[:4])}")

# Parse data columns
pattern = r'EGF_(CTRL|\d+min)_(FR\d+)_(Rep\d+)'
data_cols = [c for c in s5a.columns if re.search(pattern, str(c))]
print(f"Data columns: {len(data_cols)}")

fr_to_comp = {'FR1': 'Cyt', 'FR2': 'Cyt', 'FR3': 'Mem', 'FR4': 'Mem', 'FR5': 'Nuc', 'FR6': 'Nuc'}

# Group columns by (compartment, timepoint)
comp_tp_cols = {}
for c in data_cols:
    m = re.search(pattern, str(c))
    tp = m.group(1)
    if tp == 'CTRL':
        tp = '0min'
    fr = m.group(2)
    comp = fr_to_comp[fr]
    key = f'{comp}_{tp}'
    comp_tp_cols.setdefault(key, []).append(c)

# Coerce data columns to numeric (errors='coerce' turns non-numeric into NaN)
for c in data_cols:
    s5a[c] = pd.to_numeric(s5a[c], errors='coerce')

# Compute mean for each (gene, compartment, timepoint)
# For each row compute mean of non-NaN, set NaN if < 2 non-NaN values
abund_data = {}
for cond, cols in comp_tp_cols.items():
    sub = s5a[cols]
    n_non_nan = sub.notna().sum(axis=1)
    mean_vals = sub.mean(axis=1)
    mean_vals[n_non_nan < 2] = np.nan
    abund_data[cond] = mean_vals

abund_raw = pd.DataFrame(abund_data)
abund_raw['PG.Genes'] = s5a['PG.Genes']

# Expand rows with semicolons (one gene per row)
expanded_rows = []
for _, row in abund_raw.iterrows():
    genes = row['PG.Genes']
    if pd.isna(genes):
        continue
    for g in str(genes).split(';'):
        g = g.strip()
        if g == '':
            continue
        new_row = {'Gene': g}
        for cond in cond_order:
            new_row[cond] = row[cond]
        expanded_rows.append(new_row)

abund_expanded = pd.DataFrame(expanded_rows)

# Deduplicate: for each gene, take mean across rows (in case of ambiguous protein groups)
abund_matrix = abund_expanded.groupby('Gene')[cond_order].mean()

print(f"\nProtein abundance matrix shape: {abund_matrix.shape}")
print(f"Genes with at least 1 non-NaN: {(abund_matrix.notna().any(axis=1)).sum()}")
nan_pct = abund_matrix.isna().mean().mean() * 100
print(f"Overall NaN rate: {nan_pct:.1f}%")
print(f"\nFirst 5 rows:")
print(abund_matrix.head())

abund_matrix.to_csv(os.path.join(outdir, 'protein_abundance_matrix.csv'))
print(f"\nSaved: protein_abundance_matrix.csv")

# ============================================================
# STEP 2: Kinase abundance sub-table + heatmap
# ============================================================
print("\n" + "=" * 60)
print("STEP 2: Kinase abundance sub-table")
print("=" * 60)

kinase_df = pd.read_csv(r'C:\Users\cc\Zotero\storage\LX2UQSAS\phospho_data\kinace_ksi_source_full.csv')
kinase_df = kinase_df.dropna(subset=['Kinase Name', 'Substrate Name', 'Site'])
kinase_df = kinase_df[kinase_df['Kinase Name'].str.strip() != '']

all_kinases = kinase_df['Kinase Name'].str.strip().unique()
print(f"Unique kinase names in KinAce: {len(all_kinases)}")

kinases_found = [k for k in all_kinases if k in abund_matrix.index]
print(f"Kinases found in abundance matrix: {len(kinases_found)}/{len(all_kinases)}")

kinase_abund = abund_matrix.loc[kinases_found].copy()
# Keep only kinases with data in at least 2 conditions (any compartment/timepoint)
kinase_abund_hm = kinase_abund[kinase_abund.notna().sum(axis=1) >= 2]
print(f"Kinases with data in >=2 conditions: {len(kinase_abund_hm)}")

# Z-score per row (ignoring NaN)
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

kinase_z = kinase_abund_hm.apply(zscore_row, axis=1)
# Fill NaN with 0 for clustering
kinase_z_filled = kinase_z.fillna(0)

# Main heatmap (all kinases with data)
from scipy.cluster.hierarchy import linkage, leaves_list
if len(kinase_z_filled) > 1:
    link = linkage(kinase_z_filled.values, method='average', metric='euclidean')
    order = leaves_list(link)
    kinase_z_sorted = kinase_z.iloc[order]
else:
    kinase_z_sorted = kinase_z

fig, ax = plt.subplots(figsize=(14, max(8, len(kinase_z_sorted) * 0.25)))
sns.heatmap(kinase_z_sorted, cmap='RdBu_r', center=0, ax=ax,
            cbar_kws={'label': 'Z-score (per kinase)'},
            xticklabels=True, yticklabels=True, linewidths=0.1, linecolor='#eeeeee')
ax.set_title(f'Kinase Protein Abundance (Z-score)\n{len(kinase_z_sorted)} kinases × 15 conditions', fontsize=13)
ax.set_xlabel('Condition (Compartment_Timepoint)', fontsize=11)
ax.set_ylabel('Kinase', fontsize=11)
# Add compartment separators
for i in [5, 10]:
    ax.axvline(x=i, color='black', linewidth=2)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'step2_kinase_abundance_heatmap.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: step2_kinase_abundance_heatmap.png")

# Key kinases subset
key_kinases = ['EGFR', 'MAPK1', 'MAPK3', 'AKT1', 'AKT2', 'MTOR', 'CDK1', 'CDK2', 'CDK9',
               'CSNK2A1', 'GSK3B', 'SRC', 'PRKCA', 'PRKCD', 'RPS6KB1', 'PLK1', 'CHEK1', 'AURKB']
key_found = [k for k in key_kinases if k in kinase_abund.index]
print(f"\nKey EGF-pathway kinases found: {len(key_found)}/{len(key_kinases)}")
missing = [k for k in key_kinases if k not in kinase_abund.index]
if missing:
    print(f"Missing: {missing}")

if key_found:
    key_z = kinase_abund.loc[key_found].apply(zscore_row, axis=1)
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(key_z, cmap='RdBu_r', center=0, ax=ax,
                cbar_kws={'label': 'Z-score'}, annot=True, fmt='.2f',
                linewidths=0.5, linecolor='white', annot_kws={'size': 7})
    ax.set_title('Key EGF-Pathway Kinases — Protein Abundance (Z-score)', fontsize=13)
    ax.set_xlabel('Condition', fontsize=11)
    ax.set_ylabel('Kinase', fontsize=11)
    for i in [5, 10]:
        ax.axvline(x=i, color='black', linewidth=2)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'step2_key_kinases_heatmap.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: step2_key_kinases_heatmap.png")

# ============================================================
# STEP 3: Substrate gene list + kinase-substrate pairs
# ============================================================
print("\n" + "=" * 60)
print("STEP 3: Substrate list from KinAce")
print("=" * 60)

all_substrates = kinase_df['Substrate Name'].str.strip().unique()
print(f"Unique substrate names: {len(all_substrates)}")

ks_pairs = kinase_df[['Kinase Name', 'Substrate Name', 'Site']].drop_duplicates()
ks_pairs['Kinase Name'] = ks_pairs['Kinase Name'].str.strip()
ks_pairs['Substrate Name'] = ks_pairs['Substrate Name'].str.strip()
print(f"Unique kinase-substrate-site rows: {len(ks_pairs)}")

# Unique kinase-substrate pairs (ignoring site)
ks_unique = ks_pairs[['Kinase Name', 'Substrate Name']].drop_duplicates()
print(f"Unique kinase-substrate pairs: {len(ks_unique)}")

# ============================================================
# STEP 4: Substrate abundance sub-table
# ============================================================
print("\n" + "=" * 60)
print("STEP 4: Substrate abundance sub-table")
print("=" * 60)

substrates_found = [s for s in all_substrates if s in abund_matrix.index]
print(f"Substrates found in abundance matrix: {len(substrates_found)}/{len(all_substrates)}")

substrate_abund = abund_matrix.loc[substrates_found].copy()
substrate_abund.to_csv(os.path.join(outdir, 'substrate_abundance_matrix.csv'))
print(f"Saved: substrate_abundance_matrix.csv (shape {substrate_abund.shape})")

# ============================================================
# STEP 5: Kinase-substrate correlation
# ============================================================
print("\n" + "=" * 60)
print("STEP 5: Kinase-substrate correlation")
print("=" * 60)

cyt_cols = [c for c in cond_order if c.startswith('Cyt')]
mem_cols = [c for c in cond_order if c.startswith('Mem')]
nuc_cols = [c for c in cond_order if c.startswith('Nuc')]

def safe_pearson(x, y):
    """Pearson r, p-value. Returns (nan, nan) if insufficient data."""
    mask = (~np.isnan(x)) & (~np.isnan(y))
    if mask.sum() < 3:
        return np.nan, np.nan
    xm = x[mask]
    ym = y[mask]
    if np.std(xm) == 0 or np.std(ym) == 0:
        return np.nan, np.nan
    try:
        r, p = pearsonr(xm, ym)
        return r, p
    except Exception:
        return np.nan, np.nan

# Use kinase-substrate-site rows to retain Site information
ks_rows = ks_pairs[
    ks_pairs['Kinase Name'].isin(kinase_abund.index) &
    ks_pairs['Substrate Name'].isin(substrate_abund.index)
].copy()
print(f"Pairs with both kinase and substrate in matrix: {len(ks_rows)}")

# Precompute kinase/substrate vectors
kinase_vecs = {k: kinase_abund.loc[k, cond_order].values.astype(float)
               for k in kinase_abund.index}
sub_vecs = {s: substrate_abund.loc[s, cond_order].values.astype(float)
            for s in substrate_abund.index}

results = []
# Deduplicate so we compute correlation once per (kinase, substrate) pair
unique_pairs = ks_rows[['Kinase Name', 'Substrate Name']].drop_duplicates()
# Map pair -> list of sites
pair_sites = ks_rows.groupby(['Kinase Name', 'Substrate Name'])['Site'].apply(
    lambda x: ','.join(sorted(set(x)))
).to_dict()

cyt_idx = [cond_order.index(c) for c in cyt_cols]
mem_idx = [cond_order.index(c) for c in mem_cols]
nuc_idx = [cond_order.index(c) for c in nuc_cols]

for _, row in unique_pairs.iterrows():
    k = row['Kinase Name']
    s = row['Substrate Name']
    kv = kinase_vecs[k]
    sv = sub_vecs[s]

    r_all, p_all = safe_pearson(kv, sv)
    r_cyt, p_cyt = safe_pearson(kv[cyt_idx], sv[cyt_idx])
    r_mem, p_mem = safe_pearson(kv[mem_idx], sv[mem_idx])
    r_nuc, p_nuc = safe_pearson(kv[nuc_idx], sv[nuc_idx])

    sites = pair_sites.get((k, s), '')

    results.append({
        'Kinase': k, 'Substrate': s, 'Sites': sites,
        'Overall_Pearson_r': r_all, 'Overall_Pearson_p': p_all,
        'Cyt_r': r_cyt, 'Cyt_p': p_cyt,
        'Mem_r': r_mem, 'Mem_p': p_mem,
        'Nuc_r': r_nuc, 'Nuc_p': p_nuc,
    })

corr_df = pd.DataFrame(results)
print(f"\nCorrelations computed: {len(corr_df)}")

# Filter for valid correlations
valid = corr_df.dropna(subset=['Overall_Pearson_r'])
print(f"Valid (non-NaN) correlations: {len(valid)}")

sig = valid[valid['Overall_Pearson_p'] < 0.05]
print(f"Significant (p<0.05): {len(sig)}")
print(f"  Positive: {(sig['Overall_Pearson_r'] > 0).sum()}")
print(f"  Negative: {(sig['Overall_Pearson_r'] < 0).sum()}")
print(f"Mean r: {valid['Overall_Pearson_r'].mean():.3f}")
print(f"Median r: {valid['Overall_Pearson_r'].median():.3f}")

print("\nTop 10 positively correlated pairs:")
top_pos = valid.nlargest(10, 'Overall_Pearson_r')[['Kinase', 'Substrate', 'Overall_Pearson_r', 'Overall_Pearson_p']]
print(top_pos.to_string(index=False))

print("\nTop 10 negatively correlated pairs:")
top_neg = valid.nsmallest(10, 'Overall_Pearson_r')[['Kinase', 'Substrate', 'Overall_Pearson_r', 'Overall_Pearson_p']]
print(top_neg.to_string(index=False))

corr_df.to_csv(os.path.join(outdir, 'kinase_substrate_correlation.csv'), index=False)
print(f"\nSaved: kinase_substrate_correlation.csv")

# --- Plot A: distribution ---
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(valid['Overall_Pearson_r'], bins=60, color='#4a90d9', edgecolor='white', alpha=0.8)
mean_r = valid['Overall_Pearson_r'].mean()
median_r = valid['Overall_Pearson_r'].median()
ax.axvline(mean_r, color='red', linestyle='--', linewidth=2, label=f'Mean = {mean_r:.3f}')
ax.axvline(median_r, color='green', linestyle='--', linewidth=2, label=f'Median = {median_r:.3f}')
ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
ax.set_xlabel('Overall Pearson r (kinase vs substrate abundance)', fontsize=12)
ax.set_ylabel('Number of pairs', fontsize=12)
ax.set_title(f'Distribution of Kinase-Substrate Protein Abundance Correlations\n(n={len(valid)} valid pairs, {len(sig)} significant at p<0.05)', fontsize=12)
ax.legend(fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'step5_correlation_distribution.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: step5_correlation_distribution.png")

# --- Plot B: top pairs heatmap ---
# Sort by p-value (ascending) and take top 20
top20 = valid.nsmallest(20, 'Overall_Pearson_p').reset_index(drop=True)
# Build a heatmap: rows = pair label, cols = 15 conditions, side-by-side kinase/substrate
# We stack: for each pair, two rows (kinase z, substrate z)
hm_rows = []
hm_labels = []
for _, row in top20.iterrows():
    k = row['Kinase']
    s = row['Substrate']
    kv = kinase_abund.loc[k, cond_order].values.astype(float)
    sv = substrate_abund.loc[s, cond_order].values.astype(float)
    # z-score
    def z(v):
        mask = ~np.isnan(v)
        if mask.sum() < 2 or np.nanstd(v) == 0:
            return np.zeros_like(v)
        return (v - np.nanmean(v)) / np.nanstd(v)
    kz = z(kv)
    sz = z(sv)
    hm_rows.append(kz)
    hm_labels.append(f"[K] {k} → {s} (r={row['Overall_Pearson_r']:.2f})")
    hm_rows.append(sz)
    hm_labels.append(f"[S] {k} → {s}")

hm_mat = pd.DataFrame(hm_rows, index=hm_labels, columns=cond_order)
fig, ax = plt.subplots(figsize=(16, 12))
sns.heatmap(hm_mat, cmap='RdBu_r', center=0, ax=ax,
            cbar_kws={'label': 'Z-score'}, linewidths=0.3, linecolor='#eeeeee')
ax.set_title('Top 20 Most Significant Kinase-Substrate Pairs\n(Kinase [K] and Substrate [S] z-scores paired)', fontsize=13)
ax.set_xlabel('Condition', fontsize=11)
for i in [5, 10]:
    ax.axvline(x=i, color='black', linewidth=2)
# Thick separators between pairs
for i in range(2, len(hm_labels), 2):
    ax.axhline(y=i, color='black', linewidth=1)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=8, rotation=0)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'step5_top_pairs_heatmap.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: step5_top_pairs_heatmap.png")

# --- Plot C: scatter examples ---
top3_pos = valid.nlargest(3, 'Overall_Pearson_r')
top3_neg = valid.nsmallest(3, 'Overall_Pearson_r')
examples = pd.concat([top3_pos, top3_neg]).reset_index(drop=True)

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
comp_colors = {'Cyt': '#1f77b4', 'Mem': '#ff7f0e', 'Nuc': '#2ca02c'}

for idx, (_, row) in enumerate(examples.iterrows()):
    ax = axes[idx // 3][idx % 3]
    k = row['Kinase']
    s = row['Substrate']
    kv = kinase_abund.loc[k, cond_order].values.astype(float)
    sv = substrate_abund.loc[s, cond_order].values.astype(float)

    for i, cond in enumerate(cond_order):
        comp, tp = cond.split('_')
        if np.isnan(kv[i]) or np.isnan(sv[i]):
            continue
        ax.scatter(kv[i], sv[i], color=comp_colors[comp], s=120, alpha=0.8,
                   edgecolors='black', linewidth=0.5)
        ax.annotate(tp.replace('min',''), (kv[i], sv[i]), fontsize=8,
                    xytext=(5, 5), textcoords='offset points')

    # Legend for compartments
    for comp, c in comp_colors.items():
        ax.scatter([], [], color=c, s=80, label=comp)
    ax.legend(loc='best', fontsize=8)

    sign = '(+)' if row['Overall_Pearson_r'] > 0 else '(-)'
    ax.set_title(f"{sign} {k} → {s}\nr={row['Overall_Pearson_r']:.3f}, p={row['Overall_Pearson_p']:.1e}",
                 fontsize=11)
    ax.set_xlabel(f'{k} abundance', fontsize=10)
    ax.set_ylabel(f'{s} abundance', fontsize=10)
    ax.grid(True, alpha=0.3)

plt.suptitle('Top Kinase-Substrate Correlation Examples (Top 3 +/-)', fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'step5_scatter_examples.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: step5_scatter_examples.png")

# ============================================================
# STEP 6: Cross-reference kinase presence with site detection
# ============================================================
print("\n" + "=" * 60)
print("STEP 6: Cross-reference kinase presence vs site detection")
print("=" * 60)

km = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_results\v2\kinase_matched_sites_v2.csv')
s5c_bin = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_results\v2\s5c_binary_detection.csv')

case_proteins = ['MARCKS', 'RB1', 'NOLC1', 'GYS1', 'LARP1', 'SUDS3']

# S5C binary uses column names like "0min_Cyt"; reformat to match cond_order
# cond_order is Cyt_0min ... we'll use (comp, tp) keys
s5c_cols_map = {f'{tp}_{comp}': f'{comp}_{tp}' for tp in timepoints for comp in compartments}

# Build kinase presence: kinase is "present" in condition if abundance is non-NaN
def kinase_present(kinase_name):
    if kinase_name not in abund_matrix.index:
        return np.full(15, np.nan)
    vals = abund_matrix.loc[kinase_name, cond_order].values.astype(float)
    return (~np.isnan(vals)).astype(int)

fig, axes = plt.subplots(2, 3, figsize=(18, 14))

cmap_01 = ListedColormap(['#ffffff', '#2d6a4f'])

interesting_count = 0  # kinase present + site not detected
surprising_count = 0   # kinase absent + site detected

protein_summaries = {}

for idx, gene in enumerate(case_proteins):
    ax = axes[idx // 3][idx % 3]

    gene_km = km[(km['Gene'] == gene) & (km['MobilityLevel'] == 'High')]
    if len(gene_km) == 0:
        ax.text(0.5, 0.5, f'{gene}: no high-mob matches', ha='center', va='center')
        ax.axis('off')
        continue

    # Get sites grouped by (kinase, site) pairs
    rows_data = []
    row_labels = []
    cell_categories = []

    # For each unique (site, kinase) pair
    site_kinase_pairs = gene_km[['PTM_collapse_key', 'AA', 'Position', 'Kinase Name']].drop_duplicates()
    site_kinase_pairs = site_kinase_pairs.sort_values(['Position', 'Kinase Name'])

    prot_interesting = 0
    prot_surprising = 0

    for _, sk in site_kinase_pairs.iterrows():
        ptm = sk['PTM_collapse_key']
        kinase = sk['Kinase Name']
        aa = sk['AA']
        pos = int(sk['Position'])

        # Kinase presence vector (15 values, 1=present 0=absent)
        k_presence = kinase_present(kinase)

        # Site detection vector from s5c_bin
        site_row = s5c_bin[s5c_bin['PTM_collapse_key'] == ptm]
        if len(site_row) == 0:
            site_det = np.zeros(15, dtype=int)
        else:
            site_det = np.array([int(site_row.iloc[0][f'{tp}_{comp}'])
                                 for comp in compartments for tp in timepoints])

        # Two rows per (site, kinase): kinase presence, site detection
        rows_data.append(k_presence)
        row_labels.append(f"[K] {kinase}")
        rows_data.append(site_det)
        row_labels.append(f"[S] {aa}{pos}")

        # Count categories (treat kinase NaN as absent=0 for counting)
        kp_clean = np.where(np.isnan(k_presence), 0, k_presence).astype(int)
        for i in range(15):
            kp = kp_clean[i]
            sd = site_det[i]
            if kp == 1 and sd == 0:
                prot_interesting += 1
                interesting_count += 1
            elif kp == 0 and sd == 1:
                prot_surprising += 1
                surprising_count += 1

    protein_summaries[gene] = {
        'n_site_kinase_pairs': len(site_kinase_pairs),
        'interesting': prot_interesting,
        'surprising': prot_surprising,
    }

    if not rows_data:
        ax.axis('off')
        continue

    mat = pd.DataFrame(rows_data, index=row_labels, columns=cond_order)
    # Fill NaN with -1 for visualization (gray)
    mat_plot = mat.fillna(0)

    sns.heatmap(mat_plot, annot=True, fmt='.0f', cmap=cmap_01, ax=ax,
                cbar=False, linewidths=0.5, linecolor='#cccccc',
                vmin=0, vmax=1, annot_kws={'size': 7})
    ax.set_title(f'{gene} ({len(site_kinase_pairs)} site-kinase pairs)\n'
                 f'K=present kinase, S=site detected | interesting={prot_interesting}',
                 fontsize=10, fontweight='bold')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7, rotation=0)
    ax.set_xticklabels([c.replace('min','').replace('Cyt','C').replace('Mem','M').replace('Nuc','N').replace('_','-')
                         for c in cond_order], fontsize=7, rotation=45, ha='right')
    for i in [5, 10]:
        ax.axvline(x=i, color='black', linewidth=1.5)
    # Row pair separators
    for i in range(2, len(row_labels), 2):
        ax.axhline(y=i, color='black', linewidth=0.8)

plt.suptitle('Kinase Protein Presence vs Phosphosite Detection\n'
             'Top row [K] = kinase present in condition | Bottom row [S] = site detected in condition',
             fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'step6_kinase_site_crossref.png'), dpi=200, bbox_inches='tight')
plt.close()
print("Saved: step6_kinase_site_crossref.png")

print(f"\nPer-protein counts (interesting = kinase present + site not detected):")
for g, s in protein_summaries.items():
    print(f"  {g}: {s['n_site_kinase_pairs']} (site,kinase) pairs, "
          f"interesting={s['interesting']}, surprising={s['surprising']}")

print(f"\nTotal interesting cells (kinase present, site not detected): {interesting_count}")
print(f"Total surprising cells (kinase absent, site detected): {surprising_count}")

# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("FINAL SUMMARY")
print("=" * 60)
print(f"Step 1: Protein abundance matrix = {abund_matrix.shape} ({nan_pct:.1f}% NaN)")
print(f"Step 2: Kinases in matrix = {len(kinases_found)}/{len(all_kinases)} "
      f"({len(kinase_abund_hm)} with data in >=2 conditions)")
print(f"Step 3: Unique kinase-substrate-site rows = {len(ks_pairs)}, "
      f"unique pairs = {len(ks_unique)}")
print(f"Step 4: Substrates in matrix = {len(substrates_found)}/{len(all_substrates)}")
print(f"Step 5: Valid correlations = {len(valid)}, significant = {len(sig)} "
      f"({(sig['Overall_Pearson_r']>0).sum()} pos, {(sig['Overall_Pearson_r']<0).sum()} neg)")
print(f"Step 6: Interesting cases (kinase present, site not detected) = {interesting_count}")
print(f"         Surprising cases (kinase absent, site detected) = {surprising_count}")

print(f"\nAll files saved to: {outdir}")
for f in sorted(os.listdir(outdir)):
    fsize = os.path.getsize(os.path.join(outdir, f))
    print(f"  {f} ({fsize/1024:.0f} KB)")
