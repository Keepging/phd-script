"""
Cell-cycle ANOVA on PXD011836 phosphosites (PerseusOutput),
then cross-reference with CDK substrates from kinase_substrate_correlation.csv.
"""
import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests
import os

OUTDIR = r'D:\博士\Protein contour\Phospho\cellcycle_analysis'
os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# Load PerseusOutput — skip 2 Perseus annotation rows
# ============================================================
print('Loading PerseusOutput...')
df = pd.read_csv(
    r'D:\博士\Protein contour\Phospho_HeLa_CellCycle\PerseusOutput Phospho(STY).txt',
    sep='\t', skiprows=[1, 2])
print(f'Raw rows: {len(df)}')

# Convert quantitative columns to numeric
phase_cols = ['G1_1', 'G1_2', 'G1_3', 'S_1', 'S_2', 'S_3',
              'G2_1', 'G2_2', 'G2_3']
for c in phase_cols:
    df[c] = pd.to_numeric(df[c], errors='coerce')

df['Localization prob'] = pd.to_numeric(df['Localization prob'], errors='coerce')

# ============================================================
# Step 1: Filter Class I (loc prob >= 0.75)
# ============================================================
class1 = df[df['Localization prob'] >= 0.75].copy()
print(f'Class I sites (loc prob >= 0.75): {len(class1)}')

# ============================================================
# Step 2: Drop rows with any NaN in the 9 phase columns
# ============================================================
complete = class1.dropna(subset=phase_cols).copy()
print(f'Complete cases (no NaN in 9 phase cols): {len(complete)}')

# ============================================================
# Step 3: One-way ANOVA per site (G1 vs S vs G2)
# ============================================================
print('Running ANOVA...')
f_stats = []
p_values = []
g1_means = []
s_means = []
g2_means = []

for _, row in complete.iterrows():
    g1 = [row['G1_1'], row['G1_2'], row['G1_3']]
    s = [row['S_1'], row['S_2'], row['S_3']]
    g2 = [row['G2_1'], row['G2_2'], row['G2_3']]
    f_stat, p_val = f_oneway(g1, s, g2)
    f_stats.append(f_stat)
    p_values.append(p_val)
    g1_means.append(np.mean(g1))
    s_means.append(np.mean(s))
    g2_means.append(np.mean(g2))

complete = complete.copy()
complete['F_stat'] = f_stats
complete['p_value'] = p_values
complete['G1_mean'] = g1_means
complete['S_mean'] = s_means
complete['G2_mean'] = g2_means

# ============================================================
# Step 4: BH FDR correction
# ============================================================
print('BH FDR correction...')
reject, fdr, _, _ = multipletests(complete['p_value'].values, method='fdr_bh')
complete['FDR'] = fdr
complete['is_CC_responsive'] = reject

# ============================================================
# Step 5: Determine max phase
# ============================================================
def max_phase(row):
    vals = {'G1': row['G1_mean'], 'S': row['S_mean'], 'G2': row['G2_mean']}
    return max(vals, key=vals.get)

complete['max_phase'] = complete.apply(max_phase, axis=1)

# Parse gene (take first gene if semicolon-separated)
complete['Gene'] = (complete['Gene names'].fillna('').astype(str)
                    .apply(lambda x: x.split(';')[0].strip()))

# Build output
out = complete[['Gene', 'Position', 'Sequence window',
                'F_stat', 'p_value', 'FDR', 'is_CC_responsive',
                'G1_mean', 'S_mean', 'G2_mean', 'max_phase']].copy()
out.columns = ['Gene', 'Position', 'Sequence_window',
               'F_stat', 'p_value', 'FDR', 'is_CC_responsive',
               'G1_mean', 'S_mean', 'G2_mean', 'max_phase']

# Round numeric columns
for c in ['F_stat', 'p_value', 'FDR', 'G1_mean', 'S_mean', 'G2_mean']:
    out[c] = out[c].round(6)

out = out.sort_values('FDR').reset_index(drop=True)

anova_path = os.path.join(OUTDIR, 'cellcycle_anova_results.csv')
out.to_csv(anova_path, index=False)
print(f'\nSaved: {anova_path}')
print(f'  Total Class I sites tested: {len(out)}')

n_resp = int(out['is_CC_responsive'].sum())
print(f'  CC responsive (FDR < 0.05): {n_resp} ({n_resp/len(out)*100:.1f}%)')
print(f'  Non-responsive: {len(out) - n_resp}')

# Phase distribution among responsive
resp = out[out['is_CC_responsive']]
phase_dist = resp['max_phase'].value_counts()
print(f'\n  Phase distribution (responsive sites, by max_phase):')
for phase, count in phase_dist.items():
    print(f'    {phase}: {count} ({count/n_resp*100:.1f}%)')

# Top 10 most significant
print(f'\n  Top 10 most significant sites:')
for _, row in out.head(10).iterrows():
    print(f'    {row["Gene"]:12s} pos={row["Position"]}  '
          f'F={row["F_stat"]:8.2f}  FDR={row["FDR"]:.2e}  '
          f'G1={row["G1_mean"]:+.3f} S={row["S_mean"]:+.3f} G2={row["G2_mean"]:+.3f}  '
          f'max={row["max_phase"]}')

# ============================================================
# Step 6-8: Cross-reference with CDK substrates
# ============================================================
print('\n' + '=' * 60)
print('CDK SUBSTRATE CROSS-REFERENCE')
print('=' * 60)

ks = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\kinase_substrate_correlation.csv')
ns = ks[ks['Kinase'] != ks['Substrate']].copy()
cdk_pairs = ns[ns['Kinase'].apply(lambda x: isinstance(x, str) and x.upper().startswith('CDK'))]
cdk_substrates = sorted(cdk_pairs['Substrate'].unique())
print(f'CDK unique substrates (from K-S correlation): {len(cdk_substrates)}')

# For each substrate, count sites in ANOVA results
summary_rows = []
for sub in cdk_substrates:
    sub_sites = out[out['Gene'] == sub]
    n_sites = len(sub_sites)
    n_responsive = int(sub_sites['is_CC_responsive'].sum()) if n_sites > 0 else 0
    pct = round(n_responsive / n_sites * 100, 1) if n_sites > 0 else np.nan
    summary_rows.append({
        'Substrate': sub,
        'N_sites_in_cellcycle': n_sites,
        'N_CC_responsive': n_responsive,
        'Pct_responsive': pct,
    })

summary = pd.DataFrame(summary_rows)
summary = summary.sort_values('N_CC_responsive', ascending=False).reset_index(drop=True)

summary_path = os.path.join(OUTDIR, 'cdk_substrates_cc_responsive_summary.csv')
summary.to_csv(summary_path, index=False)
print(f'Saved: {summary_path}')

n_with_data = int((summary['N_sites_in_cellcycle'] > 0).sum())
n_with_resp = int((summary['N_CC_responsive'] > 0).sum())
print(f'\nSummary:')
print(f'  CDK unique substrates:                {len(cdk_substrates)}')
print(f'  With >=1 site in cell-cycle data:     {n_with_data} ({n_with_data/len(cdk_substrates)*100:.1f}%)')
print(f'  With >=1 CC responsive site:          {n_with_resp} ({n_with_resp/len(cdk_substrates)*100:.1f}%)')

# Top 15 CDK substrates by responsive sites
print(f'\nTop 15 CDK substrates by N_CC_responsive:')
for _, row in summary.head(15).iterrows():
    if row['N_sites_in_cellcycle'] == 0:
        continue
    print(f'  {row["Substrate"]:12s}  sites={int(row["N_sites_in_cellcycle"]):3d}  '
          f'responsive={int(row["N_CC_responsive"]):3d}  ({row["Pct_responsive"]:.0f}%)')

print(f'\nAll outputs in: {OUTDIR}')
