import sys; sys.path.insert(0, r'D:\博士\Protein contour\Phospho\martinez_network_check')
import pandas as pd, numpy as np
import warnings; warnings.filterwarnings('ignore')

OUT = r'D:\博士\Protein contour\Phospho\martinez_network_check'

# --- pair-level motif from stage1 site scores ---
site = pd.read_csv(f'{OUT}/stage1_site_motif_scores.csv')
site = site.dropna(subset=['matched_percentile'])
pair_motif = site.groupby(['Kinase Name','Gene']).agg(
    motif_pct_max=('matched_percentile','max'),
    motif_pct_mean=('matched_percentile','mean'),
    motif_log2_max=('matched_log2','max'),
    n_sites_scored=('matched_percentile','size'),
    best_site=('Position', lambda s: s.iloc[0])
).reset_index().rename(columns={'Kinase Name':'Kinase','Gene':'Substrate'})
# attach the seq of the best (max-percentile) site
idx = site.loc[site.groupby(['Kinase Name','Gene'])['matched_percentile'].idxmax()]
bestseq = idx[['Kinase Name','Gene','Position','AA','seq15']].rename(
    columns={'Kinase Name':'Kinase','Gene':'Substrate','Position':'best_site_pos','AA':'best_site_aa','seq15':'best_seq'})
pair_motif = pair_motif.merge(bestseq, on=['Kinase','Substrate'], how='left')

# --- correlation table (per-compartment) ---
corr = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\kinase_substrate_correlation.csv')
corr = corr[corr['Kinase'] != corr['Substrate']].copy()

# --- Step 4: remove CDK ---
is_cdk = corr['Kinase'].astype(str).str.upper().str.startswith('CDK')
corr_nocdk = corr[~is_cdk].copy()
print('='*60); print('STAGE 2'); print('='*60)
print(f'Correlation non-self pairs: {len(corr)}  (CDK rows: {is_cdk.sum()})')
print(f'After removing CDK: {len(corr_nocdk)} pairs, {corr_nocdk["Kinase"].nunique()} unique kinases')

# --- Step 5: join motif ---
merged = corr_nocdk.merge(pair_motif, on=['Kinase','Substrate'], how='left')
has_motif = merged['motif_pct_max'].notna()
print(f'\nPairs with a motif score (joined): {has_motif.sum()}/{len(merged)} ({has_motif.mean()*100:.1f}%)')
print('  (pairs without motif = kinase unmatched to atlas, or no scored site for that pair)')
merged.to_csv(f'{OUT}/stage2_pairs_with_motif.csv', index=False)

# --- thresholds (stated) ---
R_THR = 0.80     # strong per-compartment co-abundance (conventional "strong" Pearson)
PCT_THR = 90.0   # Kinase Library standard cutoff for a high-confidence predicted substrate
print(f'\nThresholds: per-compartment r >= {R_THR} (strong co-abundance) AND motif_pct_max >= {PCT_THR} (atlas high-confidence cutoff)')

m = merged[has_motif].copy()
survivors_all = []
for comp, rcol in [('Cyt','Cyt_r'),('Mem','Mem_r'),('Nuc','Nuc_r')]:
    sub = m[(m[rcol] >= R_THR) & (m['motif_pct_max'] >= PCT_THR)].copy()
    sub['compartment'] = comp
    sub['comp_r'] = sub[rcol]
    # combined score: equal-weight blend of correlation and motif (both 0-1)
    sub['combined'] = 0.5*sub['comp_r'] + 0.5*(sub['motif_pct_max']/100.0)
    sub = sub.sort_values('combined', ascending=False)
    survivors_all.append(sub)
    print(f'\n--- {comp}: {len(sub)} survivors (r>={R_THR} & motif>={PCT_THR}) ---')
    cols = ['Kinase','Substrate','comp_r','motif_pct_max','best_site_aa','best_site_pos','combined']
    print(sub.head(12)[cols].round(3).to_string(index=False))

surv = pd.concat(survivors_all, ignore_index=True)
surv.to_csv(f'{OUT}/stage2_survivors_by_compartment.csv', index=False)
print(f'\nSaved stage2_survivors_by_compartment.csv ({len(surv)} rows across compartments)')
print('Unique survivor pairs:', surv[['Kinase','Substrate']].drop_duplicates().shape[0])
