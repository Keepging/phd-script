import sys; sys.path.insert(0, r'D:\博士\Protein contour\Phospho\martinez_network_check')
import pandas as pd, numpy as np, kinase_library as kl
from kl_core import cut_window, map_to_atlas
import warnings; warnings.filterwarnings('ignore')

OUT = r'D:\博士\Protein contour\Phospho\martinez_network_check'
ms = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_results\v2\kinase_matched_sites_v2.csv')
print(f'matched_sites_v2 rows: {len(ms)}')

# --- map kinase + cut window ---
maprows = ms['Kinase Name'].astype(str).apply(lambda k: pd.Series(map_to_atlas(k), index=['atlas_name','kin_type','map_status']))
ms = pd.concat([ms, maprows], axis=1)
ms['seq15'], ms['center_aa'] = zip(*ms.apply(lambda r: cut_window(r['uniprot_id'], int(r['Position'])), axis=1))

# sanity: center should equal AA
mism = ms[(ms['seq15'].notna()) & (ms['center_aa'] != ms['AA'])]
print(f'center-AA mismatches: {len(mism)} (printed below if any)')
if len(mism): print(mism[['Gene','AA','Position','center_aa','seq15']].head(20).to_string())
noseq = ms['seq15'].isna().sum()
print(f'rows with no sequence window: {noseq}')

# residue routing summary
print('\nphospho-residue routing:')
print(ms['AA'].value_counts().to_string())

# --- batch score ST and TY separately ---
def score_block(block, kin_type):
    if len(block)==0: return None,None
    df = block[['seq15']].rename(columns={'seq15':'SITE'})
    pps = kl.PhosphoProteomics(df, seq_col='SITE', suppress_warnings=True)
    pct = pps.percentile(kin_type=kin_type)
    log2 = pps.score(kin_type=kin_type)
    return pct, log2

ms_valid = ms[ms['seq15'].notna() & ms['kin_type'].notna()].copy()
matched_pct, matched_log2 = [], []
for kin_type in ['ser_thr','tyrosine']:
    block = ms_valid[ms_valid['kin_type']==kin_type]
    pct, log2 = score_block(block, kin_type)
    if pct is None: continue
    pct = pct.reset_index(drop=True); log2 = log2.reset_index(drop=True)
    block = block.reset_index()
    for i,(_,r) in enumerate(block.iterrows()):
        a = r['atlas_name']
        try:
            mp = pct.iloc[i][a] if a in pct.columns else np.nan
            ml = log2.iloc[i][a] if a in log2.columns else np.nan
        except Exception:
            mp, ml = np.nan, np.nan
        matched_pct.append((r['index'], mp)); matched_log2.append((r['index'], ml))

pct_map = dict(matched_pct); log2_map = dict(matched_log2)
ms['matched_percentile'] = ms.index.map(lambda i: pct_map.get(i, np.nan))
ms['matched_log2'] = ms.index.map(lambda i: log2_map.get(i, np.nan))
ms['is_cdk'] = ms['Kinase Name'].astype(str).str.upper().str.startswith('CDK')

# save per-site table
keep = ['PTM_collapse_key','Gene','AA','Position','uniprot_id','Kinase','Kinase Name',
        'kinase_family','atlas_name','kin_type','map_status','seq15','center_aa',
        'matched_percentile','matched_log2','is_cdk']
ms[keep].to_csv(f'{OUT}/stage1_site_motif_scores.csv', index=False)
print(f'\nSaved stage1_site_motif_scores.csv ({len(ms)} rows)')

# coverage
scored = ms['matched_percentile'].notna()
print(f'rows scored (kinase in atlas + seq found): {scored.sum()}/{len(ms)} ({scored.mean()*100:.1f}%)')
print(f'  unmatched-kinase rows (annotated, NaN score): {(ms.map_status=="UNMATCHED").sum()}')

# ============ QC: CDK vs non-CDK ============
sc = ms[scored].copy()
cdk = sc[sc['is_cdk']]; noncdk = sc[~sc['is_cdk']]
print('\n' + '='*60)
print('STAGE 1 QC: CDK vs non-CDK matched-kinase motif percentile')
print('='*60)
def dist(x):
    return f'n={len(x):4d}  median={x.median():5.1f}  mean={x.mean():5.1f}  Q25={x.quantile(.25):5.1f}  Q75={x.quantile(.75):5.1f}  %>=90={ (x>=90).mean()*100:4.1f}%'
print('CDK     :', dist(cdk['matched_percentile']))
print('non-CDK :', dist(noncdk['matched_percentile']))

# proline-directed check: fraction of CDK sites with +1 == P
def plus1(seq):
    return seq[8] if isinstance(seq,str) and len(seq)>8 else None
cdk['plus1'] = cdk['seq15'].apply(plus1)
noncdk['plus1'] = noncdk['seq15'].apply(plus1)
print(f"\nCDK sites with Proline at +1: {(cdk['plus1']=='P').mean()*100:.1f}%  (n={len(cdk)})")
print(f"non-CDK sites with Proline at +1: {(noncdk['plus1']=='P').mean()*100:.1f}%")
print(f"CDK sites +1=P median percentile: {cdk[cdk['plus1']=='P']['matched_percentile'].median():.1f}")
print(f"CDK sites +1!=P median percentile: {cdk[cdk['plus1']!='P']['matched_percentile'].median():.1f}")

from scipy.stats import mannwhitneyu
u,p = mannwhitneyu(cdk['matched_percentile'], noncdk['matched_percentile'], alternative='greater')
print(f'\nMann-Whitney U (CDK > non-CDK percentile): p={p:.2e}')
print('\nQC VERDICT:', 'PASS - CDK enriched for high proline-directed motif scores' if cdk['matched_percentile'].median() > noncdk['matched_percentile'].median() and p<0.05 else 'FAIL - investigate')
