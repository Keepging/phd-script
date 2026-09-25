import kinase_library as kl
import pandas as pd

ms = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_results\v2\kinase_matched_sites_v2.csv')
our_kinases = sorted(ms['Kinase Name'].dropna().astype(str).unique())
print(f'Unique kinase names in matched_sites_v2: {len(our_kinases)}')

st_atlas = set(kl.get_kinase_list(kin_type='ser_thr'))
ty_atlas = set(kl.get_kinase_list(kin_type='tyrosine'))
print(f'ST atlas: {len(st_atlas)}, Tyr atlas: {len(ty_atlas)}')

# Manual alias map (atlas non-gene-symbol names -> handled by mapping our gene symbol to atlas token)
ALIAS = {
    'ABL1': 'ABL', 'ABL2': 'ARG',
    'MAPK1': 'ERK2', 'MAPK3': 'ERK1', 'MAPK8': 'JNK1', 'MAPK9': 'JNK2', 'MAPK10': 'JNK3',
    'MAPK14': 'P38A', 'MAPK11': 'P38B', 'MAPK12': 'P38G', 'MAPK13': 'P38D',
    'MAPK7': 'ERK5', 'MAPK15': 'ERK7',
    'RPS6KB1': 'P70S6K', 'RPS6KB2': 'P70S6KB',
    'RPS6KA1': 'RSK1', 'RPS6KA2': 'RSK3', 'RPS6KA3': 'RSK2', 'RPS6KA4': 'MSK2', 'RPS6KA5': 'MSK1',
    'CHEK1': 'CHK1', 'CHEK2': 'CHK2',
    'MAPKAPK2': 'MAPKAPK2', 'MAPKAPK3': 'MAPKAPK3', 'MAPKAPK5': 'MK5',
    'PRKAA1': 'AMPKA1', 'PRKAA2': 'AMPKA2',
    'PRKACA': 'PKACA', 'PRKACB': 'PKACB', 'PRKACG': 'PKACG',
    'PRKCA': 'PKCA', 'PRKCB': 'PKCB', 'PRKCD': 'PKCD', 'PRKCE': 'PKCE', 'PRKCG': 'PKCG',
    'PRKCH': 'PKCH', 'PRKCI': 'PKCI', 'PRKCQ': 'PKCT', 'PRKCZ': 'PKCZ',
    'PRKD1': 'PKD1', 'PRKD2': 'PKD2', 'PRKD3': 'PKD3',
    'PRKG1': 'PKG1', 'PRKG2': 'PKG2',
    'MTOR': 'MTOR', 'PRKDC': 'DNAPK', 'SMG1': 'SMG1',
    'CSNK1A1': 'CK1A', 'CSNK1D': 'CK1D', 'CSNK1E': 'CK1E', 'CSNK1G1': 'CK1G1',
    'CSNK1G2': 'CK1G2', 'CSNK1G3': 'CK1G3', 'CSNK2A1': 'CK2A1', 'CSNK2A2': 'CK2A2',
    'CAMK2A': 'CAMK2A', 'CAMK2B': 'CAMK2B', 'CAMK2D': 'CAMK2D', 'CAMK2G': 'CAMK2G',
    'CAMK4': 'CAMK4', 'CAMK1': 'CAMK1A',
    'AURKA': 'AURA', 'AURKB': 'AURB', 'AURKC': 'AURC',
    'STK11': 'LKB1', 'MKNK1': 'MNK1', 'MKNK2': 'MNK2',
    'PAK1': 'PAK1', 'PAK2': 'PAK2', 'PAK4': 'PAK4',
    'GSK3A': 'GSK3A', 'GSK3B': 'GSK3B',
    'PDPK1': 'PDHK1', 'EEF2K': 'EEF2K',
    'EGFR': 'EGFR', 'SRC': 'SRC', 'FYN': 'FYN', 'MET': 'MET', 'FER': 'FER',
    'AKT1': 'AKT1', 'AKT2': 'AKT2', 'AKT3': 'AKT3',
    'ROCK1': 'ROCK1', 'ROCK2': 'ROCK2',
    'PLK1': 'PLK1', 'PLK2': 'PLK2', 'PLK3': 'PLK3',
    'PKMYT1': 'MYT1', 'MARK2': 'MARK2',
}

def map_to_atlas(name):
    """Return (atlas_name, kin_type, status)."""
    # direct ST
    if name in st_atlas:
        return name, 'ser_thr', 'direct'
    if name in ty_atlas:
        return name, 'tyrosine', 'direct'
    # alias
    if name in ALIAS:
        al = ALIAS[name]
        if al in st_atlas:
            return al, 'ser_thr', 'alias'
        if al in ty_atlas:
            return al, 'tyrosine', 'alias'
    return None, None, 'UNMATCHED'

rows = []
for k in our_kinases:
    a, t, s = map_to_atlas(k)
    rows.append({'our_name': k, 'atlas_name': a, 'kin_type': t, 'status': s})
mp = pd.DataFrame(rows)
matched = mp[mp['status'] != 'UNMATCHED']
unmatched = mp[mp['status'] == 'UNMATCHED']
print(f'\nMapped: {len(matched)}/{len(mp)}  (direct={sum(mp.status=="direct")}, alias={sum(mp.status=="alias")})')
print(f'UNMATCHED: {len(unmatched)}')
print('\n=== UNMATCHED kinase names (for manual review) ===')
for k in unmatched['our_name']:
    # how many sites use this kinase
    n = (ms['Kinase Name'] == k).sum()
    print(f'  {k}  ({n} site-rows)')

mp.to_csv(r'D:\博士\Protein contour\Phospho\martinez_network_check\kinase_name_mapping.csv', index=False)
print('\nSaved mapping -> kinase_name_mapping.csv')

# coverage in terms of site-rows
ms2 = ms.merge(mp, left_on='Kinase Name', right_on='our_name', how='left')
nrows_mapped = (ms2['status'] != 'UNMATCHED').sum()
print(f'\nSite-rows with a mappable kinase: {nrows_mapped}/{len(ms)} ({nrows_mapped/len(ms)*100:.1f}%)')
print('CDK site-rows:', ms['Kinase Name'].astype(str).str.upper().str.startswith('CDK').sum())
