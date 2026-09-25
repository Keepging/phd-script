# %% Step 2 — classify + Step 3 — attach earlyFolding from json
import pandas as pd, numpy as np, json, os
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

JSON_DIR = r'D:\博士\Protein contour\Phospho\biophys_json'
df = pd.read_csv(r'D:\博士\Protein contour\Phospho\success_df_M1M2_cleaned.csv')
print(f'Loaded: {len(df)} sites')

# Step 2: LocalizationGroup already encodes the 7 classes
groups = df['LocalizationGroup'].value_counts().sort_values(ascending=False)
print('\n=== Step 2: Presence classes (LocalizationGroup) ===')
for g, n in groups.items():
    print(f'  {g:8s}: {n:5d} ({n/len(df)*100:4.1f}%)')

# Step 3: supplement earlyFolding from json
# existing 6 features in df
FEATS_IN_DF = ['backbone_dynamics','sidechain_dynamics','disorder_propensity',
               'helix_propensity','sheet_propensity','coil_propensity']
# need earlyFolding from json at exact site position

_cache = {}
def load_json(acc):
    if acc in _cache: return _cache[acc]
    p = os.path.join(JSON_DIR, f'{acc}.json')
    if not os.path.exists(p):
        _cache[acc] = None; return None
    d = json.load(open(p))
    # build pos -> residue dict
    res = {r['seqpos']: r for r in d['residues']}
    _cache[acc] = res; return res

ef_vals = []
for _, row in df.iterrows():
    acc = row['uniprot_id']
    pos = int(row['Position'])
    res = load_json(acc) if pd.notna(acc) else None
    if res and pos in res:
        ef_vals.append(res[pos].get('earlyFolding', np.nan))
    else:
        ef_vals.append(np.nan)

df['earlyFolding'] = ef_vals

ALL_FEATS = FEATS_IN_DF + ['earlyFolding']
print('\n=== Step 3: Feature coverage (non-null counts) ===')
for f in ALL_FEATS:
    nn = df[f].notna().sum()
    print(f'  {f:25s}: {nn:5d}/{len(df)} ({nn/len(df)*100:.1f}%)')

# %% Step 4 — Kruskal-Wallis across 7 groups, BH correction, epsilon-squared
GROUP_ORDER = ['C-only','M-only','N-only','C&M','C&N','M&N','C&M&N']
N = len(df)

stats_rows = []
for feat in ALL_FEATS:
    sub = df[['LocalizationGroup', feat]].dropna()
    samples = [sub[sub['LocalizationGroup']==g][feat].values for g in GROUP_ORDER]
    # only groups with data
    valid = [(g, s) for g, s in zip(GROUP_ORDER, samples) if len(s) > 0]
    H, p = kruskal(*[s for _, s in valid])
    # epsilon-squared = (H - k + 1) / (n - k)
    k = len(valid)
    n_total = sum(len(s) for _, s in valid)
    eps2 = (H - k + 1) / (n_total - k) if n_total > k else np.nan
    medians = {g: np.median(s) for g, s in valid}
    stats_rows.append({'feature': feat, 'H': H, 'p': p, 'eps2': eps2,
                       'n_tested': n_total, **{f'med_{g}': medians.get(g, np.nan) for g in GROUP_ORDER}})

S = pd.DataFrame(stats_rows)
_, S['p_BH'], _, _ = multipletests(S['p'], method='fdr_bh')
S['sig'] = S['p_BH'] < 0.05
S = S[['feature','H','p','p_BH','sig','eps2','n_tested'] +
      [f'med_{g}' for g in GROUP_ORDER]]

print('\n=== Step 4: Kruskal-Wallis stats (7 features × 7 groups, BH-corrected) ===')
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 200)
pd.set_option('display.float_format', lambda x: f'{x:.4f}')
print(S.to_string(index=False))

# save intermediates
OUT = r'D:\博士\Protein contour\Phospho\martinez_network_check'
os.makedirs(os.path.join(OUT, 'outputs'), exist_ok=True)
os.makedirs(os.path.join(OUT, 'figures'), exist_ok=True)
df[['PTM_collapse_key','Gene','Site','AA','Position','uniprot_id',
    'LocalizationGroup'] + ALL_FEATS].to_csv(
    os.path.join(OUT, 'outputs', 'presence_pattern_sites.csv'), index=False)
S.to_csv(os.path.join(OUT, 'outputs', 'presence_pattern_stats.csv'), index=False)
print('\nSaved outputs/presence_pattern_sites.csv and outputs/presence_pattern_stats.csv')
print('\n*** STOPPING HERE — review stats before generating figure ***')
