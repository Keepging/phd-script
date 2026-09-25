# Step-2 scaffold finding: confounder control (length + overall disorder).
# READ-ONLY on project data. Outputs -> scaffold_check/. Reconstructs the per-protein
# table exactly as 0331.py did (it was never saved to CSV).
import os, json, csv
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

ROOT = r'D:\博士\Protein contour\Phospho'
OUT  = os.path.join(ROOT, 'scaffold_check')
BIO  = os.path.join(ROOT, 'biophys_json')
rng  = np.random.default_rng(7)

g2u = json.load(open(os.path.join(ROOT,'gene_to_uniprot.json')))
H = pd.read_csv(os.path.join(ROOT,'high_mobility_sites.csv'))
L = pd.read_csv(os.path.join(ROOT,'low_mobility_sites.csv'))
high_genes = set(H['Gene'].dropna()); low_genes = set(L['Gene'].dropna())
both = high_genes & low_genes                 # mobile/hub proteins
low_only = low_genes - high_genes             # static-only proteins
print(f'both(mobile/hub)={len(both)}  low_only(static)={len(low_only)}')

JKEYS = {'backbone':'Backbone','sidechain':'Sidechain','disoMine':'Disorder',
         'helix':'Helix','sheet':'Sheet','coil':'Coil','earlyFolding':'EarlyFolding'}
FEATS = list(JKEYS.values())

rows=[]
for gene in (both | low_only):
    uid = g2u.get(gene)
    if not uid: continue
    fp = os.path.join(BIO, f'{uid}.json')
    if not os.path.exists(fp): continue
    res = json.load(open(fp))['residues']
    if not res: continue
    e = {'gene':gene,'group':'mobile' if gene in both else 'static',
         'length':len(res)}
    for jk,lab in JKEYS.items():
        v=[r.get(jk) for r in res if r.get(jk) is not None]
        e['avg_'+lab]=np.mean(v) if v else np.nan
        if jk=='disoMine':
            e['pct_dis']=100*sum(1 for x in v if x>0.5)/len(res)
    rows.append(e)
df=pd.DataFrame(rows)
df.to_csv(os.path.join(OUT,'per_protein_table.csv'),index=False)
print('reconstructed proteins:', len(df), df['group'].value_counts().to_dict())

def cliffs(b,l):
    # delta = 2*U/(n1 n2) - 1 ; positive => mobile > static
    U,p = mannwhitneyu(b,l,alternative='two-sided')
    d = 2*U/(len(b)*len(l)) - 1
    return d,p

mob=df[df.group=='mobile']; sta=df[df.group=='static']
print('\n=== CONFOUNDERS THEMSELVES (mobile vs static) ===')
for c in ['length','pct_dis']:
    b=mob[c].dropna(); l=sta[c].dropna(); d,p=cliffs(b.values,l.values)
    print(f'  {c:10s} mobile_med={np.median(b):.1f} static_med={np.median(l):.1f} Cliff_d={d:+.3f} p={p:.2e}')

# ---------- RAW ----------
print('\n=== RAW vs CONTROLLED (Cliff delta; + = mobile>static) ===')
print(f"{'feature':14s} {'raw_d':>8s} {'raw_p':>10s} | {'resid_d':>8s} {'resid_p':>10s} | {'match_d':>8s} {'match_p':>10s} {'perm_p':>8s}")
results=[]

# residualization design: control length(log)+pct_dis (disorder feature -> length only)
df['loglen']=np.log10(df['length'])
def residualize(col, ctrl):
    sub=df[[col]+ctrl].dropna()
    X=np.column_stack([np.ones(len(sub))]+[sub[c].values for c in ctrl])
    y=sub[col].values
    beta,_,_,_=np.linalg.lstsq(X,y,rcond=None)
    r=y-X@beta
    out=pd.Series(np.nan,index=df.index); out.loc[sub.index]=r
    return out

# greedy nearest-neighbor match on standardized [loglen,pct_dis]
def matched_indices():
    feats=['loglen','pct_dis']
    d2=df.dropna(subset=feats)
    z=(d2[feats]-d2[feats].mean())/d2[feats].std()
    mob_idx=d2[d2.group=='mobile'].index.tolist()
    sta_pool=d2[d2.group=='static'].index.tolist()
    zsta=z.loc[sta_pool].values; used=set(); pairs=[]
    for mi in mob_idx:
        zv=z.loc[mi].values
        dist=np.sqrt(((zsta-zv)**2).sum(1))
        order=np.argsort(dist)
        for oi in order:
            si=sta_pool[oi]
            if si in used: continue
            if dist[oi]>0.25: break          # caliper
            used.add(si); pairs.append((mi,si)); break
    return [m for m,_ in pairs],[s for _,s in pairs]

mob_m, sta_m = matched_indices()
print(f'(matched pairs within caliper: {len(mob_m)})')

for lab in FEATS:
    col='avg_'+lab
    b=df.loc[df.group=='mobile',col].dropna().values
    l=df.loc[df.group=='static',col].dropna().values
    raw_d,raw_p=cliffs(b,l)
    # residualized
    ctrl=['loglen'] if lab=='Disorder' else ['loglen','pct_dis']
    rcol=residualize(col,ctrl)
    rb=rcol[df.group=='mobile'].dropna().values
    rl=rcol[df.group=='static'].dropna().values
    res_d,res_p=cliffs(rb,rl)
    # matched
    mb=df.loc[mob_m,col].dropna().values
    ml=df.loc[sta_m,col].dropna().values
    mat_d,mat_p=cliffs(mb,ml)
    # permutation null on raw Cliff d
    allv=np.concatenate([b,l]); n1=len(b)
    null=[]
    for _ in range(1000):
        rng.shuffle(allv)
        d,_=cliffs(allv[:n1],allv[n1:]); null.append(d)
    null=np.array(null); perm_p=(np.sum(np.abs(null)>=abs(raw_d))+1)/(len(null)+1)
    results.append({'feature':lab,'raw_d':raw_d,'raw_p':raw_p,'resid_d':res_d,
                    'resid_p':res_p,'match_d':mat_d,'match_p':mat_p,'perm_p':perm_p,'ctrl':'+'.join(ctrl)})
    print(f"{lab:14s} {raw_d:+8.3f} {raw_p:10.2e} | {res_d:+8.3f} {res_p:10.2e} | {mat_d:+8.3f} {mat_p:10.2e} {perm_p:8.3f}  ctrl={'+'.join(ctrl)}")

pd.DataFrame(results).to_csv(os.path.join(OUT,'scaffold_confounder_results.csv'),index=False)
print('\nSaved per_protein_table.csv + scaffold_confounder_results.csv in scaffold_check/')
print('Cliff delta guide: |0.11| small, |0.28| medium, |0.43| large')
