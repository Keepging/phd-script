import sys, json, os
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import warnings; warnings.filterwarnings('ignore')

OUT=r'D:\博士\Protein contour\Phospho\martinez_network_check'
JSON=r'D:\博士\Protein contour\Phospho\biophys_json'
TPS=[0,2,8,20,90]; COMPS=['Cyt','Mem','Nuc']
FEATS=['backbone','sidechain','disoMine','helix','sheet','coil','earlyFolding']

pam=pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\protein_abundance_matrix.csv')
for c in pam.columns:
    if c!='Gene': pam[c]=pd.to_numeric(pam[c],errors='coerce')
pam=pam.set_index('Gene')
ms=pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_results\v2\kinase_matched_sites_v2.csv')

def fracs(gene):
    if gene not in pam.index: return None
    M=np.zeros((5,3))
    for i,t in enumerate(TPS):
        v=np.array([pam.loc[gene,f'{c}_{t}min'] for c in COMPS],dtype=float)
        if np.any(np.isnan(v)) or np.nansum(v)<=0: return None
        M[i]=v/v.sum()
    return M

# ---- 1. global fraction trajectory (all complete proteins) ----
allF=[]; genes_complete=[]
for g in pam.index:
    Fm=fracs(g)
    if Fm is not None: allF.append(Fm); genes_complete.append(g)
allF=np.array(allF)              # N x 5 x 3
GFRAC=np.median(allF,axis=0)     # 5 x 3 global median fraction trajectory
print('Global median fraction trajectory (rows=tp, cols=Cyt/Mem/Nuc):')
for i,t in enumerate(TPS): print(f'  {t:>2}min: '+'  '.join(f'{c}={GFRAC[i,j]:.3f}' for j,c in enumerate(COMPS)))

def reloc(Fm, residual):
    if Fm is None: return np.nan, None
    X=Fm-GFRAC if residual else Fm
    tv=[np.abs(X[i]-X[0]).sum() for i in range(1,5)]
    k=int(np.argmax(tv));
    d=X[TPS.index(TPS[1:][k])]-X[0]
    return max(tv), f'{COMPS[int(np.argmin(d))]}->{COMPS[int(np.argmax(d))]}'

# ---- pool = phosphoproteins with json + complete fractions ----
pool=[]
for g,acc in ms[['Gene','uniprot_id']].drop_duplicates().itertuples(index=False):
    if g not in pam.index: continue
    if not os.path.exists(os.path.join(JSON,f'{acc}.json')): continue
    Fm=fracs(g)
    if Fm is None: continue
    rr,dr=reloc(Fm,True); rraw,draw=reloc(Fm,False)
    pool.append({'Gene':g,'acc':acc,'reloc_resid':rr,'dir_resid':dr,'reloc_raw':rraw,'dir_raw':draw})
P=pd.DataFrame(pool).dropna(subset=['reloc_resid'])
print(f'\nPhospho-protein pool with complete data: {len(P)}')

# ---- 2. mover/non-mover on reloc_resid (tertiles) ----
q1,q2=P['reloc_resid'].quantile([1/3,2/3])
P['grp']=np.where(P['reloc_resid']>=q2,'mover',np.where(P['reloc_resid']<=q1,'non-mover','mid'))
print(f'reloc_resid tertile cuts: low<= {q1:.3f}, high>= {q2:.3f}  (raw reloc median {P.reloc_raw.median():.3f} -> resid median {P.reloc_resid.median():.3f})')
print('Direction one-sidedness AFTER de-global (movers):')
print(P[P.grp=="mover"]['dir_resid'].value_counts().to_string())
print('(compare raw direction, movers):')
print(P[P.grp=="mover"]['dir_raw'].value_counts().to_string())
mov=P[P.grp=='mover']; non=P[P.grp=='non-mover']
print(f'\nGroups: mover n={len(mov)}, non-mover n={len(non)} (tertiles; reason: balanced high/low-contrast groups from one phospho pool)')

# ---- 4. biophysics features ----
_cache={}
def feats(acc):
    if acc in _cache: return _cache[acc]
    d=json.load(open(os.path.join(JSON,f'{acc}.json')))
    res=sorted(d['residues'],key=lambda x:x['seqpos'])
    arr={f:np.array([r.get(f,np.nan) for r in res],dtype=float) for f in FEATS}
    _cache[acc]=(arr,len(res)); return _cache[acc]

def protein_feats(g,acc):
    arr,L=feats(acc)
    full={f'full_{f}':np.nanmean(arr[f]) for f in FEATS}
    # local +/-5 around each phospho site
    pos=ms[ms['Gene']==g]['Position'].dropna().astype(int).unique()
    locvals={f:[] for f in FEATS}
    for p in pos:
        lo,hi=max(0,p-1-5),min(L,p+5)
        for f in FEATS: locvals[f].append(np.nanmean(arr[f][lo:hi]))
    local={f'local_{f}':(np.nanmean(locvals[f]) if locvals[f] else np.nan) for f in FEATS}
    ab=np.nanmean(pam.loc[g].values.astype(float))
    return {**full,**local,'length':L,'n_sites':len(pos),'abundance':ab}

rows=[]
for _,r in P[P.grp.isin(['mover','non-mover'])].iterrows():
    rows.append({'Gene':r['Gene'],'grp':r['grp'],**protein_feats(r['Gene'],r['acc'])})
B=pd.DataFrame(rows)

# ---- 6a. confounds first ----
print('\n=== Confound check (mover vs non-mover) ===')
for c in ['length','abundance','n_sites']:
    a=B[B.grp=='mover'][c]; b=B[B.grp=='non-mover'][c]
    u,p=mannwhitneyu(a,b)
    print(f'  {c:10s}: mover med={a.median():.1f}  non med={b.median():.1f}  MW p={p:.3f}')

# ---- 5. feature comparison ----
def rbc(a,b):  # rank-biserial effect size
    u,_=mannwhitneyu(a,b,alternative='two-sided'); return 1-2*u/(len(a)*len(b))
def compare(df,cols,label):
    res=[]
    for c in cols:
        a=df[df.grp=='mover'][c].dropna(); b=df[df.grp=='non-mover'][c].dropna()
        u,p=mannwhitneyu(a,b,alternative='two-sided')
        res.append({'feature':c,'mover_med':a.median(),'non_med':b.median(),
                    'effect_rbc':rbc(a,b),'p':p})
    R=pd.DataFrame(res)
    R['p_BH']=multipletests(R['p'],method='fdr_bh')[1]
    R['sig']=R['p_BH']<0.05
    print(f'\n=== {label}: mover vs non-mover (BH-corrected) ===')
    print(R.round(4).to_string(index=False))
    return R
cols=[f'local_{f}' for f in FEATS]+[f'full_{f}' for f in FEATS]
Rfull=compare(B,cols,'ALL features')

# ---- 6b. length+abundance matched subanalysis ----
print('\n=== Length+abundance-matched subanalysis (1:1 nearest, 300 iters) ===')
rng=np.random.default_rng(0)
movb=B[B.grp=='mover'].copy(); nonb=B[B.grp=='non-mover'].copy()
# z-scale length & log-abundance for matching distance
for c in ['length','abundance']:
    mu,sd=B[c].mean(),B[c].std()
    movb[c+'_z']=(movb[c]-mu)/sd; nonb[c+'_z']=(nonb[c]-mu)/sd
sigcols=[c for c in cols if Rfull.set_index('feature').loc[c,'sig']]
if not sigcols: sigcols=cols  # if none sig, test all anyway
agg={c:[] for c in sigcols}
for _ in range(300):
    picks=[]
    pool_non=nonb.sample(frac=1,random_state=int(rng.integers(1e9)))
    used=set()
    for _,m in movb.iterrows():
        d=((pool_non['length_z']-m['length_z'])**2+(pool_non['abundance_z']-m['abundance_z'])**2)
        d=d[~pool_non.index.isin(used)]
        if len(d)==0: break
        j=d.idxmin(); used.add(j); picks.append(j)
    matched=nonb.loc[picks]
    for c in sigcols:
        a=movb[c].dropna(); b=matched[c].dropna()
        if len(a)>2 and len(b)>2:
            _,p=mannwhitneyu(a,b,alternative='two-sided'); agg[c].append(p)
print('Matched non-movers length med vs mover length med (last iter): '
      f'{matched["length"].median():.0f} vs {movb["length"].median():.0f}')
for c in sigcols:
    ps=np.array(agg[c]);
    print(f'  {c:22s}: median matched p={np.median(ps):.4f}  frac sig(<0.05)={np.mean(ps<0.05):.2f}')

B.to_csv(f'{OUT}/stage6_biophysics_features.csv',index=False)
Rfull.to_csv(f'{OUT}/stage6_feature_comparison.csv',index=False)
P.to_csv(f'{OUT}/stage6_reloc_resid_pool.csv',index=False)
print('\nSaved stage6_*.csv')
