# Anomaly-first probes on REAL data. Read-only; outputs -> anomaly_probe/.
import os, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
ROOT=r'D:\博士\Protein contour\Phospho'; OUT=os.path.join(ROOT,'anomaly_probe')
rng=np.random.default_rng(7)
TPS=['0min','2min','8min','20min','90min']; COMPS=['Cyt','Mem','Nuc']

# ============ PROBE 1: why is per-compartment co-variation/baseline so different? ============
# Nuc random baseline ~0.46, Cyt ~0.09 (from fig_corr_vs_random). Test: shared temporal trend.
pam=pd.read_csv(os.path.join(ROOT,'kinase_correlation','protein_abundance_matrix.csv')).drop_duplicates('Gene').set_index('Gene')
for c in pam.columns: pam[c]=pd.to_numeric(pam[c],errors='coerce')
print("=== PROBE 1: per-compartment temporal structure (explains baseline gap) ===")
print(f"{'comp':4s} {'n_var':>6s} {'medSD':>7s} {'PC1var%':>8s} {'randMed|r|':>10s} {'medMonotone':>11s}")
p1=[]
for comp in COMPS:
    cols=[f'{comp}_{t}' for t in TPS]
    M=pam[cols].values.astype(float)
    # keep proteins with >=3 finite & nonzero variance across the 5 timepoints
    ok=[]
    for row in M:
        m=np.isfinite(row)
        if m.sum()>=3 and np.nanstd(row[m])>0: ok.append(row)
    A=np.array(ok)
    sd=np.array([np.nanstd(r[np.isfinite(r)]) for r in A])
    # z-score each protein across timepoints (fill NaN with row-mean for PCA)
    Z=[]
    for r in A:
        m=np.isfinite(r); mu=r[m].mean(); s=r[m].std()
        rr=np.where(m,r,mu); Z.append((rr-mu)/(s if s>0 else 1))
    Z=np.array(Z)
    # PC1 variance ratio
    U,S,Vt=np.linalg.svd(Z-Z.mean(0),full_matrices=False)
    pc1=(S[0]**2)/np.sum(S**2)*100
    # direct random-pair median |pearson r| (sample 5000)
    idx=np.arange(len(A)); rs=[]
    for _ in range(5000):
        i,j=rng.integers(0,len(A),2)
        a,b=A[i],A[j]; m=np.isfinite(a)&np.isfinite(b)
        if m.sum()>=3 and a[m].std()>0 and b[m].std()>0:
            rs.append(np.corrcoef(a[m],b[m])[0,1])
    rs=np.array(rs)
    # fraction of proteins monotone over time (|spearman vs time|>0.9)
    mono=np.mean([abs(spearmanr(np.arange(m.sum()), r[m])[0])>0.9 for r,m in
                  [(r,np.isfinite(r)) for r in A]])
    p1.append((comp,len(A),np.median(sd),pc1,np.median(rs),mono))
    print(f"{comp:4s} {len(A):6d} {np.median(sd):7.3f} {pc1:8.1f} {np.median(np.abs(rs)):10.3f} {mono:11.3f}")
pd.DataFrame(p1,columns=['comp','n_var','medSD','PC1var_pct','rand_med_absr','frac_monotone']).to_csv(os.path.join(OUT,'probe1_compartment_structure.csv'),index=False)

# ============ PROBE 2: swap Δr within-family vs cross-family ============
print("\n=== PROBE 2: swap Δr by family (is 'kinase-agnostic' really 'family-agnostic'?) ===")
sm=pd.read_csv(os.path.join(ROOT,'martinez_network_check','swap_test_matrix.csv'),index_col=0)
fam={'CDK2':'CMGC','MAPK1':'CMGC','GSK3B':'CMGC','AKT1':'AGC','PRKACA':'AGC','CSNK2A1':'acidophilic','SRC':'Tyr'}
anchors={'CDK2 top5':'CDK2','MAPK1 top5':'MAPK1','AKT1 top5':'AKT1'}
within=[]; cross=[]; ck2=[]
for rowname,anch in anchors.items():
    orig=sm.loc[rowname,anch]
    for k in sm.columns:
        if k==anch: continue
        dr=orig-sm.loc[rowname,k]
        if k=='CSNK2A1': ck2.append(dr)
        elif fam[k]==fam[anch]: within.append(dr)
        else: cross.append(dr)
print(f"  within-family  Δr: median={np.median(within):+.3f} (n={len(within)})")
print(f"  cross-family(excl CK2) Δr: median={np.median(cross):+.3f} (n={len(cross)})")
print(f"  to CSNK2A1/CK2 Δr: median={np.median(ck2):+.3f} (n={len(ck2)})")

# ============ PROBE 3: scaffold feature collapse vs correlation with confounders ============
print("\n=== PROBE 3: does a feature's collapse-after-control track its correlation with disorder/length? ===")
pt=pd.read_csv(os.path.join(ROOT,'scaffold_check','per_protein_table.csv'))
pt['loglen']=np.log10(pt['length'])
feat_collapse={'Backbone':(-0.293,-0.160),'Sidechain':(-0.272,-0.166),'Coil':(0.254,0.127),
               'Sheet':(-0.196,-0.140),'Helix':(-0.147,-0.054),'EarlyFolding':(-0.207,-0.023),'Disorder':(0.115,0.020)}
print(f"{'feature':12s} {'raw_d':>7s} {'matched_d':>9s} {'retain%':>8s} {'|r_disorder|':>12s} {'|r_loglen|':>10s}")
p3=[]
for f,(raw,matched) in feat_collapse.items():
    col='avg_'+f
    rd=abs(spearmanr(pt[col],pt['pct_dis'],nan_policy='omit')[0])
    rl=abs(spearmanr(pt[col],pt['loglen'],nan_policy='omit')[0])
    retain=abs(matched)/abs(raw)*100 if raw else np.nan
    p3.append((f,raw,matched,retain,rd,rl))
    print(f"{f:12s} {raw:+7.3f} {matched:+9.3f} {retain:7.0f}% {rd:12.3f} {rl:10.3f}")
pd.DataFrame(p3,columns=['feature','raw_d','matched_d','retain_pct','abs_r_disorder','abs_r_loglen']).to_csv(os.path.join(OUT,'probe3_collapse_vs_confounder.csv'),index=False)

# ============ PROBE 4: movement_score vs biophysics at SITE level (why movers NS but presence sig?) ============
print("\n=== PROBE 4: site-level corr(movement_score, feature) — does mobility carry biophysical signal? ===")
H=pd.read_csv(os.path.join(ROOT,'high_mobility_sites.csv')); L=pd.read_csv(os.path.join(ROOT,'low_mobility_sites.csv'))
feats_site=['backbone_dynamics','sidechain_dynamics','disorder_propensity','helix_propensity','sheet_propensity','coil_propensity']
S=pd.concat([H,L],ignore_index=True)
print(f"  combined sites with movement_score: {S['movement_score'].notna().sum()}")
for f in feats_site:
    rho,p=spearmanr(S['movement_score'],S[f],nan_policy='omit')
    print(f"  {f:22s} rho={rho:+.3f} p={p:.2e}")

# ============ PROBE 5: scaffold threshold sensitivity (high-threshold sweep) ============
print("\n=== PROBE 5: scaffold backbone Cliff d under high-mobility threshold sweep (low=<5 fixed) ===")
g2u=json.load(open(os.path.join(ROOT,'gene_to_uniprot.json'))); BIO=os.path.join(ROOT,'biophys_json')
_bb={}
def backbone_mean(gene):
    if gene in _bb: return _bb[gene]
    uid=g2u.get(gene); v=np.nan
    if uid and os.path.exists(os.path.join(BIO,uid+'.json')):
        res=json.load(open(os.path.join(BIO,uid+'.json')))['residues']
        vals=[r.get('backbone') for r in res if r.get('backbone') is not None]
        v=np.mean(vals) if vals else np.nan
    _bb[gene]=v; return v
low_genes=set(L['Gene'].dropna())
for thr in [10,15,20,25]:
    hg=set(H[H['movement_score']>=thr]['Gene'].dropna())
    both=hg & low_genes; lowonly=low_genes-hg
    bvals=[backbone_mean(g) for g in both]; bvals=[x for x in bvals if np.isfinite(x)]
    lvals=[backbone_mean(g) for g in lowonly]; lvals=[x for x in lvals if np.isfinite(x)]
    U,p=mannwhitneyu(bvals,lvals,alternative='two-sided')
    d=2*U/(len(bvals)*len(lvals))-1
    print(f"  high>={thr:2d}: both={len(both):4d} static={len(lowonly):4d}  backbone Cliff d={d:+.3f} p={p:.1e}")
print("\nSaved probe1/probe3 CSVs in anomaly_probe/")
