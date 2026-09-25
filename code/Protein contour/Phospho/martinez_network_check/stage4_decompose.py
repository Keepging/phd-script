import sys; sys.path.insert(0, r'D:\博士\Protein contour\Phospho\martinez_network_check')
import pandas as pd, numpy as np
from scipy.stats import pearsonr
import warnings; warnings.filterwarnings('ignore')

OUT = r'D:\博士\Protein contour\Phospho\martinez_network_check'
TPS = [0,2,8,20,90]
COMPS = ['Cyt','Mem','Nuc']

pam = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\protein_abundance_matrix.csv')
for c in pam.columns:
    if c!='Gene': pam[c]=pd.to_numeric(pam[c],errors='coerce')
pam=pam.set_index('Gene')

def raw(gene,comp):
    if gene not in pam.index: return None
    return pam.loc[gene,[f'{comp}_{t}min' for t in TPS]].values.astype(float)
def log2fc(gene,comp):
    v=raw(gene,comp)
    if v is None or np.isnan(v[0]) or v[0]<=0: return None
    out=np.full(5,np.nan)
    for i in range(5):
        if not np.isnan(v[i]) and v[i]>0: out[i]=np.log2(v[i]/v[0])
    return out

# ---------- A1: global compartment trajectories ----------
print('='*64); print('A1. GLOBAL COMPARTMENT TRAJECTORIES (median log2FC vs 0min)'); print('='*64)
GLOBAL={}
for comp in COMPS:
    mat=[]
    for g in pam.index:
        f=log2fc(g,comp)
        if f is not None: mat.append(f)
    mat=np.array(mat)
    med=np.nanmedian(mat,axis=0)
    GLOBAL[comp]=med
    print(f'{comp:3s} (n_proteins={len(mat)}): ' + '  '.join(f'{t}min={m:+.3f}' for t,m in zip(TPS,med)))
mt=GLOBAL['Mem']
trough=TPS[1:][np.argmin(mt[1:])]
print(f'\nMembrane shape check: trough at {trough}min (min log2FC={mt[1:].min():+.3f}), 90min={mt[4]:+.3f}')
print('=> membrane "2-8min drop, then partial recovery" CONFIRMED' if (mt[1]<0 and mt[2]<0 and trough in (2,8) and mt[4]>mt[1:].min()) else '=> shape differs, inspect')

# ---------- A2/A3: residualize + recompute r ----------
def pair_r(K,S,comp,residual):
    fk=log2fc(K,comp); fs=log2fc(S,comp)
    if fk is None or fs is None: return np.nan,0,np.nan,np.nan
    if residual:
        fk=fk-GLOBAL[comp]; fs=fs-GLOBAL[comp]
    mask=~(np.isnan(fk)|np.isnan(fs))
    n=int(mask.sum())
    if n<2 or np.nanstd(fk[mask])==0 or np.nanstd(fs[mask])==0: return np.nan,n,np.nanvar(fk[mask]),np.nanvar(fs[mask])
    r,_=pearsonr(fk[mask],fs[mask])
    return r,n,float(np.var(fk[mask])),float(np.var(fs[mask]))

surv=pd.read_csv(f'{OUT}/stage2_survivors_by_compartment.csv')  # r>=0.8 survivors
# we also need the r>=0.7 set: rebuild from pairs_with_motif
pwm=pd.read_csv(f'{OUT}/stage2_pairs_with_motif.csv')
pwm=pwm[pwm['motif_pct_max'].notna()].copy()

def build_survivors(rthr):
    rows=[]
    for comp,rcol in [('Cyt','Cyt_r'),('Mem','Mem_r'),('Nuc','Nuc_r')]:
        sub=pwm[(pwm[rcol]>=rthr)&(pwm['motif_pct_max']>=90)].copy()
        sub['compartment']=comp; sub['comp_r_pub']=sub[rcol]
        rows.append(sub)
    return pd.concat(rows,ignore_index=True)

def annotate(df):
    out=[]
    for _,r in df.iterrows():
        K,S,comp=r['Kinase'],r['Substrate'],r['compartment']
        r_raw,n_raw,vk,vs=pair_r(K,S,comp,False)
        r_res,n_res,_,_=pair_r(K,S,comp,True)
        att=r_raw-r_res if (not np.isnan(r_raw) and not np.isnan(r_res)) else np.nan
        # QC
        low_n = n_raw<=3
        sk=np.sqrt(vk) if not np.isnan(vk) else np.nan
        ss=np.sqrt(vs) if not np.isnan(vs) else np.nan
        low_var = (not np.isnan(sk) and sk<0.1) or (not np.isnan(ss) and ss<0.1)
        qc_pass = (n_raw>=4) and (not low_var)
        # decomposition class
        if np.isnan(r_res): dclass='unscorable'
        elif r_res>=0.6: dclass='pair-specific'
        elif r_res>=0.3: dclass='partial'
        else: dclass='global-rider'
        out.append({**r.to_dict(),'r_raw':r_raw,'r_residual':r_res,'attenuation':att,
                    'n_pts':n_raw,'kin_std':sk,'sub_std':ss,'low_n':low_n,'low_var':low_var,
                    'qc_pass':qc_pass,'decomp_class':dclass})
    return pd.DataFrame(out)

S08=annotate(build_survivors(0.8))
S07=annotate(build_survivors(0.7))
S08.to_csv(f'{OUT}/stage4_survivors_r080_decomposed.csv',index=False)
S07.to_csv(f'{OUT}/stage4_survivors_r070_decomposed.csv',index=False)

# ---------- B: QC focus AAK1->AP2M1 ----------
print('\n'+'='*64); print('B. CORRELATION QC METADATA'); print('='*64)
aak=S08[(S08['Kinase']=='AAK1')&(S08['Substrate']=='AP2M1')]
if len(aak):
    a=aak.iloc[0]
    print(f"AAK1->AP2M1 ({a['compartment']}): pub_r={a['comp_r_pub']:.3f} r_raw={a['r_raw']:.3f} "
          f"n={a['n_pts']} kin_std={a['kin_std']:.3f} sub_std={a['sub_std']:.3f} "
          f"low_n={a['low_n']} low_var={a['low_var']} -> QC {'PASS' if a['qc_pass'] else 'FAIL'}")
print('\nQC-failing survivors @r>=0.8 (n<=3 or low variance):')
qf=S08[~S08['qc_pass']][['Kinase','Substrate','compartment','comp_r_pub','n_pts','kin_std','sub_std','low_n','low_var']]
print(f'  {len(qf)}/{len(S08)} rows fail QC')
print(qf.round(3).to_string(index=False) if len(qf) else '  (none)')

# ---------- C: threshold comparison ----------
print('\n'+'='*64); print('C. THRESHOLD COMPARISON r>=0.8 vs r>=0.7'); print('='*64)
print(f'survivor-rows: r>=0.8 -> {len(S08)} | r>=0.7 -> {len(S07)}  (delta {len(S07)-len(S08)})')
for comp in COMPS:
    print(f'  {comp}: 0.8={len(S08[S08.compartment==comp])}  0.7={len(S07[S07.compartment==comp])}')
key=['Kinase','Substrate','compartment']
extra=S07.merge(S08[key],on=key,how='left',indicator=True)
extra=extra[extra['_merge']=='left_only']
print(f'\nPairs entering only at 0.7 (the loosening): {len(extra)}')
print('  decomp_class of the 0.7-only pairs:')
print('   '+extra['decomp_class'].value_counts().to_string().replace('\n','\n   '))
print(f'  of these, QC pass: {extra["qc_pass"].sum()}/{len(extra)}; pair-specific: {(extra.decomp_class=="pair-specific").sum()}')

# ---------- D: final ----------
def final(df):
    return df[(df['qc_pass'])&(df['decomp_class']=='pair-specific')].sort_values(['compartment','r_residual'],ascending=[True,False])
F08=final(S08); F07=final(S07)
print('\n'+'='*64); print('D. FINAL: real, pair-specific, QC-pass, EGF/spatial signals'); print('='*64)
cols=['Kinase','Substrate','compartment','comp_r_pub','r_raw','r_residual','attenuation','motif_pct_max','n_pts']
print(f'\n--- @ r>=0.8 : {len(F08)} pairs ---')
print(F08[cols].round(3).to_string(index=False))
print(f'\n--- @ r>=0.7 : {len(F07)} pairs (adds {len(F07)-len(F08)}) ---')
add=F07.merge(F08[key],on=key,how='left',indicator=True); add=add[add['_merge']=='left_only']
print('extra pairs gained at 0.7 that are real+pair-specific+QC-pass:')
print(add[cols].round(3).to_string(index=False) if len(add) else '  (none)')
F08.to_csv(f'{OUT}/stage4_FINAL_r080.csv',index=False); F07.to_csv(f'{OUT}/stage4_FINAL_r070.csv',index=False)
print('\nSaved stage4 decomposed + FINAL csvs.')
