import pandas as pd, numpy as np
from scipy.stats import spearmanr
import warnings; warnings.filterwarnings('ignore')

OUT = r'D:\博士\Protein contour\Phospho\martinez_network_check'
TPS=[0,2,8,20,90]; COMPS=['Cyt','Mem','Nuc']

pam=pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\protein_abundance_matrix.csv')
for c in pam.columns:
    if c!='Gene': pam[c]=pd.to_numeric(pam[c],errors='coerce')
pam=pam.set_index('Gene')

F=pd.read_csv(f'{OUT}/stage4_FINAL_r080.csv')   # 69 pair-specific QC-pass survivors
subs=F[['Substrate']].drop_duplicates()
print(f'FINAL r>=0.8 pair-specific QC-pass rows: {len(F)}; unique substrates: {len(subs)}')

def fracs(gene):
    """Return 5x3 array of compartment fractions per timepoint (rows=tp, cols=Cyt/Mem/Nuc)."""
    if gene not in pam.index: return None
    M=np.zeros((5,3))
    for i,t in enumerate(TPS):
        vals=np.array([pam.loc[gene,f'{c}_{t}min'] for c in COMPS],dtype=float)
        if np.any(np.isnan(vals)) or np.nansum(vals)<=0: return None
        M[i]=vals/vals.sum()
    return M

rows=[]
for g in subs['Substrate']:
    Fm=fracs(g)
    if Fm is None:
        rows.append({'Substrate':g,'reloc_index':np.nan,'note':'incomplete fractions'}); continue
    base=Fm[0]
    # reloc index = max over later timepoints of sum|frac(t)-frac(0)|  (total variation distance *2)
    tv=[np.abs(Fm[i]-base).sum() for i in range(1,5)]
    reloc=max(tv); tmax=TPS[1:][int(np.argmax(tv))]
    # direction: which compartment gains most at tmax
    delta=Fm[TPS.index(tmax)]-base
    gain=COMPS[int(np.argmax(delta))]; lose=COMPS[int(np.argmin(delta))]
    # also temporal SD of each compartment fraction, summed
    sd=np.sum(np.std(Fm,axis=0))
    rows.append({'Substrate':g,'reloc_index':round(reloc,3),'reloc_SD':round(sd,3),
                 'peak_time':tmax,'dir':f'{lose}->{gain}',
                 'Cyt0':round(base[0],2),'Mem0':round(base[1],2),'Nuc0':round(base[2],2)})
R=pd.DataFrame(rows)

# merge pair-level evidence (take strongest pair per substrate)
agg=F.groupby('Substrate').agg(r_residual=('r_residual','max'),motif=('motif_pct_max','max'),
                               compartment=('compartment','first')).reset_index()
R=R.merge(agg,on='Substrate',how='left')

# classify movers. Threshold: reloc_index>=0.20 => >=20% of total fraction mass shifted compartments
THR=0.20
R['mover']=R['reloc_index']>=THR
moved=R[R['mover']].sort_values('reloc_index',ascending=False)
print(f'\nReloc index = max_t sum_comp |frac(t)-frac(0)|  (0=no shift, 2=complete swap).')
print(f'Mover threshold: reloc_index >= {THR} (>=20% of fractional mass redistributes). Reason: 10% is the project mobility cutoff per compartment; summed across 3 comps ~20% marks a clear shift above noise.')
print(f'\nMovers: {R["mover"].sum()}/{R["reloc_index"].notna().sum()} substrates with complete fractions')
print('\n=== Substrates showing cross-compartment movement ===')
print(moved[['Substrate','reloc_index','peak_time','dir','Cyt0','Mem0','Nuc0','r_residual','motif']].to_string(index=False))

print('\nDirection summary (movers):')
print(moved['dir'].value_counts().to_string())
print('Peak-time summary (movers):')
print(moved['peak_time'].value_counts().sort_index().to_string())

# Q4 cross relationship
val=R.dropna(subset=['reloc_index','r_residual'])
rho_r,p_r=spearmanr(val['reloc_index'],val['r_residual'])
rho_m,p_m=spearmanr(val['reloc_index'],val['motif'])
print(f'\n=== Q4: does pair-specificity / motif relate to movement? (n={len(val)}) ===')
print(f'Spearman reloc_index vs r_residual: rho={rho_r:.3f} p={p_r:.3f}')
print(f'Spearman reloc_index vs motif_pct : rho={rho_m:.3f} p={p_m:.3f}')
R.to_csv(f'{OUT}/stage5_relocalization.csv',index=False)
print('\nSaved stage5_relocalization.csv')
