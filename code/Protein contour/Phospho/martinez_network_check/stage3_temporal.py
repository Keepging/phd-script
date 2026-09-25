import sys; sys.path.insert(0, r'D:\博士\Protein contour\Phospho\martinez_network_check')
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings('ignore')

OUT = r'D:\博士\Protein contour\Phospho\martinez_network_check'
def save_tiff(fig, fn, dpi=300):
    p = fn.rsplit('.',1)[0]+'.tiff'; fig.savefig(p, dpi=dpi, bbox_inches='tight', format='tiff', pil_kwargs={'compression':'tiff_lzw'})

pam = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\protein_abundance_matrix.csv')
for c in pam.columns:
    if c!='Gene': pam[c]=pd.to_numeric(pam[c], errors='coerce')
pam = pam.set_index('Gene')
TPS=[0,2,8,20,90]
def prof(gene, comp):
    if gene not in pam.index: return None
    return pam.loc[gene, [f'{comp}_{t}min' for t in TPS]].values.astype(float)

def log2fc(p):
    if p is None or np.isnan(p[0]) or p[0]<=0: return None
    return np.log2(p/p[0])

surv = pd.read_csv(f'{OUT}/stage2_survivors_by_compartment.csv')
print('='*60); print('STAGE 3: temporal structure of survivors'); print('='*60)

rows=[]
for _,r in surv.iterrows():
    comp=r['compartment']; K=r['Kinase']; S=r['Substrate']
    sp=prof(S,comp); kp=prof(K,comp)
    sfc=log2fc(sp); kfc=log2fc(kp)
    # substrate temporal descriptors in this compartment
    if sfc is not None:
        absfc=np.abs(sfc[1:])  # at 2,8,20,90
        maxfc=absfc.max();
        fc8=abs(sfc[2])  # 8min
        early=(fc8 >= 0.5*maxfc) and (maxfc>=0.3)  # half the change by 8min & meaningful (>=0.3 log2 ~1.23x)
        direction='up' if sfc[np.argmax(absfc)+1]>0 else 'down'
    else:
        maxfc=np.nan; early=False; direction='NA'
    # membrane-specificity: substrate max|log2fc| across compartments
    rng={}
    for cc in ['Cyt','Mem','Nuc']:
        f=log2fc(prof(S,cc)); rng[cc]=np.abs(f[1:]).max() if f is not None else np.nan
    memspec = (not np.isnan(rng['Mem'])) and (rng['Mem']>=rng['Cyt']) and (rng['Mem']>=rng['Nuc'])
    rows.append({**r.to_dict(),'sub_maxlog2fc':maxfc,'sub_dir':direction,'early_2_8min':early,
                 'mem_range':rng['Mem'],'cyt_range':rng['Cyt'],'nuc_range':rng['Nuc'],'mem_specific':memspec})
S3=pd.DataFrame(rows)

# classify
def cls(r):
    if np.isnan(r['sub_maxlog2fc']) or r['sub_maxlog2fc']<0.3: return 'stable_background'
    if r['compartment']=='Mem' and r['mem_specific'] and r['early_2_8min']: return 'EGF-like (mem+early)'
    if r['early_2_8min'] and r['sub_maxlog2fc']>=0.3: return 'early_dynamic'
    return 'late/other_dynamic'
S3['class']=S3.apply(cls,axis=1)
S3.to_csv(f'{OUT}/stage3_survivors_temporal.csv', index=False)

print('\nClass distribution (all 111 survivor-rows):')
print(S3['class'].value_counts().to_string())
print('\nBy compartment x class:')
print(pd.crosstab(S3['compartment'], S3['class']).to_string())

print('\n=== EGF-like survivors (membrane + early 2-8min co-jump) ===')
egf=S3[S3['class']=='EGF-like (mem+early)'].sort_values('combined',ascending=False)
cols=['Kinase','Substrate','comp_r','motif_pct_max','best_site_aa','best_site_pos','sub_maxlog2fc','sub_dir']
print(egf[cols].round(3).to_string(index=False))

# also show Mem early_dynamic (mem survivors that jump early even if not strictly mem-max)
print('\n=== Other early-dynamic Mem survivors ===')
memearly=S3[(S3['compartment']=='Mem')&(S3['class']=='early_dynamic')].sort_values('combined',ascending=False)
print(memearly[cols].round(3).head(10).to_string(index=False))

# ---- figure: top Mem survivors temporal profiles ----
memsurv=S3[S3['compartment']=='Mem'].sort_values('combined',ascending=False).head(9)
fig,axes=plt.subplots(3,3,figsize=(14,11))
for ax,(_,r) in zip(axes.flat, memsurv.iterrows()):
    K,S=r['Kinase'],r['Substrate']
    for g,col,lw in [(K,'black',2.5),(S,'#e74c3c',1.8)]:
        f=log2fc(prof(g,'Mem'))
        if f is not None: ax.plot(TPS,f,marker='o',color=col,lw=lw,label=g)
    ax.axhline(0,color='#bbb',lw=0.6); ax.axvspan(2,8,color='#f9e79f',alpha=0.3)
    ax.set_title(f'{K} -> {S}  (Mem r={r["comp_r"]:.2f}, motif={r["motif_pct_max"]:.0f})',fontsize=9)
    ax.set_xlabel('min'); ax.set_ylabel('log2 FC vs 0min'); ax.legend(fontsize=8); ax.set_xticks(TPS)
fig.suptitle('Top 9 Membrane survivors — abundance trajectories (yellow = early 2-8 min window)',fontsize=12)
plt.tight_layout()
fig.savefig(f'{OUT}/stage3_mem_survivor_trajectories.png',dpi=300,bbox_inches='tight',facecolor='white')
save_tiff(fig,f'{OUT}/stage3_mem_survivor_trajectories.png'); plt.close(fig)
print('\nSaved stage3_mem_survivor_trajectories.png + .tiff and stage3_survivors_temporal.csv')
