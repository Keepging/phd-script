import pandas as pd, numpy as np, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings; warnings.filterwarnings('ignore')

base=r'D:\博士\Protein contour\Phospho\kinase_correlation'
OUT=r'D:\博士\Protein contour\Phospho\martinez_network_check'
TPS=[0,2,8,20,90]
COMP_DEF=[('Cytosol',['Cyt_%dmin'%t for t in TPS],'Cyt_r'),
          ('Membrane',['Mem_%dmin'%t for t in TPS],'Mem_r'),
          ('Nucleus',['Nuc_%dmin'%t for t in TPS],'Nuc_r'),
          ('Overall',['%s_%dmin'%(c,t) for c in ['Cyt','Mem','Nuc'] for t in TPS],'Overall_Pearson_r')]

prot=pd.read_csv(os.path.join(base,'protein_abundance_matrix.csv'))
for c in prot.columns:
    if c!='Gene': prot[c]=pd.to_numeric(prot[c],errors='coerce')
prot=prot.drop_duplicates('Gene').set_index('Gene')
corr=pd.read_csv(os.path.join(base,'kinase_substrate_correlation.csv'))

kinases=[k for k in corr['Kinase'].unique() if k in prot.index]
substrate_set=set(corr['Substrate'].unique())
nonsub=[g for g in prot.index if g not in substrate_set]   # non-substrate partner pool
print(f'kinases in matrix: {len(kinases)}, non-substrate partner pool: {len(nonsub)}')

def vpearson(A,B):
    """row-wise pearson for two arrays (n_pairs x k); NaN-safe per row needing >=3 valid shared points."""
    out=np.full(A.shape[0],np.nan)
    for i in range(A.shape[0]):
        a,b=A[i],B[i]; m=~(np.isnan(a)|np.isnan(b))
        if m.sum()>=3:
            aa,bb=a[m],b[m]
            if aa.std()>0 and bb.std()>0:
                out[i]=np.corrcoef(aa,bb)[0,1]
    return out

rng=np.random.default_rng(7)
N_SHUF=100
observed={}; randomized={}
kin_arr_cache={}
for name,cols,rcol in COMP_DEF:
    obs=corr[rcol].dropna().values
    observed[name]=obs
    n=len(obs)
    Kmat=prot.loc[kinases,cols].values
    Pmat=prot.loc[nonsub,cols].values
    rand_vals=[]
    for _ in range(N_SHUF):
        ki=rng.integers(0,len(kinases),size=n)
        pi=rng.integers(0,len(nonsub),size=n)
        rand_vals.append(vpearson(Kmat[ki],Pmat[pi]))
    randomized[name]=np.concatenate(rand_vals)
    randomized[name]=randomized[name][~np.isnan(randomized[name])]
    print(f'{name}: observed n={n} (med {np.median(obs):.2f}) | random n={len(randomized[name])} (med {np.median(randomized[name]):.2f})')

# ---------------- figure ----------------
plt.rcParams.update({'font.size':13})
fig,ax=plt.subplots(figsize=(12,7))
C_OBS='#1f6fb4'; C_RND='#d98c3f'
labels=[d[0] for d in COMP_DEF]
xc=np.arange(len(labels)); off=0.21; w=0.34
def draw(data,center,color):
    vp=ax.violinplot(data,positions=[center],widths=w,showextrema=False)
    for b in vp['bodies']:
        b.set_facecolor(color); b.set_alpha(0.75); b.set_edgecolor('white'); b.set_linewidth(0.8)
    med=np.median(data)
    ax.hlines(med,center-w/2,center+w/2,color='black',lw=2.2,zorder=5)
    return med
for i,name in enumerate(labels):
    draw(observed[name],xc[i]-off,C_OBS)
    draw(randomized[name],xc[i]+off,C_RND)

ax.set_xticks(xc); ax.set_xticklabels(labels,fontsize=14)
ax.tick_params(axis='y',labelsize=13); ax.tick_params(axis='x',labelsize=14)
ax.set_xlabel('Compartment',fontsize=16)
ax.set_ylabel('Abundance correlation (Pearson r)',fontsize=16)
ax.set_title('Kinase–substrate abundance correlation vs randomized pairs, by compartment',fontsize=17,pad=14)
ax.axhline(0,color='#999',lw=1,ls='--',zorder=0)
ax.set_ylim(-1.05,1.15)
ax.legend(handles=[Patch(facecolor=C_OBS,label='observed'),
                   Patch(facecolor=C_RND,label='randomized baseline'),
                   plt.Line2D([0],[0],color='black',lw=2.2,label='median')],
          fontsize=13,loc='lower left',frameon=True,edgecolor='#cccccc')
ax.grid(axis='y',alpha=0.25)
plt.tight_layout()
p=os.path.join(OUT,'kinase_substrate_correlation_vs_random.png')
fig.savefig(p,dpi=300,bbox_inches='tight',facecolor='white')
plt.close(fig)
print('Saved',p)
