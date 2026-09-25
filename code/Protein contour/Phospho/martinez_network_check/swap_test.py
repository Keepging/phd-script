import pandas as pd, numpy as np
from scipy.stats import pearsonr
import warnings; warnings.filterwarnings('ignore')

OUT = r'D:\博士\Protein contour\Phospho\martinez_network_check'
pam = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\protein_abundance_matrix.csv')
cond = [c for c in pam.columns if c!='Gene']
for c in cond: pam[c]=pd.to_numeric(pam[c],errors='coerce')
pam=pam.set_index('Gene')
corr=pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\kinase_substrate_correlation.csv')
corr=corr[corr['Kinase']!=corr['Substrate']]

def vec(g):
    return pam.loc[g, cond].values.astype(float) if g in pam.index else None
def r_pair(a,b):
    va,vb=vec(a),vec(b)
    if va is None or vb is None: return np.nan
    m=~(np.isnan(va)|np.isnan(vb))
    if m.sum()<3 or np.std(va[m])==0 or np.std(vb[m])==0: return np.nan
    return pearsonr(va[m],vb[m])[0]
def mean_r(kin, subs):
    rs=[r_pair(kin,s) for s in subs]; rs=[r for r in rs if not np.isnan(r)]
    return np.mean(rs) if rs else np.nan

# anchors (substrate-group owners) and the kinases swapped in (incl cross-family/type)
anchors = ['CDK2','MAPK1','AKT1']
swap_in = ['CDK2','MAPK1','AKT1','GSK3B','CSNK2A1','PRKACA','SRC']
families = {'CDK2':'CMGC/proline','MAPK1':'CMGC/proline','AKT1':'AGC/basophilic',
           'GSK3B':'CMGC/proline','CSNK2A1':'acidophilic','PRKACA':'AGC/basophilic','SRC':'Tyr'}

groups={}
for a in anchors:
    top5=corr[corr['Kinase']==a].nlargest(5,'Overall_Pearson_r')['Substrate'].tolist()
    groups[a]=top5
    print(f'{a} top-5 substrates: {top5}')

# matrix rows=substrate groups, cols=swapped kinase
mat=pd.DataFrame(index=[f'{a} top5' for a in anchors], columns=swap_in, dtype=float)
for a in anchors:
    for k in swap_in:
        mat.loc[f'{a} top5', k]=mean_r(k, groups[a])
print('\n=== SWAP MATRIX: mean r(kinase, substrate-group) over 15 conditions ===')
print(mat.round(3).to_string())

# quantify: original (diagonal) vs swapped
print('\n=== Specificity quantification ===')
deltas=[]
for a in anchors:
    row=mat.loc[f'{a} top5']
    orig=row[a]
    others=row.drop(a)
    drop=orig-others
    deltas.extend(drop.values)
    print(f'{a}: original r={orig:.3f} | swapped mean={others.mean():.3f} '
          f'(min {others.min():.3f}, max {others.max():.3f}) | mean drop={drop.mean():.3f}')
    # cross-family specifically
    cross=others[[k for k in others.index if families[k]!=families[a]]]
    print(f'    cross-family swapped mean r={cross.mean():.3f} (drop {orig-cross.mean():.3f})')

deltas=np.array(deltas)
print(f'\nOverall: mean r drop when swapping kinase = {deltas.mean():.3f} (median {np.median(deltas):.3f})')
mat.to_csv(f'{OUT}/swap_test_matrix.csv')
print('Saved swap_test_matrix.csv')
print('\nVERDICT:', 'LOW specificity - correlation barely distinguishes kinases (swap drop small)'
      if abs(deltas.mean())<0.15 else 'SOME specificity - original kinase scores notably higher than swapped')
