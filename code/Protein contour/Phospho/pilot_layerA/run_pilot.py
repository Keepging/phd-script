# Layer A pilot: HPA compartment x Scop3P phosphosites x b2bTools features
# READ-ONLY w.r.t. project files. All outputs stay in pilot_layerA/.
import csv, json, os, time, urllib.request, urllib.error
import numpy as np
from collections import Counter, defaultdict
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

ROOT = r'D:\博士\Protein contour\Phospho'
PILOT = os.path.join(ROOT, 'pilot_layerA')
BIO = os.path.join(ROOT, 'biophys_json')
CACHE = os.path.join(PILOT, 'scop3p_pilot_cache.json')
CAP_PER_COMP = 350          # bound API load for the pilot
FEATS = ['backbone','sidechain','disoMine','helix','sheet','coil','earlyFolding']
# Layer-B-aligned 3 compartments via clean single-compartment HPA main location
COMP_TOKENS = {
 'Nucleus': {'Nucleoplasm','Nucleoli','Nucleoli fibrillar center','Nucleoli rim',
             'Nuclear bodies','Nuclear speckles','Nuclear membrane','Kinetochore','Mitotic chromosome'},
 'Cytosol': {'Cytosol'},
 'Membrane': {'Plasma membrane','Endoplasmic reticulum'},
}

# ---- gene -> uniprot ----
g2u = {}
with open(os.path.join(ROOT,'gene_uniprot_map.csv'),encoding='utf-8') as f:
    for r in csv.DictReader(f):
        if r['gene'] and r['uniprot'] and r['gene'] not in g2u:
            g2u[r['gene']] = r['uniprot']

def has_bio(acc): return os.path.exists(os.path.join(BIO, acc+'.json'))

# ---- HPA clean assignment ----
assign = {}  # uniprot -> comp
counts_stage = Counter()
with open(os.path.join(PILOT,'subcellular_location.tsv'),encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['Reliability']=='Uncertain': continue
        main = [t.strip() for t in r['Main location'].split(';') if t.strip()]
        if not main: continue
        comp = None
        for c,toks in COMP_TOKENS.items():
            if set(main) <= toks:   # ALL main tokens within one compartment => clean
                comp = c; break
        if comp is None: continue
        counts_stage['hpa_clean']+=1
        gn = r['Gene name']
        acc = g2u.get(gn)
        if not acc: counts_stage['no_uniprot']+=1; continue
        if not has_bio(acc): counts_stage['no_biophys']+=1; continue
        assign.setdefault(acc, comp)
print('HPA clean single-compartment proteins:', counts_stage['hpa_clean'])
print('  dropped no-uniprot:',counts_stage['no_uniprot'],' no-biophys:',counts_stage['no_biophys'])
bycomp = defaultdict(list)
for acc,c in assign.items(): bycomp[c].append(acc)
for c in COMP_TOKENS: print(f'  mapped+biophys {c}: {len(bycomp[c])} proteins')

# deterministic cap
sel = {}
for c in COMP_TOKENS:
    accs = sorted(bycomp[c])[:CAP_PER_COMP]
    for a in accs: sel[a]=c
print('selected proteins total (capped):', len(sel))

# ---- Scop3P live pull (cached) ----
cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
def scop3p(acc):
    if acc in cache: return cache[acc]
    url='https://iomics.ugent.be/scop3p/api/modifications?accession='+acc
    try:
        req=urllib.request.Request(url,headers={'Accept':'application/json'})
        d=json.loads(urllib.request.urlopen(req,timeout=25).read().decode())
        mods=d.get('modifications') or []
        out=[{'position':m.get('position'),'residue':m.get('residue'),'name':m.get('name')}
             for m in mods if m.get('name') and 'phospho' in str(m.get('name')).lower()]
    except Exception as e:
        out={'__error__':repr(e)}
    cache[acc]=out
    return out

n=0; errs=0
for acc in sel:
    r=scop3p(acc); n+=1
    if isinstance(r,dict) and '__error__' in r: errs+=1
    if n%50==0:
        json.dump(cache,open(CACHE,'w'))
        print(f'  scop3p pulled {n}/{len(sel)} (errors {errs})', flush=True)
    time.sleep(0.12)
json.dump(cache,open(CACHE,'w'))
print(f'scop3p done: {n} proteins, {errs} errors')

# ---- map sites -> features ----
_bcache={}
def residues(acc):
    if acc in _bcache: return _bcache[acc]
    d=json.load(open(os.path.join(BIO,acc+'.json')))
    res={int(r['seqpos']):r for r in d['residues']}
    _bcache[acc]=res; return res

rows=[]
for acc,comp in sel.items():
    sites=cache.get(acc)
    if not isinstance(sites,list): continue
    res=residues(acc)
    seen=set()
    for s in sites:
        pos=s.get('position')
        if pos is None or pos in seen: continue
        seen.add(pos)
        rr=res.get(int(pos))
        if not rr: continue
        if rr.get('aa') not in ('S','T','Y'): continue   # phospho-acceptor sanity
        row={'uniprot':acc,'comp':comp,'pos':pos,'aa':rr['aa']}
        ok=all(rr.get(f) is not None for f in FEATS)
        if not ok: continue
        for f in FEATS: row[f]=float(rr[f])
        rows.append(row)

print('\n=== SITE COUNTS PER COMPARTMENT (mapped to features) ===')
cc=Counter(r['comp'] for r in rows)
for c in COMP_TOKENS: print(f'  {c}: {cc[c]} sites')
print('  total sites:', len(rows))
nprot=Counter();
for r in rows: nprot[(r['comp'],r['uniprot'])]=1
pcount=Counter(c for (c,_) in nprot)
for c in COMP_TOKENS: print(f'  {c}: {pcount[c]} proteins contributing sites')

# save site table
with open(os.path.join(PILOT,'pilot_sites.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['uniprot','comp','pos','aa']+FEATS); w.writeheader()
    for r in rows: w.writerow(r)

# ---- KW + BH + eps2 ----
groups=list(COMP_TOKENS.keys())
stat_rows=[]
for feat in FEATS:
    samples=[[r[feat] for r in rows if r['comp']==g] for g in groups]
    samples=[s for s in samples if len(s)>0]
    H,p=kruskal(*samples)
    k=len(samples); ntot=sum(len(s) for s in samples)
    eps2=(H-k+1)/(ntot-k) if ntot>k else float('nan')
    meds={g:np.median([r[feat] for r in rows if r['comp']==g]) for g in groups}
    stat_rows.append({'feature':feat,'H':H,'p':p,'eps2':eps2,'n':ntot,
                      **{f'med_{g}':meds[g] for g in groups}})
pvals=[s['p'] for s in stat_rows]
_,pbh,_,_=multipletests(pvals,method='fdr_bh')
for s,pb in zip(stat_rows,pbh): s['p_BH']=pb; s['sig']=pb<0.05

print('\n=== KRUSKAL-WALLIS per feature (3 compartments, BH-corrected) ===')
hdr=f"{'feature':16s} {'H':>9s} {'p':>11s} {'p_BH':>11s} {'sig':>5s} {'eps2':>9s}  medians(Nuc/Cyt/Mem)"
print(hdr)
for s in sorted(stat_rows,key=lambda x:-x['eps2']):
    print(f"{s['feature']:16s} {s['H']:9.3f} {s['p']:11.2e} {s['p_BH']:11.2e} {str(s['sig']):>5s} "
          f"{s['eps2']:9.5f}  {s['med_Nucleus']:.3f}/{s['med_Cytosol']:.3f}/{s['med_Membrane']:.3f}")
maxe=max(s['eps2'] for s in stat_rows)
print(f"\nHEADLINE max eps2 = {maxe:.5f}  (Layer B presence-pattern max eps2 = 0.01377)")

with open(os.path.join(PILOT,'pilot_stats.csv'),'w',newline='') as f:
    cols=['feature','H','p','p_BH','sig','eps2','n']+[f'med_{g}' for g in groups]
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader()
    for s in stat_rows: w.writerow({k:s[k] for k in cols})
print('\nSaved pilot_sites.csv, pilot_stats.csv, scop3p_pilot_cache.json in pilot_layerA/')
