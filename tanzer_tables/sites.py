"""
PXD021366 (Tanzer 2021) crude fractionation: collapse Spectronaut reports to site / protein tables.

Inputs (taken out of the PRIDE zips, not kept in the repo):
    phos: "20200224 crude fractionation phospho TNF_Report.txt"   (20200224crudefractionationphosphoTNF_Report.zip)
    prot: "20200303 TNF crude fractionation proteome _Report.txt"  (20200303_082317_200301TNFcrudefractionationproteomeagainstcombinedlibrary.zip)

Rules:
    site      = phospho residue with probability >= 0.75 in EG.PTMLocalizationProbabilities
    protein   = first Swiss-Prot entry in PG.UniProtIds (PG.ProteinNames not starting with the accession), else the first
    position  = PEP.PeptidePosition of that protein + offset in peptide - 1
    quantity  = sum of EG.TotalQuantity (Settings) over precursors carrying the site in a run (> 0 only)
    loc_prob  = max localization probability of the site over those precursors
    protein quantity: PG.Quantity > 0; "Filtered" / NaN -> missing

Usage:
    python sites.py phos_report.txt prot_report.txt [--out DIR]
"""
import argparse
import collections
import re

import numpy as np
import pandas as pd

CUT = 0.75
COMP = ['cyto', 'membrane', 'nucleus']
CONDS = ('DMSO', 'TNF', 'TNF_7Oxo')
fn_re = re.compile(r'U937_(cyto|membrane|nucleus|total)_phospho_(DMSO|TNF_7Oxo|TNF)_(\d)$')
rre = re.compile(r'U937_(cyto|membrane|nucleus|total)_(DMSO|TNF_7Oxo|TNF)_(\d)\.raw')
tok = re.compile(r'([A-Z])(\[[^\]]*\])?')


def sites_of(locstr):
    # returns list of (offset1based, residue, prob) with phospho prob >= CUT
    s = locstr.strip('_'); out = []; pos = 0
    for m in tok.finditer(s):
        pos += 1
        mod = m.group(2)
        if mod and mod.startswith('[Phospho'):
            p = float(mod.split(':')[1].strip(' %]')) / 100
            if p >= CUT: out.append((pos, m.group(1), p))
    return out


def rep_index(acc, names):
    for i, (a, n) in enumerate(zip(acc, names)):
        if not n.startswith(a.split('-')[0]): return i
    return 0


def read_sites(phos):
    cols = ['R.FileName', 'PG.UniProtIds', 'PG.ProteinNames', 'PG.Genes', 'PEP.PeptidePosition',
            'EG.PTMLocalizationProbabilities', 'EG.TotalQuantity (Settings)']
    quant = collections.defaultdict(float)   # (site, run) -> summed quantity
    prob = collections.defaultdict(float)    # (site, run) -> max localization probability
    nprec = collections.Counter()            # (site, run) -> contributing precursors
    siteinfo = {}
    for ch in pd.read_csv(phos, sep='\t', usecols=cols, chunksize=200000, dtype=str, keep_default_na=False):
        for fn, ids, names, genes, pp, loc, q in ch[cols].itertuples(index=False):
            if 'Phospho' not in loc: continue
            q = float(q)
            if not (q > 0): continue
            acc = ids.split(';'); nm = names.split(';'); gn = genes.split(';'); pos = pp.split(';')
            i = rep_index(acc, nm)
            if len(pos) != len(acc): i = 0
            start = int(pos[i].split(',')[0])
            run = fn_re.search(fn.strip()).groups()
            for off, aa, p in sites_of(loc):
                site = f'{acc[i]}_{aa}{start + off - 1}'
                quant[(site, run)] += q
                prob[(site, run)] = max(prob[(site, run)], p)
                nprec[(site, run)] += 1
                siteinfo[site] = (acc[i], gn[i] if i < len(gn) else '', aa, start + off - 1)
    rows = [(site, *siteinfo[site], f'{c}_{cond}_{rep}', c, cond, int(rep), quant[k], prob[k], nprec[k])
            for k in quant for site, (c, cond, rep) in [k]]
    df = pd.DataFrame(rows, columns=['site_key', 'accession', 'gene', 'residue', 'position', 'run',
                                     'compartment', 'condition', 'replicate', 'quantity', 'loc_prob', 'n_precursors'])
    for c in ('compartment', 'condition', 'run'): df[c] = df[c].astype('category')
    return df.sort_values(['site_key', 'compartment', 'condition', 'replicate'], ignore_index=True)


def read_proteins(prot):
    pt = pd.read_csv(prot, sep='\t', dtype=str)
    qcols = [c for c in pt.columns if c.endswith('.PG.Quantity')]
    names = ['_'.join(rre.search(c).groups()) for c in qcols]
    q = pt[qcols].replace({'Filtered': np.nan, 'NaN': np.nan}).astype(float)
    q = q.where(q > 0)
    q.columns = names
    meta = pt[['PG.ProteinGroups', 'PG.Genes', 'PG.UniProtIds']].set_axis(['protein_group', 'genes', 'uniprot_ids'], axis=1)
    return pd.concat([meta, q], axis=1)


def report(sites, prot):
    pres = sites.groupby('site_key').apply(lambda d: set(zip(d['compartment'], d['condition'], d['replicate'].astype(str))))
    siteprot = sites.drop_duplicates('site_key').set_index('site_key')['accession']
    qcols = prot.columns[3:]
    runs = [tuple(c.rsplit('_', 1)[0].split('_', 1)) + (c.rsplit('_', 1)[1],) for c in qcols]
    ppres = collections.defaultdict(set)
    for accs, vals in zip(prot['uniprot_ids'], prot[qcols].to_numpy()):
        got = {r for r, v in zip(runs, vals) if v == v}
        for a in accs.split(';'): ppres[a] |= got

    def any_run(rs, comp): return any(r[0] == comp for r in rs)
    def valid3(rs, comp):  # >=3 of 4 replicates in at least one condition
        return any(sum(1 for r in rs if r[0] == comp and r[1] == cond) >= 3 for cond in CONDS)
    def all12(rs, comp): return sum(1 for r in rs if r[0] == comp) == 12

    print('localized sites (>=0.75) quantified in >=1 run:', len(pres))
    for label, f in [('>=1 run', any_run), ('>=3/4 reps in >=1 condition', valid3), ('all 12 runs', all12)]:
        s3 = [s for s in pres.index if all(f(pres[s], c) for c in COMP)]
        s3p_any = [s for s in s3 if ppres.get(siteprot[s])]
        s3p3 = [s for s in s3 if all(f(ppres.get(siteprot[s], set()), c) for c in COMP)]
        print(f'[{label}] site in all 3 compartments: {len(s3)}; + parent protein quantified anywhere in proteome: '
              f'{len(s3p_any)}; + parent protein in all 3 compartments (same rule): {len(s3p3)}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('phos'); ap.add_argument('prot')
    ap.add_argument('--out', help='directory for tanzer_sites.parquet / tanzer_proteins.parquet')
    a = ap.parse_args()
    sites = read_sites(a.phos)
    prot = read_proteins(a.prot)
    report(sites, prot)
    if a.out:
        sites.to_parquet(f'{a.out}/tanzer_sites.parquet', index=False, compression='zstd')
        prot.to_parquet(f'{a.out}/tanzer_proteins.parquet', index=False, compression='zstd')
