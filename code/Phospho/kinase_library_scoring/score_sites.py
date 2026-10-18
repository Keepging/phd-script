"""Kinase Library scoring of phosphosites.

Inputs
  site_meta.tsv          site_key, site_id (ACC_[STY]pos), protein, gene, ...
  uniprot_sequences.tsv  output of fetch_uniprot.py (canonical sequences)
Outputs (in out_dir)
  site_sequences.tsv         site_id, protein, position, residue_expected, residue_in_uniprot,
                             residue_check, window_-5_+4, substrate_15mer, uniprot_entry, sequence_version
  kinase_scores_long.parquet site_id, kinase, kinase_family, score, percentile
Scoring: kinase-library (PyPI) PhosphoProteomics.score / .percentile with package defaults
  (st_fav=True, pp=False, Ochoa phosphoproteome percentiles). S/T sites -> all ser_thr matrices,
  Y sites -> all tyrosine matrices incl. the dual-specificity *_TYR matrices (non_canonical=True).
"""
import argparse, csv, json, re, sys, warnings
warnings.filterwarnings('ignore')
import pandas as pd
import kinase_library as kl
from kinase_library.modules import data as kdata

ap = argparse.ArgumentParser()
ap.add_argument('site_meta'); ap.add_argument('uniprot_tsv'); ap.add_argument('out_dir')
a = ap.parse_args()

meta = pd.read_csv(a.site_meta, sep='\t', dtype=str)
sites = meta[['site_id', 'protein']].drop_duplicates('site_id').reset_index(drop=True)
m = sites['site_id'].str.extract(r'^(?P<acc>[A-Z0-9]+(?:-\d+)?)_(?P<res>[STY])(?P<pos>\d+)$')
assert m.notna().all().all(), 'unparseable site_id present'
sites['residue_expected'] = m['res']; sites['position'] = m['pos'].astype(int)
assert (m['acc'] == sites['protein']).all(), 'site_id accession differs from protein column'

uni = pd.read_csv(a.uniprot_tsv, sep='\t', dtype=str, keep_default_na=False)
uni = uni.set_index('requested_accession')
seqs = uni['sequence'].to_dict()

PAD = 7
def window(seq, pos, left, right):
    """residues pos-left .. pos+right (1-based pos), '_' outside the protein."""
    p = '_' * PAD + seq + '_' * PAD
    i = pos - 1 + PAD
    return p[i - left:i + right + 1]

rows = []
for r in sites.itertuples(index=False):
    seq = seqs.get(r.protein, '')
    e = uni.loc[r.protein] if r.protein in uni.index else None
    if not seq:
        check, actual, w10, w15 = 'accession_not_found', '', '', ''
    elif r.position > len(seq):
        check, actual, w10, w15 = 'position_out_of_range', '', '', ''
    else:
        actual = seq[r.position - 1]
        w10 = window(seq, r.position, 5, 4)
        w15 = window(seq, r.position, 7, 7)
        if actual == r.residue_expected:
            check = 'match'
        elif actual in 'STY':
            check = 'mismatch_STY'       # different phosphoacceptor at this position
        else:
            check = 'mismatch_non_STY'   # not a phosphoacceptor -> cannot be scored
    rows.append(dict(site_id=r.site_id, protein=r.protein, position=r.position,
                     residue_expected=r.residue_expected, residue_in_uniprot=actual, residue_check=check,
                     **{'window_-5_+4': w10, 'substrate_15mer': w15},
                     uniprot_entry=e['uniprot_entry'] if e is not None else '',
                     sequence_version=e['sequence_version'] if e is not None else '',
                     sequence_length=len(seq) if seq else ''))
ss = pd.DataFrame(rows)
ss.to_csv(f'{a.out_dir}/site_sequences.tsv', sep='\t', index=False)
print(ss['residue_check'].value_counts().to_string())

# ---- scoring (sites with S/T/Y at the position; residue taken from the UniProt sequence)
scorable = ss[ss['residue_check'].isin(['match', 'mismatch_STY'])].copy()
scorable['Sequence'] = scorable['substrate_15mer']
pp = kl.PhosphoProteomics(scorable[['site_id', 'Sequence']], seq_col='Sequence', suppress_warnings=True)
assert len(pp.omited_entries) == 0, pp.omited_entries

info = kdata.get_kinome_info().set_index('MATRIX_NAME')
parts = []
for kin_type, kw in (('ser_thr', {}), ('tyrosine', {'non_canonical': True})):
    kins = kdata.get_kinase_list(kin_type, **kw)
    sub = getattr(pp, kin_type + '_data')
    if len(sub) == 0:
        continue
    sc = pp.score(kin_type=kin_type, kinases=kins, values_only=True, **kw)
    pc = pp.percentile(kin_type=kin_type, kinases=kins, values_only=True, **kw)
    sc.index = sub['site_id'].values; pc.index = sub['site_id'].values
    long = (sc.stack().rename('score').to_frame()
              .join(pc.stack().rename('percentile')).reset_index()
              .rename(columns={'level_0': 'site_id', 'level_1': 'kinase'}))
    long['kinase_family'] = long['kinase'].map(info['FAMILY'])
    parts.append(long)
out = pd.concat(parts, ignore_index=True)[['site_id', 'kinase', 'kinase_family', 'score', 'percentile']]
out['score'] = out['score'].astype('float64'); out['percentile'] = out['percentile'].astype('float64')
out.to_parquet(f'{a.out_dir}/kinase_scores_long.parquet', index=False)
print(out.shape, out['kinase'].nunique(), 'kinases;', out['site_id'].nunique(), 'sites scored')
print(out.head())
