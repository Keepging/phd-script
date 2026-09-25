import json, os
import pandas as pd

JSON_DIR = r'D:\博士\Protein contour\Phospho\biophys_json'

_seq_cache = {}
def get_seq(acc):
    if acc in _seq_cache:
        return _seq_cache[acc]
    p = os.path.join(JSON_DIR, f'{acc}.json')
    if not os.path.exists(p):
        _seq_cache[acc] = None
        return None
    d = json.load(open(p))
    seq = ''.join(r['aa'] for r in sorted(d['residues'], key=lambda x: x['seqpos']))
    _seq_cache[acc] = seq
    return seq

def cut_window(acc, pos, up=5, down=4):
    """Cut -up/+down window around 1-based pos. Pad with '_' at termini (atlas truncation char)."""
    seq = get_seq(acc)
    if seq is None:
        return None, 'NO_JSON'
    if pos < 1 or pos > len(seq):
        return None, f'POS_OOR(len={len(seq)})'
    center = seq[pos-1]
    left = seq[max(0, pos-1-up):pos-1]
    right = seq[pos:pos+down]
    # pad
    left = '_' * (up - len(left)) + left
    right = right + '_' * (down - len(right))
    window = left + center + right
    return window, center

ms = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_results\v2\kinase_matched_sites_v2.csv')
print(f'matched_sites_v2: {len(ms)} rows')
print(f'unique uniprot_id: {ms["uniprot_id"].nunique()}, null uniprot_id: {ms["uniprot_id"].isna().sum()}')

# Test on a few CDK substrate sites + a couple random
test_genes = [('RRM1', 559), ('FEN1', 187), ('MAP4', 696), ('ABCF1', 109)]
print('\n=== Window samples (up=5, down=4; center = phospho residue) ===')
for gene, pos in test_genes:
    rows = ms[(ms['Gene'] == gene) & (ms['Position'] == pos)]
    if len(rows) == 0:
        print(f'  {gene} pos{pos}: not found in matched_sites')
        continue
    row = rows.iloc[0]
    acc = row['uniprot_id']
    aa_expected = row['AA']
    win, center = cut_window(acc, pos)
    marker = ' '*5 + '^' + ' '*4  # points at center (index 5)
    print(f'  {gene} {aa_expected}{pos}  acc={acc}')
    print(f'    window: {win}   (expected center AA={aa_expected}, got={center})')
    print(f'            {marker}')
