"""Shared core: sequence windows + kinase name mapping + batch motif scoring."""
import json, os
import pandas as pd
import kinase_library as kl

JSON_DIR = r'D:\博士\Protein contour\Phospho\biophys_json'
ST_ATLAS = set(kl.get_kinase_list(kin_type='ser_thr'))
TY_ATLAS = set(kl.get_kinase_list(kin_type='tyrosine'))

ALIAS = {
    'ABL1':'ABL','ABL2':'ARG','YES1':'YES','ERBB2':'HER2','MERTK':'MER','FLT4':'VEGFR3',
    'PTK2':'FAK','PTK2B':'PYK2','PTK6':'BRK','TNK2':'ACK','TNK1':'TNK1','FRK':'FRK','FER':'FER',
    'MAPK1':'ERK2','MAPK3':'ERK1','MAPK8':'JNK1','MAPK9':'JNK2','MAPK10':'JNK3',
    'MAPK14':'P38A','MAPK11':'P38B','MAPK12':'P38G','MAPK13':'P38D','MAPK7':'ERK5','MAPK15':'ERK7','MAPK6':'ERK3',
    'MAP2K1':'MEK1','MAP2K2':'MEK2','MAP2K5':'MEK5',
    'MAP3K1':'MEKK1','MAP3K2':'MEKK2','MAP3K3':'MEKK3','MAP3K5':'ASK1','MAP3K7':'TAK1',
    'MAP3K8':'COT','MAP3K9':'MLK1','MAP3K10':'MLK2','MAP3K11':'MLK3','MAP3K12':'DLK','MAP3K15':'MAP3K15',
    'RPS6KB1':'P70S6K','RPS6KB2':'P70S6KB','RPS6KA1':'P90RSK','RPS6KA2':'RSK3','RPS6KA3':'RSK2',
    'RPS6KA4':'MSK2','RPS6KA5':'MSK1','RPS6KA6':'RSK4',
    'CHEK1':'CHK1','CHEK2':'CHK2','MAPKAPK5':'MK5',
    'PRKAA1':'AMPKA1','PRKAA2':'AMPKA2','PRKACA':'PKACA','PRKACB':'PKACB','PRKACG':'PKACG',
    'PRKCA':'PKCA','PRKCB':'PKCB','PRKCD':'PKCD','PRKCE':'PKCE','PRKCG':'PKCG','PRKCH':'PKCH',
    'PRKCI':'PKCI','PRKCQ':'PKCT','PRKCZ':'PKCZ','PRKD1':'PKD1','PRKD2':'PKD2','PRKD3':'PKD3',
    'PRKG1':'PKG1','PRKG2':'PKG2','PRKDC':'DNAPK',
    'CSNK1A1':'CK1A','CSNK1D':'CK1D','CSNK1E':'CK1E','CSNK1G1':'CK1G1','CSNK1G2':'CK1G2','CSNK1G3':'CK1G3',
    'CSNK2A1':'CK2A1','CSNK2A2':'CK2A2',
    'AURKA':'AURA','AURKB':'AURB','AURKC':'AURC','STK11':'LKB1','MKNK1':'MNK1','MKNK2':'MNK2',
    'PDPK1':'PDK1','EEF2K':'EEF2K','PKMYT1':'MYT1',
    'STK38':'NDR1','STK38L':'NDR2','STK3':'MST2','STK4':'MST1','STK24':'MST3','STK26':'MST4',
    'IKBKB':'IKKB','IKBKE':'IKKE','CHUK':'IKKA','OXSR1':'OSR1',
    'SIK1':'SIK','SIK2':'SIK','SIK3':'SIK',
    'EIF2AK2':'PKR','EIF2AK1':'HRI','EIF2AK3':'PERK','EIF2AK4':'GCN2',
}

def map_to_atlas(name):
    if name in ST_ATLAS: return name, 'ser_thr', 'direct'
    if name in TY_ATLAS: return name, 'tyrosine', 'direct'
    if name in ALIAS:
        al = ALIAS[name]
        if al in ST_ATLAS: return al, 'ser_thr', 'alias'
        if al in TY_ATLAS: return al, 'tyrosine', 'alias'
    return None, None, 'UNMATCHED'

_seq_cache = {}
def get_seq(acc):
    if acc in _seq_cache: return _seq_cache[acc]
    p = os.path.join(JSON_DIR, f'{acc}.json')
    if not os.path.exists(p):
        _seq_cache[acc] = None; return None
    d = json.load(open(p))
    seq = ''.join(r['aa'] for r in sorted(d['residues'], key=lambda x: x['seqpos']))
    _seq_cache[acc] = seq
    return seq

def cut_window(acc, pos, half=7):
    """Centered (2*half+1)-mer around 1-based pos; '_' padding at termini."""
    seq = get_seq(acc)
    if seq is None: return None, 'NO_JSON'
    if pos < 1 or pos > len(seq): return None, f'POS_OOR'
    center = seq[pos-1]
    left = seq[max(0,pos-1-half):pos-1]
    right = seq[pos:pos+half]
    left = '_'*(half-len(left)) + left
    right = right + '_'*(half-len(right))
    return left+center+right, center
