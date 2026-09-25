import kinase_library as kl

print('=== kinase_library self-test ===\n')

# ST kinase list
st_kins = kl.get_kinase_list(kin_type='ser_thr')
print(f'ser_thr kinases available: {len(st_kins)}')
cdks = [k for k in st_kins if str(k).upper().startswith('CDK')]
print(f'CDK kinases in atlas: {cdks}\n')

# ---- Test 1: classic CDK proline-directed site  S-P-x-K ----
# Johnson atlas window is -5..+4 (10 residues). Center = phosphoacceptor.
# Build a 15-mer with S at center (phos_pos=8, 1-indexed).
# Sequence: ....  P at +1, K at +3  => canonical CDK consensus
seq_cdk = 'AAALLSPVKRAAAAA'   # S at index 6 (1-based), +1=P, +3=K
sub_cdk = kl.Substrate(seq_cdk, phos_pos=6)
print(f'Test1 CDK site seq={seq_cdk}  phosres={sub_cdk.substrate[5]}  (expect S, +1 P, +3 K)')
sc = sub_cdk.score(output_type='series')
top10 = sc.head(10)
print('Top 10 kinases by log2 score:')
print(top10.to_string())
for cdk in cdks:
    if cdk in sc.index:
        rank = list(sc.index).index(cdk) + 1
        print(f'  {cdk}: log2score={sc[cdk]:.3f}  rank={rank}/{len(sc)}')
print()

# ---- Test 2: basophilic non-proline site (AKT-like R-x-R-x-x-S) ----
seq_akt = 'AARQRTRSFSAAAAA'   # AKT consensus RxRxxS, no proline at +1
sub_akt = kl.Substrate(seq_akt, phos_pos=8)
print(f'Test2 AKT-like site seq={seq_akt}  phosres={sub_akt.substrate[7]}')
sc2 = sub_akt.score(output_type='series')
print('Top 10 kinases by log2 score:')
print(sc2.head(10).to_string())
for cdk in cdks:
    if cdk in sc2.index:
        rank = list(sc2.index).index(cdk) + 1
        print(f'  {cdk}: log2score={sc2[cdk]:.3f}  rank={rank}/{len(sc2)}')
print()

# ---- Try percentile (needs background phosphoproteome DB) ----
print('=== percentile test (needs background DB) ===')
try:
    pct = sub_cdk.percentile(kinases=cdks)
    print('CDK percentile for CDK site:')
    print(pct.to_string())
except Exception as e:
    print(f'percentile() failed: {type(e).__name__}: {e}')
