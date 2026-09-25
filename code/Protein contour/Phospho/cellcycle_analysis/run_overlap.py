"""
Cell-cycle overlap analysis: Martinez-Val (PXD023690) vs PXD011836
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import os

OUTDIR = r'D:\博士\Protein contour\Phospho\cellcycle_analysis'
os.makedirs(OUTDIR, exist_ok=True)

def save_tiff(fig, filename, dpi=300):
    path = filename if filename.endswith('.tiff') else filename.rsplit('.', 1)[0] + '.tiff'
    fig.savefig(path, dpi=dpi, bbox_inches='tight', format='tiff', pil_kwargs={'compression': 'tiff_lzw'})

# ============================================================
# Step 1: Martinez-Val genes
# ============================================================
pam = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\protein_abundance_matrix.csv')
martinez_genes = set(pam['Gene'].dropna().astype(str).str.strip())
print(f'Martinez-Val genes (protein_abundance_matrix): {len(martinez_genes)}')

# ============================================================
# Step 2: PXD011836 phospho genes
# ============================================================
phos = pd.read_csv(r'D:\博士\Protein contour\Phospho_HeLa_CellCycle\Phospho (STY)Sites.txt',
                    sep='\t', low_memory=False)
phos = phos[phos['Reverse'] != '+']
phos = phos[phos['Potential contaminant'] != '+']
cc_genes_phos = set()
for g in phos['Gene names'].dropna().astype(str):
    for part in g.split(';'):
        part = part.strip()
        if part:
            cc_genes_phos.add(part)
print(f'PXD011836 phospho genes (Phospho(STY)Sites.txt): {len(cc_genes_phos)}')

# Also proteinGroups
pg = pd.read_csv(r'D:\博士\Protein contour\Phospho_HeLa_CellCycle\proteinGroups.txt',
                  sep='\t', low_memory=False)
pg = pg[pg['Reverse'] != '+']
pg = pg[pg['Potential contaminant'] != '+']
cc_genes_pg = set()
for g in pg['Gene names'].dropna().astype(str):
    for part in g.split(';'):
        part = part.strip()
        if part:
            cc_genes_pg.add(part)
print(f'PXD011836 protein genes (proteinGroups.txt): {len(cc_genes_pg)}')

cellcycle_genes = cc_genes_phos

# ============================================================
# Step 3: Overlap
# ============================================================
intersection = martinez_genes & cellcycle_genes
martinez_only = martinez_genes - cellcycle_genes
cellcycle_only = cellcycle_genes - martinez_genes

print(f'\n=== OVERLAP STATISTICS ===')
print(f'Martinez-Val genes:      {len(martinez_genes)}')
print(f'PXD011836 phospho genes: {len(cellcycle_genes)}')
print(f'Intersection:            {len(intersection)} '
      f'({len(intersection)/len(martinez_genes)*100:.1f}% of MV, '
      f'{len(intersection)/len(cellcycle_genes)*100:.1f}% of CC)')
print(f'Martinez-only:           {len(martinez_only)}')
print(f'CellCycle-only:          {len(cellcycle_only)}')

# ============================================================
# Step 4: CDK top 20 substrates
# ============================================================
ks = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\kinase_substrate_correlation.csv')
ns = ks[ks['Kinase'] != ks['Substrate']].copy()
ns = ns.dropna(subset=['Overall_Pearson_r'])

cdk_pairs = ns[ns['Kinase'].apply(lambda x: isinstance(x, str) and x.upper().startswith('CDK'))].copy()
cdk_top20 = cdk_pairs.nlargest(20, 'Overall_Pearson_r')

print(f'\n=== CDK TOP 20 SUBSTRATES vs CELL-CYCLE DATA ===')
results = []
for _, row in cdk_top20.iterrows():
    sub = row['Substrate']
    in_cc = sub in cellcycle_genes
    results.append({
        'Kinase': row['Kinase'],
        'Substrate': sub,
        'Sites': row['Sites'],
        'Overall_Pearson_r': round(row['Overall_Pearson_r'], 4),
        'In_CellCycle_Data': in_cc
    })
    marker = 'YES' if in_cc else ' NO'
    print(f'  [{marker}] {row["Kinase"]:8s} -> {sub:12s}  '
          f'r={row["Overall_Pearson_r"]:.4f}  sites={row["Sites"]}')

n_found = sum(1 for r in results if r['In_CellCycle_Data'])
print(f'\nCDK top 20 in cell-cycle data: {n_found}/20 ({n_found/20*100:.0f}%)')

res_df = pd.DataFrame(results)
res_df.to_csv(os.path.join(OUTDIR, 'cdk_top20_cellcycle_overlap.csv'), index=False)

# ============================================================
# Step 5: Venn diagram
# ============================================================
fig, ax = plt.subplots(figsize=(8, 6))
v = venn2([martinez_genes, cellcycle_genes],
          set_labels=('Martinez-Val\n(PXD023690)\nTotal proteome',
                      'PXD011836\nCell-cycle phospho'),
          set_colors=('#3498db', '#e74c3c'),
          alpha=0.6, ax=ax)

for text in v.set_labels:
    if text:
        text.set_fontsize(12)
for text in v.subset_labels:
    if text:
        text.set_fontsize(14)

ax.set_title(
    f'Gene Overlap: Martinez-Val Proteome vs PXD011836 Cell-Cycle Phosphoproteome\n'
    f'Intersection: {len(intersection)} genes '
    f'({len(intersection)/len(cellcycle_genes)*100:.1f}% of cell-cycle phospho genes covered)',
    fontsize=11, pad=15)

plt.tight_layout()
venn_path = os.path.join(OUTDIR, 'cellcycle_overlap_venn.png')
fig.savefig(venn_path, dpi=300, bbox_inches='tight', facecolor='white')
save_tiff(fig, venn_path)
plt.close(fig)
print(f'\nSaved Venn: {venn_path} + .tiff')

# ============================================================
# Step 6: CDK substrate phase intensities from PerseusOutput
# ============================================================
print(f'\n=== STEP 6: CDK substrate phase intensities ===')

perseus = pd.read_csv(
    r'D:\博士\Protein contour\Phospho_HeLa_CellCycle\PerseusOutput Phospho(STY).txt',
    sep='\t', skiprows=[1, 2])
print(f'PerseusOutput loaded: {perseus.shape}')

phase_cols = ['G1_1', 'G1_2', 'G1_3', 'S_1', 'S_2', 'S_3',
              'G2_1', 'G2_2', 'G2_3', 'Mix']
for c in phase_cols:
    perseus[c] = pd.to_numeric(perseus[c], errors='coerce')

cdk_subs_in_cc = [r['Substrate'] for r in results if r['In_CellCycle_Data']]
print(f'CDK substrates to extract: {len(cdk_subs_in_cc)}')

perseus['gene_list'] = (perseus['Gene names'].fillna('').astype(str)
                        .apply(lambda x: [g.strip() for g in x.split(';') if g.strip()]))

out_rows = []
for gene in cdk_subs_in_cc:
    mask = perseus['gene_list'].apply(lambda gl: gene in gl)
    sub_df = perseus[mask]
    if len(sub_df) == 0:
        print(f'  WARNING: {gene} not found in PerseusOutput')
        continue
    for _, row in sub_df.iterrows():
        g1_mean = np.nanmean([row['G1_1'], row['G1_2'], row['G1_3']])
        s_mean = np.nanmean([row['S_1'], row['S_2'], row['S_3']])
        g2_mean = np.nanmean([row['G2_1'], row['G2_2'], row['G2_3']])
        out_rows.append({
            'Gene': gene,
            'UniqueID': row.get('Unique identifier', ''),
            'Position': row.get('Position', ''),
            'AminoAcid': row.get('Amino acid', ''),
            'Localization_prob': row.get('Localization prob', ''),
            'G1_1': row['G1_1'], 'G1_2': row['G1_2'], 'G1_3': row['G1_3'],
            'S_1': row['S_1'], 'S_2': row['S_2'], 'S_3': row['S_3'],
            'G2_1': row['G2_1'], 'G2_2': row['G2_2'], 'G2_3': row['G2_3'],
            'Mix': row['Mix'],
            'G1_mean': round(g1_mean, 4) if not np.isnan(g1_mean) else np.nan,
            'S_mean': round(s_mean, 4) if not np.isnan(s_mean) else np.nan,
            'G2_mean': round(g2_mean, 4) if not np.isnan(g2_mean) else np.nan,
            'S_vs_G1': round(s_mean - g1_mean, 4) if not (np.isnan(s_mean) or np.isnan(g1_mean)) else np.nan,
            'G2_vs_G1': round(g2_mean - g1_mean, 4) if not (np.isnan(g2_mean) or np.isnan(g1_mean)) else np.nan,
        })

phase_df = pd.DataFrame(out_rows)
phase_path = os.path.join(OUTDIR, 'cellcycle_cdk_substrates_by_phase.csv')
phase_df.to_csv(phase_path, index=False)
print(f'Saved: {phase_path}')
print(f'  Rows: {len(phase_df)} phosphosites across {phase_df["Gene"].nunique()} genes')

print(f'\nPer-gene summary (mean log2 intensity by phase):')
gene_summary = phase_df.groupby('Gene')[['G1_mean', 'S_mean', 'G2_mean', 'S_vs_G1', 'G2_vs_G1']].mean()
for gene, row in gene_summary.iterrows():
    print(f'  {gene:12s}  G1={row["G1_mean"]:+.3f}  S={row["S_mean"]:+.3f}  '
          f'G2={row["G2_mean"]:+.3f}  S-G1={row["S_vs_G1"]:+.3f}  G2-G1={row["G2_vs_G1"]:+.3f}')

print(f'\nAll outputs in: {OUTDIR}')
