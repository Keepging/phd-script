"""
4 sub-tasks on the kinase-substrate correlation table.
Input: kinase_correlation/kinase_substrate_correlation.csv
Output dir: multi_kinase_analysis/
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ============================================================
# Setup
# ============================================================
SRC = r'D:\博士\Protein contour\Phospho\kinase_correlation\kinase_substrate_correlation.csv'
OUTDIR = r'D:\博士\Protein contour\Phospho\multi_kinase_analysis'
os.makedirs(OUTDIR, exist_ok=True)

plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 150


def save_tiff(fig, filename, dpi=300):
    """User-provided helper. Saves as .tiff with LZW compression."""
    path = filename if filename.endswith('.tiff') else filename.rsplit('.', 1)[0] + '.tiff'
    fig.savefig(path, dpi=dpi, bbox_inches='tight',
                format='tiff', pil_kwargs={'compression': 'tiff_lzw'})


# ============================================================
# Load + filter non-self
# ============================================================
df = pd.read_csv(SRC)
print(f'Loaded: {df.shape}')
print(f'Columns: {list(df.columns)}')

ns = df[df['Kinase'] != df['Substrate']].copy()
ns = ns.dropna(subset=['Overall_Pearson_r'])
N_NONSELF = len(ns)
print(f'Non-self pairs (with valid Overall_Pearson_r): {N_NONSELF}')


def is_cdk(name):
    return isinstance(name, str) and name.upper().startswith('CDK')


ns['is_cdk'] = ns['Kinase'].apply(is_cdk)

# ============================================================
# TASK 1 — Top 50 pairs horizontal barh
# ============================================================
print('\n' + '=' * 60)
print('TASK 1: Top 50 pairs barh')
print('=' * 60)

top50 = ns.nlargest(50, 'Overall_Pearson_r').reset_index(drop=True)
top50['label'] = top50['Kinase'] + '  →  ' + top50['Substrate']
n_cdk_top50 = int(top50['is_cdk'].sum())
print(f'CDK in top 50: {n_cdk_top50}/50 ({n_cdk_top50/50*100:.0f}%)')

fig, ax = plt.subplots(figsize=(10, 14))
# Plot bottom-up so top of list ends up at the top of the axes
top50_plot = top50.iloc[::-1].reset_index(drop=True)
colors = ['#d62728' if cdk else '#1f77b4' for cdk in top50_plot['is_cdk']]
ax.barh(top50_plot['label'], top50_plot['Overall_Pearson_r'],
        color=colors, edgecolor='white', linewidth=0.4)
ax.set_xlabel('Overall Pearson r (15 conditions, non-self)', fontsize=11)
ax.set_title(f'Top 50 Kinase → Substrate Pairs by Overall Pearson r\n'
             f'CDK family: {n_cdk_top50}/50 ({n_cdk_top50/50*100:.0f}%)',
             fontsize=12)
ax.set_xlim(0.97, 1.005)  # r values are all close to 1 in the top 50
ax.tick_params(axis='y', labelsize=7.5)
ax.tick_params(axis='x', labelsize=9)
ax.grid(True, axis='x', alpha=0.3)

# Legend
import matplotlib.patches as mpatches
ax.legend(handles=[
    mpatches.Patch(color='#d62728', label=f'CDK family ({n_cdk_top50}/50)'),
    mpatches.Patch(color='#1f77b4', label=f'Other ({50 - n_cdk_top50}/50)'),
], loc='lower right', fontsize=10, frameon=True, edgecolor='#333333')

plt.tight_layout()
png1 = os.path.join(OUTDIR, 'task1_top50_barh.png')
fig.savefig(png1, dpi=300, bbox_inches='tight', facecolor='white')
save_tiff(fig, png1)
plt.close(fig)
print(f'Saved: {png1} + .tiff')

# ============================================================
# TASK 2 — Ranked lists
# ============================================================
print('\n' + '=' * 60)
print('TASK 2: Per-kinase / per-substrate top 5 lists')
print('=' * 60)

# Per-kinase top 5 substrates
def top5_rows(grp_df, group_col, target_col):
    out = []
    for g, sub in grp_df.groupby(group_col):
        s = sub.dropna(subset=['Overall_Pearson_r']).nlargest(5, 'Overall_Pearson_r')
        for rank, (_, row) in enumerate(s.iterrows(), start=1):
            out.append({
                'group': g,
                'Rank': rank,
                'target': row[target_col],
                'r': row['Overall_Pearson_r'],
                'Sites': row.get('Sites', ''),
            })
    return pd.DataFrame(out)

kin_top5 = top5_rows(ns, 'Kinase', 'Substrate')
sub_top5 = top5_rows(ns, 'Substrate', 'Kinase')

print(f'Kinase top-5 rows: {len(kin_top5)} (across {ns["Kinase"].nunique()} kinases)')
print(f'Substrate top-5 rows: {len(sub_top5)} (across {ns["Substrate"].nunique()} substrates)')

f_kin = os.path.join(OUTDIR, 'task2_kinase_top5_substrates.csv')
f_sub = os.path.join(OUTDIR, 'task2_substrate_top5_kinases.csv')
kin_top5.to_csv(f_kin, index=False)
sub_top5.to_csv(f_sub, index=False)
print(f'Saved: {f_kin}')
print(f'Saved: {f_sub}')

# ============================================================
# TASK 3 — Multi-kinase signal scan
# ============================================================
print('\n' + '=' * 60)
print('TASK 3: Substrates with multi-kinase clusters (r>0.7, within 0.05 of top)')
print('=' * 60)

t3_rows = []
for sub, grp in ns.groupby('Substrate'):
    high = grp[grp['Overall_Pearson_r'] > 0.7].copy()
    if len(high) < 2:
        continue
    top_r = high['Overall_Pearson_r'].max()
    cluster = high[high['Overall_Pearson_r'] >= top_r - 0.05]
    if len(cluster) < 2:
        continue
    cluster_sorted = cluster.sort_values('Overall_Pearson_r', ascending=False)
    kinases = cluster_sorted['Kinase'].tolist()
    rs = cluster_sorted['Overall_Pearson_r'].tolist()
    spread = max(rs) - min(rs)
    t3_rows.append({
        'Substrate': sub,
        'N_competing_kinases': len(cluster),
        'Kinases': ';'.join(kinases),
        'Correlations': ';'.join(f'{r:.4f}' for r in rs),
        'Spread': round(spread, 4),
        'Top_r': round(top_r, 4),
    })

t3 = pd.DataFrame(t3_rows).sort_values('N_competing_kinases', ascending=False)
f_t3 = os.path.join(OUTDIR, 'task3_multi_kinase_substrates.csv')
t3.to_csv(f_t3, index=False)
print(f'Substrates with >=2 competing kinases (r>0.7, top±0.05): {len(t3)}')
if len(t3):
    print(f'Distribution of N_competing_kinases:')
    print(t3['N_competing_kinases'].value_counts().sort_index(ascending=False).to_string())
print(f'\nTop 10 substrates by N_competing_kinases:')
print(t3.head(10)[['Substrate', 'N_competing_kinases', 'Kinases', 'Spread', 'Top_r']].to_string(index=False))
print(f'\nSaved: {f_t3}')

# ============================================================
# TASK 4 — Per-kinase heatmaps (top 20 kinases by pair count)
# ============================================================
print('\n' + '=' * 60)
print('TASK 4: Per-kinase heatmaps (top 20 kinases by # substrates)')
print('=' * 60)

# Top 20 kinases by number of substrate pairs (rows in ns)
pair_counts = ns.groupby('Kinase').size().sort_values(ascending=False)
top20_kinases = pair_counts.head(20).index.tolist()
print(f'Top 20 kinases by pair count:')
for k, c in zip(top20_kinases, pair_counts.head(20).values):
    print(f'  {k:12s}: {c} pairs')

ncols = ['Cyt_r', 'Mem_r', 'Nuc_r', 'Overall_Pearson_r']
xlabs = ['Cytoplasm', 'Membrane', 'Nucleus', 'Overall']

for kin in top20_kinases:
    sub = ns[ns['Kinase'] == kin].dropna(subset=['Overall_Pearson_r'])
    sub = sub.nlargest(30, 'Overall_Pearson_r').copy()
    if len(sub) == 0:
        continue
    sub['row_label'] = sub['Substrate'] + ' (' + sub['Sites'].fillna('').astype(str) + ')'
    mat = sub.set_index('row_label')[ncols].copy()
    mat.columns = xlabs

    fig, ax = plt.subplots(figsize=(7, max(6, 0.35 * len(mat) + 2)))
    sns.heatmap(mat, cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                annot=True, fmt='.2f', annot_kws={'size': 7},
                cbar_kws={'label': 'Pearson r'},
                linewidths=0.4, linecolor='white',
                mask=mat.isna(), ax=ax)
    ax.set_title(f'{kin} — top {len(mat)} substrates × 4 conditions', fontsize=12)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='x', labelsize=10)
    plt.tight_layout()
    safe = kin.replace('/', '_').replace(' ', '_')
    p = os.path.join(OUTDIR, f'task4_{safe}_heatmap.png')
    fig.savefig(p, dpi=300, bbox_inches='tight', facecolor='white')
    save_tiff(fig, p)
    plt.close(fig)

print(f'Saved 20 heatmap pairs (PNG + TIFF) to {OUTDIR}')

# ============================================================
# SUMMARY
# ============================================================
r_gt9 = ns[ns['Overall_Pearson_r'] > 0.9]
n_r_gt9 = len(r_gt9)
n_r_gt9_cdk = int(r_gt9['is_cdk'].sum())
cdk_pct_in_high = n_r_gt9_cdk / n_r_gt9 * 100 if n_r_gt9 else 0

med_all = ns['Overall_Pearson_r'].median()
med_cdk = ns[ns['is_cdk']]['Overall_Pearson_r'].median()
med_non = ns[~ns['is_cdk']]['Overall_Pearson_r'].median()

print('\n' + '=' * 60)
print('FINAL SUMMARY')
print('=' * 60)
print(f'Total non-self pairs (valid Overall_r): {N_NONSELF}')
print(f'r > 0.9 pairs: {n_r_gt9}  '
      f'(CDK: {n_r_gt9_cdk} = {cdk_pct_in_high:.1f}% of high-r pairs)')
print(f'Median Overall_r — all non-self:  {med_all:.4f}')
print(f'Median Overall_r — CDK kinases:   {med_cdk:.4f}  '
      f'(n={int(ns["is_cdk"].sum())})')
print(f'Median Overall_r — non-CDK:       {med_non:.4f}  '
      f'(n={int((~ns["is_cdk"]).sum())})')
print(f'\nAll outputs in: {OUTDIR}')
