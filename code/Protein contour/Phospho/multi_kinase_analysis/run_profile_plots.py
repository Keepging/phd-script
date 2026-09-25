"""
Top 20 kinases × top 5 substrates: z-scored 15-condition abundance profile overlay.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

OUTDIR = r'D:\博士\Protein contour\Phospho\multi_kinase_analysis'

def save_tiff(fig, filename, dpi=300):
    path = filename if filename.endswith('.tiff') else filename.rsplit('.', 1)[0] + '.tiff'
    fig.savefig(path, dpi=dpi, bbox_inches='tight', format='tiff',
                pil_kwargs={'compression': 'tiff_lzw'})

# ============================================================
# Load
# ============================================================
pam = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\protein_abundance_matrix.csv')
cond_cols = [c for c in pam.columns if c != 'Gene']
for c in cond_cols:
    pam[c] = pd.to_numeric(pam[c], errors='coerce')
pam_idx = pam.set_index('Gene')
print(f'protein_abundance_matrix: {pam.shape}')

ks = pd.read_csv(r'D:\博士\Protein contour\Phospho\kinase_correlation\kinase_substrate_correlation.csv')
ns = ks[ks['Kinase'] != ks['Substrate']].copy()
ns = ns.dropna(subset=['Overall_Pearson_r'])
print(f'Non-self pairs: {len(ns)}')

# Ordered condition columns + short labels
ordered_cols = [f'{comp}_{t}min' for comp in ['Cyt', 'Mem', 'Nuc'] for t in [0, 2, 8, 20, 90]]
short_labels = [f'{c[0]}{t}' for c in ['C', 'M', 'N'] for t in [0, 2, 8, 20, 90]]

# Top 20 kinases by pair count
pair_counts = ns.groupby('Kinase').size().sort_values(ascending=False)
top20_kinases = pair_counts.head(20).index.tolist()

# Color palette for substrates
sub_colors = ['#e74c3c', '#2ecc71', '#3498db', '#9b59b6', '#f39c12']

for kin in top20_kinases:
    print(f'\n{"="*50}')
    print(f'{kin} ({pair_counts[kin]} pairs)')
    print(f'{"="*50}')

    # Top 5 substrates
    kin_pairs = ns[ns['Kinase'] == kin].nlargest(5, 'Overall_Pearson_r')

    # Check kinase in abundance matrix
    kin_found = kin in pam_idx.index
    if kin_found:
        kin_vals = pam_idx.loc[kin, ordered_cols].values.astype(float)
        kin_mean = np.nanmean(kin_vals)
        kin_std = np.nanstd(kin_vals, ddof=1)
        kin_z = (kin_vals - kin_mean) / kin_std if kin_std > 0 else kin_vals * 0
    else:
        kin_z = None
    print(f'  Kinase {kin}: {"FOUND" if kin_found else "NOT FOUND"}')

    fig, ax = plt.subplots(figsize=(12, 5))

    # Plot kinase
    x = np.arange(15)
    if kin_z is not None:
        ax.plot(x, kin_z, color='black', linewidth=3, marker='o', markersize=4,
                label=f'{kin} (kinase)', zorder=10)

    # Plot substrates
    missing_subs = []
    for i, (_, row) in enumerate(kin_pairs.iterrows()):
        sub = row['Substrate']
        r_val = row['Overall_Pearson_r']
        found = sub in pam_idx.index
        status = 'FOUND' if found else 'NOT FOUND'
        print(f'  #{i+1} {sub:12s}  r={r_val:.4f}  {status}')

        if found:
            sub_vals = pam_idx.loc[sub, ordered_cols].values.astype(float)
            sub_mean = np.nanmean(sub_vals)
            sub_std = np.nanstd(sub_vals, ddof=1)
            sub_z = (sub_vals - sub_mean) / sub_std if sub_std > 0 else sub_vals * 0
            ax.plot(x, sub_z, color=sub_colors[i], linewidth=1.5,
                    marker='s', markersize=3, alpha=0.85,
                    label=f'{sub} (r={r_val:.3f})')
        else:
            missing_subs.append(sub)

    # Compartment dividers
    ax.axvline(4.5, color='#7f8c8d', ls='--', lw=1, alpha=0.6)
    ax.axvline(9.5, color='#7f8c8d', ls='--', lw=1, alpha=0.6)

    # Compartment labels at top
    ax.text(2, ax.get_ylim()[1] if ax.get_ylim()[1] != 0 else 2, 'Cytoplasm',
            ha='center', va='bottom', fontsize=9, color='#7f8c8d', style='italic')
    ax.text(7, ax.get_ylim()[1] if ax.get_ylim()[1] != 0 else 2, 'Membrane',
            ha='center', va='bottom', fontsize=9, color='#7f8c8d', style='italic')
    ax.text(12, ax.get_ylim()[1] if ax.get_ylim()[1] != 0 else 2, 'Nucleus',
            ha='center', va='bottom', fontsize=9, color='#7f8c8d', style='italic')

    ax.set_xticks(x)
    ax.set_xticklabels(short_labels, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Z-scored abundance', fontsize=11)
    ax.set_xlabel('Condition (Compartment × Timepoint)', fontsize=10)
    ax.set_title(f'{kin} — Top 5 Substrates by Abundance Co-variation', fontsize=13)

    # Add missing note
    if not kin_found:
        ax.text(0.5, 0.95, f'⚠ Kinase {kin} not found in abundance matrix',
                transform=ax.transAxes, ha='center', fontsize=10, color='red',
                bbox=dict(boxstyle='round', facecolor='#ffe0e0', alpha=0.8))
    if missing_subs:
        note = 'Not in matrix: ' + ', '.join(missing_subs)
        ax.text(0.02, 0.02, note, transform=ax.transAxes, fontsize=7,
                color='#999', va='bottom')

    ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.01, 1),
              borderaxespad=0, frameon=True, edgecolor='#ccc')
    ax.grid(True, axis='y', alpha=0.2)

    plt.tight_layout()
    safe = kin.replace('/', '_').replace(' ', '_')
    p = os.path.join(OUTDIR, f'task4b_{safe}_profile.png')
    fig.savefig(p, dpi=300, bbox_inches='tight', facecolor='white')
    save_tiff(fig, p)
    plt.close(fig)

print(f'\n\nDone. Saved 20 profile plot pairs (PNG + TIFF) to {OUTDIR}')
