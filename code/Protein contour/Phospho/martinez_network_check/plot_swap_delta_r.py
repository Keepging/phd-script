import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns

# ── load ──────────────────────────────────────────────────────────────────────
p = r'D:\博士\Protein contour\Phospho\martinez_network_check\swap_test_matrix.csv'
mat = pd.read_csv(p, index_col=0)
print("Matrix shape:", mat.shape)
print(mat.round(4))

# Row labels are "CDK2 top5", "MAPK1 top5", "AKT1 top5"
# Column labels are kinase names; cognate = col whose name appears in row label
def cognate_col(row_name, cols):
    for c in cols:
        if c.upper() in row_name.upper():
            return c
    return None

deltas = []
for row_name, row in mat.iterrows():
    cog = cognate_col(row_name, mat.columns.tolist())
    if cog is None:
        print(f"  WARNING: no cognate found for {row_name!r}")
        continue
    r_cog = row[cog]
    for col in mat.columns:
        if col == cog:
            continue
        deltas.append(r_cog - row[col])

deltas = np.array(deltas)
mean_d  = deltas.mean()
med_d   = np.median(deltas)
print(f"\nTotal Δr values: {len(deltas)}")
print(f"Mean Δr = {mean_d:.4f},  Median Δr = {med_d:.4f}")

# ── figure ────────────────────────────────────────────────────────────────────
font = {'family': 'Arial', 'size': 11}
matplotlib.rc('font', **font)

fig, ax = plt.subplots(figsize=(8, 5))

bins = np.linspace(deltas.min() - 0.02, deltas.max() + 0.02, 26)

# histogram
ax.hist(deltas, bins=bins, color='steelblue', alpha=0.55,
        edgecolor='white', linewidth=0.6, density=True, label='Δr distribution')

# KDE overlay
kde = gaussian_kde(deltas, bw_method=0.35)
xkde = np.linspace(bins[0], bins[-1], 400)
ax.plot(xkde, kde(xkde), color='#1a3a6b', linewidth=2.2, label='KDE')

# vertical reference lines
ax.axvline(mean_d, color='#d62728', linestyle='--', linewidth=1.8, zorder=5)
ax.axvline(med_d,  color='#ff7f0e', linestyle='--', linewidth=1.8, zorder=5)

ymax = ax.get_ylim()[1]
ax.text(mean_d + 0.008, ymax * 0.88,
        f'mean = {mean_d:.2f}', color='#d62728',
        fontsize=11, style='italic', va='top')
ax.text(med_d  + 0.008, ymax * 0.74,
        f'median = {med_d:.2f}', color='#ff7f0e',
        fontsize=11, style='italic', va='top')

# annotation
ax.text(0.97, 0.97,
        'Swapping kinase identity\ndrops r by only ~0.09 on average',
        transform=ax.transAxes, fontsize=11, style='italic',
        ha='right', va='top', color='#444444',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#f5f5f5',
                  edgecolor='#cccccc', alpha=0.85))

ax.set_xlabel('Δr  (cognate − swapped kinase)', fontsize=13)
ax.set_ylabel('Density', fontsize=13)
ax.tick_params(labelsize=11)
ax.set_xlim(-0.12, 0.50)

ax.set_title('Swap Test: Kinase Specificity of Abundance Correlation',
             fontsize=16, fontweight='bold', pad=10)

legend = ax.legend(fontsize=11, frameon=False, loc='upper left')

sns.despine(ax=ax)
ax.set_facecolor('white'); fig.patch.set_facecolor('white')

plt.tight_layout()
out = r'D:\博士\Protein contour\Phospho\martinez_network_check\swap_test_delta_r_distribution.png'
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f"\nSaved: {out}")
