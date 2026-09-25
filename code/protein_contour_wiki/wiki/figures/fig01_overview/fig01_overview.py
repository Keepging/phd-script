"""
Fig 01 — Project Overview Data Flow
Layered DAG: raw → methods → analyses → findings, with MARCKS sub-branch.
20 nodes total. Core novelty (K-S correlation axis) highlighted.

Outputs: fig01_overview.png (300 dpi) + fig01_overview.svg
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
import os

# ============================================================
# Style constants
# ============================================================
COLORS = {
    'data':                '#4a90d9',  # blue
    'method':              '#9e9e9e',  # gray
    'analysis':            '#f39c12',  # orange
    'finding_validated':   '#27ae60',  # green
    'finding_negative':    '#c0392b',  # red
    'finding_preliminary': '#f1c40f',  # yellow
    'special':             '#1abc9c',  # teal
}
HIGHLIGHT_EDGE = '#d35400'  # bold orange border for core-novelty
EDGE_STYLES = {
    'main': dict(color='#444444', lw=1.4, linestyle='-',  arrowstyle='-|>'),
    'core': dict(color='#d35400', lw=2.6, linestyle='-',  arrowstyle='-|>'),
    'fail': dict(color='#c0392b', lw=1.4, linestyle='--', arrowstyle='-|>'),
}

# ============================================================
# Node definitions
# Each: id -> dict(label, type, x, y, w, h, highlight=False)
# Coord system: x ∈ [0, 14], y ∈ [-2, 9]; positive y is up.
# ============================================================
nodes = {
    # --- Raw data (col 0) ---
    'raw': dict(
        label='Martinez-Val 2021\nPXD023690\n5 tp × 6 fr × 4 rep\n~6,952 prot / 7,957 sites',
        type='data', x=1.0, y=4.5, w=2.4, h=1.6, highlight=False),

    # --- Methods / preprocessing (col 1) ---
    'seqfrac': dict(
        label='Sequential\nfractionation\n6 fr → 3 comp\n(Cyt | Mem | Nuc)',
        type='method', x=4.5, y=7.4, w=2.2, h=1.4, highlight=False),
    'b2b': dict(
        label='b2bTools\n7 per-residue\nfeatures',
        type='method', x=4.5, y=5.7, w=2.2, h=1.2, highlight=False),
    'mob': dict(
        label='Movement Score\n≥10 % = High\n(587 / 8,470)',
        type='method', x=4.5, y=4.0, w=2.2, h=1.3, highlight=False),
    'kinace': dict(
        label='KinAce K-S DB\n(PSP + EPSD)\n~16,400 entries',
        type='data', x=4.5, y=2.3, w=2.2, h=1.3, highlight=False),

    # --- Analyses (col 2) ---
    'a_phos': dict(
        label='Phospho vs\nnon-phospho',
        type='analysis', x=8.0, y=7.4, w=1.9, h=1.0, highlight=False),
    'a_mob': dict(
        label='High vs Low\nmobility biophys',
        type='analysis', x=8.0, y=6.1, w=1.9, h=1.0, highlight=False),
    'a_hub': dict(
        label='Hub vs Low-only\n(length-matched\n1000 iter)',
        type='analysis', x=8.0, y=4.8, w=1.9, h=1.1, highlight=False),
    'a_ks': dict(
        label='K-S Pearson r\n× 15 conditions\n(+ null baseline)',
        type='analysis', x=8.0, y=3.0, w=2.0, h=1.3, highlight=True),  # ★
    'a_bandle': dict(
        label='BANDLE\nattempt',
        type='analysis', x=8.0, y=1.3, w=1.7, h=0.9, highlight=False),

    # --- Findings (col 3) ---
    'F004': dict(
        label='F004  [V]\nPhospho proteins\nLOWER backbone-dyn,\n+disoMine\n(7/7 p<0.001)',
        type='finding_validated', x=11.6, y=7.4, w=2.4, h=1.5, highlight=False),
    'F001': dict(
        label='F001  [V]\nHigh-mob sites\ndistinct biophys\n5/6 sig  |d|≤0.33',
        type='finding_validated', x=11.6, y=5.8, w=2.4, h=1.4, highlight=False),
    'F002': dict(
        label='F002  [V]\nHub > Low-only\n6/7 features\n(length-matched)',
        type='finding_validated', x=11.6, y=4.3, w=2.4, h=1.4, highlight=False),
    'F003': dict(
        label='F003  [V]  ★\nCORE NOVELTY\nMem K-S r = 0.614\nvs random 0.216\np = 1.09e-77',
        type='finding_validated', x=11.6, y=2.5, w=2.5, h=1.7, highlight=True),  # ★
    'F005': dict(
        label='F005  [X]\nBANDLE failed:\n6-fr insufficient',
        type='finding_negative', x=11.6, y=0.8, w=2.3, h=1.1, highlight=False),
    'F006': dict(
        label='F006  [X]\nBH between-group:\nall 84 tests NS',
        type='finding_negative', x=11.6, y=-0.4, w=2.3, h=1.1, highlight=False),

    # --- MARCKS sub-branch (single y-row at y=-2.2) ---
    'marcks': dict(
        label='MARCKS\ncase study',
        type='special', x=1.0, y=-2.2, w=2.0, h=1.0, highlight=False),
    'm_l12': dict(
        label='L1 + L2:\nbiophys profile\n+ K-S coupling',
        type='analysis', x=4.5, y=-2.2, w=2.2, h=1.2, highlight=False),
    'm_l3': dict(
        label='L3 structural\n(AF / DSSP /\nSASA / NLS-NES)\nplanned',
        type='finding_preliminary', x=8.0, y=-2.2, w=2.2, h=1.3, highlight=False),
    'm_l4': dict(
        label='L4 occupancy\nBLOCKED',
        type='finding_negative', x=11.6, y=-2.2, w=2.0, h=1.0, highlight=False),
}

# ============================================================
# Edges
# (src, dst, style_key)
# ============================================================
edges = [
    # raw → methods
    ('raw', 'seqfrac', 'main'),
    ('raw', 'b2b',     'main'),
    ('raw', 'mob',     'main'),
    ('raw', 'kinace',  'main'),

    # methods → analyses
    ('b2b',     'a_phos', 'main'),
    ('b2b',     'a_mob',  'main'),
    ('b2b',     'a_hub',  'main'),
    ('mob',     'a_mob',  'main'),
    ('mob',     'a_hub',  'main'),
    ('seqfrac', 'a_ks',   'core'),     # core path
    ('kinace',  'a_ks',   'core'),     # core path
    ('seqfrac', 'a_bandle','main'),

    # analyses → findings
    ('a_phos',   'F004', 'main'),
    ('a_mob',    'F001', 'main'),
    ('a_hub',    'F002', 'main'),
    ('a_ks',     'F003', 'core'),      # core path
    ('a_bandle', 'F005', 'fail'),
    ('a_hub',    'F006', 'fail'),

    # MARCKS sub-branch
    ('raw',    'marcks', 'main'),
    ('marcks', 'm_l12',  'main'),
    ('m_l12',  'm_l3',   'main'),
    ('m_l3',   'm_l4',   'main'),
]

# ============================================================
# Drawing
# ============================================================
fig, ax = plt.subplots(figsize=(19, 12))
ax.set_xlim(-0.5, 17.0)
ax.set_ylim(-3.4, 10.5)
ax.set_aspect('equal')
ax.axis('off')

# ----- Background panels for visual grouping -----
def column_band(x_center, y_top, y_bot, label, color='#f7f7f7'):
    """Faint gray vertical band behind a column of nodes."""
    rect = mpatches.Rectangle((x_center - 1.3, y_bot), 2.6, y_top - y_bot,
                              facecolor=color, edgecolor='none', zorder=0, alpha=0.55)
    ax.add_patch(rect)
    ax.text(x_center, y_top + 0.15, label, ha='center', va='bottom',
            fontsize=11, fontweight='bold', color='#666666', zorder=1)

# Main flow column bands (y ~ -1.0 to 9.0)
column_band(1.0,  9.0, -1.05, 'INPUT')
column_band(4.5,  9.0, -1.05, 'METHODS / PREPROCESSING')
column_band(8.0,  9.0, -1.05, 'ANALYSES')
column_band(11.6, 9.0, -1.05, 'FINDINGS')

# MARCKS subgraph background band (y < -1.5)
marcks_band = mpatches.FancyBboxPatch(
    (0.1, -2.95), 13.5, 1.55,
    boxstyle="round,pad=0.08",
    facecolor='#f0fbf8', edgecolor='#1abc9c', linewidth=1.5,
    linestyle='--', zorder=0, alpha=0.7)
ax.add_patch(marcks_band)
ax.text(0.4, -1.45, 'MARCKS 4-layer case study  (side branch)',
        ha='left', va='top',
        fontsize=11, fontweight='bold', color='#0e6f5c', style='italic',
        zorder=1)

# ----- Helper: get node anchor on a given side -----
def anchor(node_id, side):
    n = nodes[node_id]
    cx, cy, w, h = n['x'], n['y'], n['w'], n['h']
    if side == 'right':  return (cx + w/2, cy)
    if side == 'left':   return (cx - w/2, cy)
    if side == 'top':    return (cx, cy + h/2)
    if side == 'bottom': return (cx, cy - h/2)
    raise ValueError(side)

# ----- Decide best anchor sides for each edge -----
def edge_anchors(src, dst):
    sx = nodes[src]['x']; dx = nodes[dst]['x']
    sy = nodes[src]['y']; dy = nodes[dst]['y']
    if dx > sx + 0.3:
        return anchor(src, 'right'), anchor(dst, 'left')
    if dx < sx - 0.3:
        return anchor(src, 'left'),  anchor(dst, 'right')
    if dy > sy:
        return anchor(src, 'top'),    anchor(dst, 'bottom')
    return     anchor(src, 'bottom'), anchor(dst, 'top')

# ----- Draw edges first (so nodes overlay) -----
for src, dst, style_key in edges:
    a, b = edge_anchors(src, dst)
    st = EDGE_STYLES[style_key]
    arrow = FancyArrowPatch(
        a, b,
        arrowstyle=st['arrowstyle'],
        mutation_scale=14,
        color=st['color'],
        lw=st['lw'],
        linestyle=st['linestyle'],
        connectionstyle='arc3,rad=0.05',
        zorder=2,
        shrinkA=0, shrinkB=2)
    ax.add_patch(arrow)

# ----- Draw nodes -----
def draw_node(nid, n):
    cx, cy, w, h = n['x'], n['y'], n['w'], n['h']
    face = COLORS[n['type']]
    edge_c = HIGHLIGHT_EDGE if n['highlight'] else '#222222'
    edge_lw = 2.8 if n['highlight'] else 0.9
    box = FancyBboxPatch(
        (cx - w/2, cy - h/2), w, h,
        boxstyle='round,pad=0.04,rounding_size=0.18',
        facecolor=face,
        edgecolor=edge_c,
        linewidth=edge_lw,
        zorder=3,
        alpha=0.92)
    ax.add_patch(box)
    # text colour: white on dark, black on yellow
    text_color = '#1a1a1a' if n['type'] in ('finding_preliminary',) else 'white'
    ax.text(cx, cy, n['label'], ha='center', va='center',
            fontsize=8.5, fontweight='bold' if n['highlight'] else 'normal',
            color=text_color, zorder=4, linespacing=1.15)

for nid, n in nodes.items():
    draw_node(nid, n)

# ============================================================
# Hypothesis overlay (right margin annotation)
# ============================================================
hyp_x = 14.5
hyp_panel = FancyBboxPatch(
    (hyp_x - 0.1, 1.6), 2.4, 5.6,
    boxstyle="round,pad=0.08", facecolor='#f5eef8',
    edgecolor='#9b59b6', linewidth=1.2, zorder=1, alpha=0.85)
ax.add_patch(hyp_panel)
ax.text(hyp_x + 1.1, 6.9, 'HYPOTHESES', ha='center', va='top',
        fontsize=10, fontweight='bold', color='#5b2c6f', zorder=2)
hyp_text = (
    'H1  [X] falsified\n'
    'site-level code\n(FWO original)\n'
    '|d|<=0.33 too small\n\n'
    'H2  [~] testing\n'
    'hardware-software\n'
    '3/3 confirmed:\n'
    'F002, F003, F004\n\n'
    'H3  [?] untested\n'
    'localization-\ndomain exposure\n'
    '(Lennart)\nneeds MARCKS L3'
)
ax.text(hyp_x + 1.1, 6.5, hyp_text, ha='center', va='top',
        fontsize=8.0, color='#5b2c6f', zorder=2, linespacing=1.25)

# ============================================================
# Legend
# ============================================================
legend_items = [
    ('Raw / curated data',           COLORS['data']),
    ('Preprocessing / method',       COLORS['method']),
    ('Analysis (active)',            COLORS['analysis']),
    ('Validated finding',            COLORS['finding_validated']),
    ('Preliminary / planned',        COLORS['finding_preliminary']),
    ('Negative / blocked',           COLORS['finding_negative']),
    ('Case-study branch (MARCKS)',   COLORS['special']),
]
legend_handles = [mpatches.Patch(facecolor=c, edgecolor='#222222',
                                  linewidth=0.7, label=lbl)
                  for lbl, c in legend_items]
# add highlight + edge type entries
legend_handles += [
    mpatches.Patch(facecolor='white', edgecolor=HIGHLIGHT_EDGE, linewidth=2.4,
                   label='★ Core novelty (highlighted)'),
    Line2D([0], [0], color='#444444', lw=1.4, label='Main data flow'),
    Line2D([0], [0], color='#d35400', lw=2.6, label='Core-novelty axis'),
    Line2D([0], [0], color='#c0392b', lw=1.4, linestyle='--', label='Failed / blocked'),
]
leg = ax.legend(handles=legend_handles, loc='lower right',
                 bbox_to_anchor=(1.00, 0.005), ncol=1, fontsize=8.5,
                 frameon=True, edgecolor='#222222', framealpha=0.95,
                 title='Node / edge encoding', title_fontsize=9.5)
leg.get_title().set_fontweight('bold')

# ============================================================
# Title
# ============================================================
ax.text(8.0, 10.20,
        'Project Overview — Data Flow from Martinez-Val Raw Data to All Validated Findings',
        ha='center', va='bottom', fontsize=15, fontweight='bold')
ax.text(8.0, 9.85,
        'Fig 01  |  20 nodes  |  Core novelty: kinase-substrate spatiotemporal correlation (a_ks -> F003)',
        ha='center', va='bottom', fontsize=10, color='#555555', style='italic')

plt.tight_layout()

# ============================================================
# Save
# ============================================================
outdir = os.path.dirname(os.path.abspath(__file__))
png_path = os.path.join(outdir, 'fig01_overview.png')
svg_path = os.path.join(outdir, 'fig01_overview.svg')

plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white')
print(f"OK Saved {png_path}")
print(f"OK Saved {svg_path}")
print(f"   Nodes: {len(nodes)}, Edges: {len(edges)}")
