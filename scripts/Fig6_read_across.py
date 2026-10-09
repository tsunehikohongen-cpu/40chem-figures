#!/usr/bin/env python3
"""Build Fig. 6 — Read-across of the 5 Unknown environmental chemicals against
Category 1 in-group and Category 0 out-group anchor sets.

Methods
-------
- In-group anchors (n=5): 5 Category 1 compounds with established developmental-
  toxicity mechanisms:
      Valproic acid, Tretinoin, Isotretinoin, Methotrexate, Thalidomide
- Out-group anchors (n=4): 4 Category 0 non-toxic reference compounds:
      L-Ascorbic acid, Saccharin Sodium Salt hydrate, D-Glucitol,
      Acetylsalicylic acid
- For each Unknown chemical: cosine similarity of its ssGSEA NES profile
  (2,692 GO BP terms) is computed against each anchor; mean similarity to
  the in-group vs out-group is reported in panel (a), and the full
  5 × 9 matrix in panel (b).

Input
-----
- data/NES_matrix_GO_BP_2023_40chemicals.csv : term × 40-chemical ssGSEA NES matrix
    (20 Cat1 | 15 Cat0 | 5 Unknown). Run from the repository root.

Output
------
- output/Fig6_read_across.png / .pdf (600 dpi)
"""
import os
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from fig_style import apply_style, DPI, CATEGORY_COLORS

apply_style()

NES_CSV = 'data/NES_matrix_GO_BP_2023_40chemicals.csv'

IN_GROUP = ['Valproic acid', 'Tretinoin', 'Isotretinoin',
            'Methotrexate', 'Thalidomide']
OUT_GROUP = ['L-Ascorbic acid', 'Saccharin Sodium Salt hydrate',
             'D-Glucitol', 'Acetylsalicylic acid']
UNKNOWN_ORDER = ['Deltamethrin', "2,2',4,4',5-Pentabromodiphenyl ether",
                 'Bisphenol A', 'Chlorpyrifos', 'imidacloprid']

IN_COLOR  = '#C0392B'   # red / Category 1
OUT_COLOR = '#1cb0b9'   # teal / Category 0 non-toxic references

# ---- Load NES matrix and compute full 40×40 similarity ----
nes = pd.read_csv(NES_CSV, index_col=0)
chems = list(nes.columns)
assert len(chems) == 40
sim = pd.DataFrame(cosine_similarity(nes.T.values), index=chems, columns=chems)

# ---- Build read-across tables ----
anchors = IN_GROUP + OUT_GROUP
unknowns = [u for u in UNKNOWN_ORDER if u in chems]
assert len(unknowns) == 5, f'missing Unknown chemicals: {set(UNKNOWN_ORDER)-set(chems)}'
mat = sim.loc[unknowns, anchors].astype(float)

means = pd.DataFrame({
    'in_mean':  sim.loc[unknowns, IN_GROUP].mean(axis=1),
    'out_mean': sim.loc[unknowns, OUT_GROUP].mean(axis=1),
})
means['delta'] = means['in_mean'] - means['out_mean']
print(means.round(3))

# ============================= figure =============================
fig = plt.figure(figsize=(18, 6))
gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 1.2], wspace=0.30,
                       left=0.08, right=0.97, top=0.90, bottom=0.14)

# ---- Panel a: grouped bars (in-group vs out-group per Unknown) ----
ax_a = fig.add_subplot(gs[0])
y = np.arange(len(unknowns))[::-1]
h = 0.38
ax_a.barh(y + h/2, means.loc[unknowns, 'in_mean'],  height=h,
          color=IN_COLOR,  edgecolor='black', lw=0.4,
          label='In-group anchors\n(established Category 1 mechanisms)')
ax_a.barh(y - h/2, means.loc[unknowns, 'out_mean'], height=h,
          color=OUT_COLOR, edgecolor='black', lw=0.4,
          label='Out-group anchors\n(non-toxic Category 0 references)')
ax_a.set_yticks(y); ax_a.set_yticklabels(unknowns, fontsize=12)
ax_a.set_xlabel('Mean cosine similarity\n(GSEA-NES profile)', fontsize=13)
ax_a.set_xlim(0, 0.9)
ax_a.legend(loc='lower center', bbox_to_anchor=(0.5, -0.28),
            ncol=2, frameon=False, fontsize=10)
ax_a.set_title('a  Read-across: Unknown chemicals vs. anchors',
               fontsize=14, loc='left', pad=10)

# ---- Panel b: 5 × 9 heatmap of individual Unknown × anchor similarity ----
gs_b = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1],
    width_ratios=[25, 0.7, 0.9], wspace=0.05)
ax_b   = fig.add_subplot(gs_b[0])
ax_grp = fig.add_subplot(gs_b[1])
cax    = fig.add_subplot(gs_b[2])

im = ax_b.imshow(mat.values, aspect='auto', cmap='RdBu_r',
                 vmin=0.40, vmax=0.90, interpolation='nearest')
# annotate values
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        v = mat.values[i, j]
        tc = 'white' if (v < 0.50 or v > 0.82) else 'black'
        ax_b.text(j, i, f'{v:.2f}', ha='center', va='center',
                  fontsize=10, color=tc)
ax_b.set_xticks(range(len(anchors)))
ax_b.set_xticklabels(anchors, rotation=45, ha='right', fontsize=11)
ax_b.set_yticks(range(len(unknowns)))
ax_b.set_yticklabels(unknowns, fontsize=11, fontweight='bold')
ax_b.set_title('b  Unknown chemicals vs. individual anchors',
               fontsize=14, loc='left', pad=10)

# anchor-group color strip
ax_grp.set_xlim(0, 1); ax_grp.set_ylim(0, len(unknowns))
for i in range(len(unknowns)):
    ax_grp.add_patch(plt.Rectangle((0, i), 0.5, 1,
        facecolor=IN_COLOR, edgecolor='white', lw=0.2))
    ax_grp.add_patch(plt.Rectangle((0.5, i), 0.5, 1,
        facecolor=OUT_COLOR, edgecolor='white', lw=0.2))
ax_grp.axis('off')
ax_grp.text(0.5, len(unknowns) + 0.3, 'anchor\ngroup',
            ha='center', va='bottom', fontsize=9, fontweight='bold')

cbar = fig.colorbar(im, cax=cax)
cbar.set_label('Cosine similarity', fontsize=11, fontweight='bold')
cbar.ax.tick_params(labelsize=10)

os.makedirs('output', exist_ok=True)
for ext in ('png', 'pdf'):
    fig.savefig(f'output/Fig6_read_across.{ext}', dpi=DPI, bbox_inches='tight',
                facecolor='white')
    print(f'saved output/Fig6_read_across.{ext}')
plt.close(fig)
