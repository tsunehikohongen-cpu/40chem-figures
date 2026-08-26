#!/usr/bin/env python3
# =====================================================================
# Figure 3a  --  Expression heatmap of the 182-gene set across chemicals
# Reproduction script (matches the published panel a).
#
# REQUIRED INPUT FILES (default: ./data/, override with --input-dir):
#   - Top1000_data.csv                              (rlog/log2FC expression, top-1000 genes)
#   - Wilcoxon_GeneSet182.csv                       (Gene column = the 182 genes)
#   - 1000_ECv_wilcoxon_nodes_LabelA_complete.tsv   (per-gene UP/DOWN direction)
# OUTPUT: Fig3A_Heatmap_bigfont.{png,tif,svg} (600 dpi) in --output-dir (default: ./output)
# Dependencies: python>=3.10, pandas, numpy, scipy, matplotlib, (font: Helvetica/Nimbus Sans/Arial)
# =====================================================================
import argparse
import os
import re

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Nimbus Sans', 'TeX Gyre Heros', 'Arial', 'DejaVu Sans'],
    'svg.fonttype': 'path',
})
COL_C1 = '#C0392B'
COL_C0 = '#2166AC'

# fixed chemical order & display names taken from the original figure
CAT1 = [('Acitretin', 'Acitretin'), ('Busulfan', 'Busulfan'), ('Carbamazepine', 'Carbamazepine'),
        ('Cyclophosphamide', 'Cyclophosphamide'), ('Fluconazole', 'Fluconazole'), ('X5Fluorouracil', '5-Fluorouracil'),
        ('Hydroxyurea', 'Hydroxyurea'), ('Isotretinoin', 'Isotretinoin'), ('Methotrexate', 'Methotrexate'),
        ('Phenytoin', 'Phenytoin'), ('Thalidomide', 'Thalidomide'), ('Topiramate', 'Topiramate'),
        ('Tretinoin', 'Tretinoin'), ('Valproicacid', 'Valproic acid'), ('Pomalidomide', 'Pomalidomide'),
        ('Aspirin', 'Aspirin'), ('Cytarabine', 'Cytarabine'), ('Ibuprofen', 'Ibuprofen'),
        ('Trimethadione', 'Trimethadione'), ('Vismodegib', 'Vismodegib')]
CAT0 = [('Ibrutinib', 'Ibrutinib'), ('Pazopanib', 'Pazopanib'), ('Ribavirin', 'Ribavirin'),
        ('Tacrolimus', 'Tacrolimus'), ('Bosentan', 'Bosentan'), ('Cisplatin', 'Cisplatin'),
        ('Dabrafenib', 'Dabrafenib'), ('Dasatinib', 'Dasatinib'), ('Imatinib', 'Imatinib'),
        ('Rotenone', 'Rotenone'), ('Acetaminophen', 'Acetaminophen'), ('Acetylsalicylicacid', '5-Acetylsalicylic acid'),
        ('DGlucitol', 'D-Glucitol'), ('LAscorbicacid', 'L-Ascorbic acid'), ('SaccharinSodiumSalthydrate', 'Saccharin Na')]
ORDER = CAT1 + CAT0
CHEMS = [k for k, _ in ORDER]
DISP = [v for _, v in ORDER]
NCAT1 = len(CAT1)


def chem(c):
    return re.sub(r'[_]?\d+$', '', c)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", default="data")
    p.add_argument("--output-dir", default="output")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    inp = args.input_dir

    d = pd.read_csv(os.path.join(inp, "Top1000_data.csv"))
    d['sym'] = d['symbol'].astype(str).str.strip()
    d = d.drop_duplicates('sym')
    numeric_cols = [c for c in d.columns[4:] if c != 'sym']
    d[numeric_cols] = d[numeric_cols].apply(pd.to_numeric, errors='coerce')
    genes_all = [g.strip() for g in pd.read_csv(os.path.join(inp, "Wilcoxon_GeneSet182.csv"))['Gene'].astype(str)]
    nodes = pd.read_csv(os.path.join(inp, "1000_ECv_wilcoxon_nodes_LabelA_complete.tsv"), sep="\t").set_index('name')

    expr_cols = [c for c in list(d.columns)[4:] if c != 'sym']
    cmap_cols = {}
    for c in expr_cols:
        cmap_cols.setdefault(chem(c), []).append(c)

    # mean across concentrations -> chemical x gene matrix
    mat = pd.DataFrame({k: d[cmap_cols[k]].mean(axis=1).values for k in CHEMS}, index=d['sym'])
    genes = [g for g in genes_all if g in mat.index]
    dirn = {g: (nodes.loc[g, 'direction'] if g in nodes.index else 'NS') for g in genes}
    up = [g for g in genes if dirn[g] == 'UP']
    dn = [g for g in genes if dirn[g] == 'DOWN']
    other = [g for g in genes if dirn[g] not in ('UP', 'DOWN')]
    gorder = up + dn + other
    nU = len(up)
    nD = len(dn)
    M = mat.loc[gorder, CHEMS].T  # rows=chemicals(ordered), cols=genes
    Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1).replace(0, 1), axis=0)

    # ---- figure (taller / more vertical) ----
    fig = plt.figure(figsize=(15, 14))
    gs = gridspec.GridSpec(2, 3, width_ratios=[0.30, 20, 0.45], height_ratios=[0.45, 20],
                            wspace=0.015, hspace=0.015)
    axL = fig.add_subplot(gs[1, 0])
    axT = fig.add_subplot(gs[0, 1])
    ax = fig.add_subplot(gs[1, 1])
    cax = fig.add_subplot(gs[1, 2])

    im = ax.imshow(Z.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2, interpolation='nearest')
    # gene names on x-axis
    ax.set_xticks(range(len(gorder)))
    ax.set_xticklabels(gorder, rotation=90, fontsize=4.0, fontfamily='monospace')
    ax.set_yticks([])  # chemical labels go on the sidebar (axL), not here
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(1.0)
    # dashed dividers
    ax.axvline(nU - 0.5, color='black', lw=1.6, ls='--')  # UP|DOWN
    ax.axhline(NCAT1 - 0.5, color='black', lw=1.8, ls='--')  # Cat1|Cat0

    # left category sidebar (chemical names placed on this axis, to the LEFT of the strip)
    catarr = np.array([[0]] * NCAT1 + [[1]] * (len(DISP) - NCAT1))
    axL.imshow(catarr, aspect='auto', cmap=mpl.colors.ListedColormap([COL_C1, COL_C0]))
    axL.set_xticks([])
    axL.set_yticks(range(len(DISP)))
    axL.set_yticklabels(DISP, fontsize=12)
    axL.tick_params(length=0)
    axL.text(0, (NCAT1 - 1) / 2, 'Cat.1  n = %d' % NCAT1, rotation=90, ha='center', va='center',
              color='white', fontsize=11, fontweight='bold')
    axL.text(0, NCAT1 + (len(DISP) - NCAT1 - 1) / 2, 'Cat.0  n = %d' % (len(DISP) - NCAT1), rotation=90, ha='center', va='center',
              color='white', fontsize=11, fontweight='bold')

    # top UP/DOWN header bar
    toparr = np.array([[0] * nU + [1] * nD])
    axT.imshow(toparr, aspect='auto', cmap=mpl.colors.ListedColormap([COL_C1, COL_C0]))
    axT.set_xticks([])
    axT.set_yticks([])
    axT.text((nU - 1) / 2, 0, 'UP-regulated  (n = %d)' % nU, ha='center', va='center',
              color='white', fontsize=15, fontweight='bold')
    axT.text(nU + (nD - 1) / 2, 0, 'DOWN-regulated  (n = %d)' % nD, ha='center', va='center',
              color='white', fontsize=15, fontweight='bold')

    cb = fig.colorbar(im, cax=cax)
    cb.set_label('Z-score (row-wise)', fontsize=14, fontweight='bold')
    cb.ax.tick_params(labelsize=12)
    cb.set_ticks([-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2])
    outpath = os.path.join(args.output_dir, "Fig3A_Heatmap_bigfont")
    for ext, kw in [('png', {'dpi': 600}),
                    ('tif', {'dpi': 600, 'pil_kwargs': {'compression': 'tiff_lzw', 'dpi': (600, 600)}}),
                    ('svg', {})]:
        fig.savefig(f"{outpath}.{ext}", bbox_inches='tight', **kw)
    plt.close(fig)
    print("chemicals", len(CHEMS), "genes", len(gorder), "up", nU, "down", nD, "Cat1", NCAT1, "| done")


if __name__ == "__main__":
    main()
