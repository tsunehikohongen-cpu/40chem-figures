#!/usr/bin/env python3
"""
Reproduce 'Pathway_Heatmap_GO_BP_Development.png' (Figure 2B) from scratch.

Part 1: Per-gene Delta-rlog (chemical vs. DMSO control mean) -> ssGSEA against
        GO_Biological_Process_2023 -> NES (normalized enrichment score) matrix
        (2,692 terms x 40 chemicals).
Part 2: Select the 30 development-related pathways and render the heatmap.
        Selection method (confirmed to exactly reproduce the published figure,
        both term membership and row order):
          1. Keyword-filter the 2,692 GO BP terms (case-insensitive substring
             match) against a fixed developmental-biology keyword list -> 291
             development-related terms.
          2. Rank those 291 terms by cross-chemical variance of NES; keep the
             top 30.
          3. Order the 30 rows by Ward-linkage / Euclidean-distance
             hierarchical clustering (scipy.cluster.hierarchy). Columns
             (chemicals) are kept in their fixed No.1-40 order (Table 1).

Inputs required (default: ./data/, override with --input-dir):
  - 40_Chem_DESeq_normalized.csv  (DESeq rlog-transformed expression matrix;
     columns = sample names, e.g. 'Acitretin_1', 'DMSO-1 0.1', etc.)

Outputs (default: ./output/, override with --output-dir):
  - ssgsea_NES_GO_Biological_Process_2023.csv  (2,692 GO BP terms x 40 chemicals;
    cached -- if this file already exists in --output-dir, Part 1 is skipped)
  - Pathway_Heatmap_GO_BP_Development.png / .pdf

Dependencies: python>=3.10, pandas, numpy, scipy, matplotlib, gseapy
  gseapy's ssGSEA call downloads the 'GO_Biological_Process_2023' gene set
  library from Enrichr on first use (internet access required unless the NES
  CSV cache already exists).
"""
import argparse
import os
import warnings

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.cluster.hierarchy import linkage, dendrogram

warnings.filterwarnings('ignore')

from matplotlib import font_manager
try:
    font_manager.findfont('Helvetica', fallback_to_default=False)
    rcParams['font.family'] = 'Helvetica'
except Exception:
    try:
        font_manager.findfont('Nimbus Sans', fallback_to_default=False)
        rcParams['font.family'] = 'Nimbus Sans'   # metric-compatible Helvetica clone (Linux/URW)
    except Exception:
        rcParams['font.family'] = 'sans-serif'
rcParams['axes.labelweight'] = 'bold'
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
rcParams['legend.frameon'] = False

# ===== 40 chemicals: (display name, data-column name), No.1-40 fixed order (Table 1) =====
HIGH_DOSE = [
    ('Acitretin', 'Acitretin_1'), ('Busulfan', 'Busulfan_1'), ('Carbamazepine', 'Carbamazepine_1'),
    ('Cyclophosphamide', 'Cyclophosphamide_1'), ('Fluconazole', 'Fluconazole_1'),
    ('5-Fluorouracil', '5-Fluorouracil_1'), ('Hydroxyurea', 'Hydroxyurea_1'),
    ('Isotretinoin', 'Isotretinoin_1'), ('Methotrexate', 'Methotrexate_1'),
    ('Phenytoin', 'Phenytoin_1'), ('Thalidomide', 'Thalidomide_1'), ('Topiramate', 'Topiramate_1'),
    ('Tretinoin', 'Tretinoin_1'), ('Valproic acid', 'Valproic acid_1'), ('Pomalidomide', 'Pomalidomide_1'),
    ('Aspirin', 'Aspirin_1'), ('Cytarabine', 'Cytarabine_0.3'), ('Ibuprofen', 'Ibuprofen_1'),
    ('Trimethadione', 'Trimethadione_1'), ('Vismodegib', 'Vismodegib_1'), ('Ibrutinib', 'Ibrutinib_1'),
    ('Pazopanib', 'Pazopanib_1'), ('Ribavirin', 'Ribavirin_1'), ('Tacrolimus', 'Tacrolimus_1'),
    ('Bosentan', 'Bosentan_1'), ('Cisplatin', 'Cisplatin_1'), ('Dabrafenib', 'Dabrafenib 1'),
    ('Dasatinib', 'Dasatinib 0.1'), ('Imatinib', 'Imatinib 1'), ('Rotenone', 'Rotenone 25'),
    ('Acetaminophen', 'Acetaminophen 1'), ('Acetylsalicylic acid', 'Acetylsalicylic acid 1'),
    ('D-Glucitol', 'D-Glucitol 1'), ('L-Ascorbic acid', 'L-Ascorbic acid 1'),
    ('Saccharin Sodium Salt hydrate', 'Saccharin Sodium Salt hydrate 1'), ('Deltamethrin', 'Deltamethrin 1'),
    ("2,2',4,4',5-Pentabromodiphenyl ether", "2,2',4,4',5-Pentabromodiphenyl ether 1"),
    ('Bisphenol A', 'Bisphenol A 1'), ('Chlorpyrifos', 'Chlorpyrifos 1'), ('imidacloprid', 'imidacloprid 1'),
]
CHEM_NAMES = [n for n, _ in HIGH_DOSE]
NO_LABELS = [f'No.{i + 1}  {c}' for i, c in enumerate(CHEM_NAMES)]

DMSO_COLS = ['DMSO-1 0.1', 'DMSO-1 0.01', 'DMSO-1 0.001',
             'DMSO-2 0.1', 'DMSO-2 0.01', 'DMSO-2 0.001']

KEYWORDS = ['development', 'differentiation', 'morphogenesis', 'embryo', 'organogenesis',
            'patterning', 'specification', 'commitment', 'stem cell', 'neural', 'cardiac',
            'angiogenesis', 'vasculo', 'germ layer', 'somite', 'axis specif', 'body plan',
            'gastrulation', 'neurogenesis', 'cardiogenesis', 'hematopoiesis']


def is_dev_related(term):
    t = term.lower()
    return any(k in t for k in KEYWORDS)


def shorten(name, maxlen=60):
    idx = name.find(' (GO:')
    if idx > 0:
        name = name[:idx]
    return name[:maxlen]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", default="data")
    p.add_argument("--output-dir", default="output")
    p.add_argument("--input-csv", default="40_Chem_DESeq_normalized.csv")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    nes_csv = os.path.join(args.output_dir, "ssgsea_NES_GO_Biological_Process_2023.csv")

    # ===== Part 1: compute ssGSEA NES matrix (skip if cached CSV already exists) =====
    if not os.path.exists(nes_csv):
        import gseapy as gp
        df = pd.read_csv(os.path.join(args.input_dir, args.input_csv))
        df = df.rename(columns={df.columns[0]: 'Gene'}).set_index('Gene')
        dmso_mean = df[DMSO_COLS].mean(axis=1)

        # Delta-rlog = rlog(chemical) - mean(rlog(DMSO_1..DMSO_6))  (approximates log2 FC)
        fc_mat = pd.DataFrame(index=df.index)
        for name, col in HIGH_DOSE:
            fc_mat[name] = df[col] - dmso_mean
        # exclude genes with |Delta-rlog| < 0.01 across ALL 40 chemicals (~21,494 genes retained)
        fc_mat = fc_mat[(fc_mat.abs() >= 0.01).any(axis=1)]

        ss = gp.ssgsea(
            data=fc_mat,
            gene_sets='GO_Biological_Process_2023',
            outdir=None,
            sample_norm_method='rank',
            no_plot=True,
            threads=4,
        )
        nes_df = ss.res2d.pivot(index='Term', columns='Name', values='NES')
        nes_df.reindex(columns=CHEM_NAMES).astype(float).to_csv(nes_csv)

    nes_df = pd.read_csv(nes_csv, index_col=0)

    # ===== Part 2a: keyword-filter to development-related GO BP terms =====
    dev_terms = [t for t in nes_df.index if is_dev_related(t)]
    print(f"Development-related GO BP terms (keyword filter): {len(dev_terms)}")  # -> 291

    # ===== Part 2b: rank by cross-chemical NES variance, keep top 30 =====
    dev_mat = nes_df.loc[dev_terms].astype(float)
    top30_terms = dev_mat.var(axis=1, skipna=True).nlargest(30).index.tolist()

    # ===== Part 2c: order rows by Ward / Euclidean hierarchical clustering =====
    mat30 = nes_df.reindex(columns=CHEM_NAMES).loc[top30_terms].astype(float)
    Z = linkage(mat30.values, method='ward', metric='euclidean')
    leaf_order = dendrogram(Z, no_plot=True)['leaves']
    row_order = [top30_terms[i] for i in leaf_order]

    mat = nes_df.reindex(columns=CHEM_NAMES).loc[row_order].astype(float)
    mat.index = [shorten(p) for p in mat.index]

    n_p, n_c = len(mat), len(mat.columns)
    fig, ax = plt.subplots(figsize=(24, 15),
                            gridspec_kw=dict(left=0.30, right=0.92, top=0.90, bottom=0.27))

    vabs = np.percentile(np.abs(mat.values), 95) or 1.0
    im = ax.imshow(mat.values, cmap='RdBu_r', vmin=-vabs, vmax=vabs,
                   aspect='auto', interpolation='nearest')

    ax.set_xticks(range(n_c))
    ax.set_xticklabels(NO_LABELS, rotation=45, ha='right', fontsize=15, fontweight='bold')
    ax.set_yticks(range(n_p))
    ax.set_yticklabels(mat.index, fontsize=16, fontweight='bold')
    ax.set_xlabel('Chemical (No.1-40)', fontweight='bold', fontsize=20)
    ax.set_ylabel('GO BP 2023', fontweight='bold', fontsize=20)
    ax.set_title(
        'ssGSEA enrichment - GO BP 2023 (high dose, 40 chemicals)\n'
        'Top 30 development-related pathways by cross-chemical variance;  NES = normalized enrichment score',
        fontweight='bold', fontsize=16, pad=12)

    ax.set_xticks(np.arange(-0.5, n_c, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_p, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=0.4)
    ax.tick_params(which='minor', bottom=False, left=False)

    cb = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.015, shrink=0.65)
    cb.set_label('NES', fontweight='bold', fontsize=15)
    cb.ax.tick_params(labelsize=13)

    outpath = os.path.join(args.output_dir, "Pathway_Heatmap_GO_BP_Development.png")
    fig.savefig(outpath, dpi=600, bbox_inches='tight', facecolor='white')
    fig.savefig(outpath.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    print(f"Saved: {outpath}")


if __name__ == "__main__":
    main()
