#!/usr/bin/env python3
"""Build Fig. 2 (composite): (a) PCA of DESeq-normalized expression,
   (b) within- vs between-category pairwise NES-profile correlation violin.

Inputs (run from the repository root)
------
- data/40_Chem_DESeq_normalized.csv : gene x sample DESeq-normalized (rlog-scale)
    expression matrix (126 samples: 40 chemicals x doses + 6 DMSO vehicles).
    Only the 40 high-concentration samples (one per chemical, see HIGH_DOSE
    below) are used; DMSO controls are excluded.
- data/NES_matrix_GO_BP_2023_40chemicals.csv : term x 40-chemical ssGSEA NES matrix
    (produced by compute_ssGSEA_NES_matrix.py). Column order is
    20 Category 1 | 15 Category 0 | 5 Unknown chemicals.

Method (panel a, as in the manuscript Methods)
------
The values are already normalized and transformed (they contain negative
values), so NO additional log transformation is applied. Genes with zero
variance across the 40 samples are removed, each gene is centred across
samples (not scaled), and PCA is applied (21,435 genes; PC1 23.09%, PC2 11.29%).

Output
------
- output/Fig2_PCA_and_NES_correlation.png / .pdf (600 dpi; panels a / b stacked)
"""
import itertools
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Ellipse, Polygon
from matplotlib.lines import Line2D
from fig_style import apply_style, DPI, CATEGORY_COLORS, WITHIN_COLOR, BETWEEN_COLOR

apply_style()

DESEQ_CSV = 'data/40_Chem_DESeq_normalized.csv'
NES_CSV   = 'data/NES_matrix_GO_BP_2023_40chemicals.csv'
OUT_DIR   = 'output'

N_CAT1, N_CAT0, N_UNK = 20, 15, 5
CATEGORY = ['Cat1']*N_CAT1 + ['Cat0']*N_CAT0 + ['Unknown']*N_UNK

# 40 chemicals: (display name, high-concentration sample column), No.1-40 fixed
# order (Category 1 | Category 0 | Unknown)
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


# =============================================================================
# (a) PCA of DESeq-normalized expression (high-concentration samples; DMSO
#     controls excluded)
# =============================================================================
def panel_a_pca(ax, expr_csv):
    expr_all = pd.read_csv(expr_csv, index_col=0)             # genes x samples
    expr = expr_all[[col for _, col in HIGH_DOSE]]            # 40 high-dose samples
    cats = CATEGORY

    # Gene-wise centred PCA on the normalized values (no extra log transform)
    expr = expr.loc[expr.var(axis=1) > 0]                     # drop zero-variance genes
    X = expr.values.T.astype(float)                           # samples x genes
    X -= X.mean(axis=0, keepdims=True)                        # gene-wise centre
    pca = PCA(n_components=2)
    pc = pca.fit_transform(X)
    var = pca.explained_variance_ratio_ * 100
    print(f'Panel a: {expr.shape[0]} genes, PC1={var[0]:.2f}%, PC2={var[1]:.2f}%')

    MARK = {'Cat1': ('o', 55), 'Cat0': ('s', 55), 'Unknown': ('^', 70)}
    LABEL = {'Cat1': 'Category 1', 'Cat0': 'Category 0', 'Unknown': 'Unknown'}
    for g in ['Cat1', 'Cat0', 'Unknown']:
        idx = [i for i, c in enumerate(cats) if c == g]
        if not idx: continue
        m, s = MARK[g]
        ax.scatter(pc[idx, 0], pc[idx, 1],
                   marker=m, s=s, c=CATEGORY_COLORS[g],
                   edgecolor='black', lw=0.5, zorder=3,
                   label=f'{LABEL[g]} (n={len(idx)})')

        # 1-SD confidence ellipse per group
        xs, ys = pc[idx, 0], pc[idx, 1]
        if len(idx) >= 3:
            cov = np.cov(xs, ys)
            vals, vecs = np.linalg.eigh(cov)
            order = vals.argsort()[::-1]; vals, vecs = vals[order], vecs[:, order]
            theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
            w, h = 2 * np.sqrt(vals)
            ax.add_patch(Ellipse((xs.mean(), ys.mean()), w, h, angle=theta,
                                 facecolor=CATEGORY_COLORS[g], alpha=0.18,
                                 edgecolor='none', zorder=1))
            # convex hull (dashed)
            try:
                from scipy.spatial import ConvexHull
                hull = ConvexHull(np.c_[xs, ys])
                poly = Polygon(np.c_[xs, ys][hull.vertices], closed=True,
                               facecolor=CATEGORY_COLORS[g], alpha=0.10,
                               edgecolor=CATEGORY_COLORS[g], lw=1.0, ls='--', zorder=2)
                ax.add_patch(poly)
            except Exception:
                pass

    ax.axhline(0, color='grey', lw=0.6); ax.axvline(0, color='grey', lw=0.6)
    ax.set_xlabel(f'PC1 ({var[0]:.1f}%)')
    ax.set_ylabel(f'PC2 ({var[1]:.1f}%)')
    ax.set_title('a  PCA — DESeq-normalized expression (DMSO excluded)',
                 fontsize=14, loc='left', pad=10)
    ax.legend(loc='upper left', fontsize=10, frameon=False)


# =============================================================================
# (b) Within- vs between-category pairwise NES-profile Pearson correlation
#     (35 Category 1 + Category 0 chemicals; Unknown excluded)
# =============================================================================
def panel_b_violin(ax, nes_csv):
    nes = pd.read_csv(nes_csv, index_col=0)
    assert list(nes.columns) == CHEM_NAMES, \
        'NES matrix columns must be the 40 chemicals in the fixed order Category 1 | Category 0 | Unknown'
    X = nes.values
    lab_idx = [i for i, c in enumerate(CATEGORY) if c != 'Unknown']
    lab_cat = [CATEGORY[i] for i in lab_idx]
    corr = np.corrcoef(X[:, lab_idx].T)

    within, between = [], []
    for a, b in itertools.combinations(range(len(lab_idx)), 2):
        (within if lab_cat[a] == lab_cat[b] else between).append(corr[a, b])
    within = np.asarray(within); between = np.asarray(between)
    dmean = within.mean() - between.mean()
    _, p = mannwhitneyu(within, between, alternative='two-sided')
    print(f'Panel b: within n={len(within)} mean={within.mean():.3f} | '
          f'between n={len(between)} mean={between.mean():.3f} | '
          f'Δmean={dmean:+.3f}  MWU p={p:.4f}')

    parts = ax.violinplot([within, between], positions=[0, 1],
                          showmeans=False, showextrema=False, widths=0.80)
    for pc, col in zip(parts['bodies'], [WITHIN_COLOR, BETWEEN_COLOR]):
        pc.set_facecolor(col); pc.set_edgecolor('none'); pc.set_alpha(0.9)
    for i, d in enumerate([within, between]):
        ax.plot([i - 0.3, i + 0.3], [d.mean()] * 2, 'k-', lw=2.5)

    ax.axhline(0, ls='--', color='grey', lw=0.9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'Within-category\n(n={len(within)})',
                        f'Between-category\n(n={len(between)})'], fontsize=13)
    ax.set_ylabel('Pairwise NES-profile correlation r', fontsize=14)
    ax.set_title(f'b  Δmean = {dmean:+.3f}  |  MWU P = {p:.4f}',
                 fontsize=14, loc='left', pad=10)


# =============================================================================
# Compose
# =============================================================================
fig = plt.figure(figsize=(8, 12))
gs = gridspec.GridSpec(2, 1, height_ratios=[1.0, 1.1], hspace=0.30,
                       left=0.11, right=0.97, top=0.96, bottom=0.06)
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1])
panel_a_pca(ax_a, DESEQ_CSV)
panel_b_violin(ax_b, NES_CSV)

os.makedirs(OUT_DIR, exist_ok=True)
for ext in ('png', 'pdf'):
    path = f'{OUT_DIR}/Fig2_PCA_and_NES_correlation.{ext}'
    fig.savefig(path, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f'saved {path}')
plt.close(fig)
