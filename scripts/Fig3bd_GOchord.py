#!/usr/bin/env python3
# =====================================================================
# Figure 3b/3d -- GO Biological Process "shared-gene" chord diagram
#   panel b: Gene Set 1 (Wilcoxon, 182 genes)
#   panel d: Gene Set 2 (INGOR, 150 genes)
# Reproduction script (matches the published panels). Self-contained: it
# performs the GO over-representation analysis (ORA) and then draws the
# chord diagram.
#
# REQUIRED INPUT FILES (default: ./data/, override with --input-dir):
#   - Wilcoxon_GeneSet182.csv / INGOR_GeneSet150.csv  (query gene list)
#   - GO_Biological_Process_2026.txt                  (Enrichr GMT-format library, tab-separated)
# OUTPUT: Fig_GOchord_Wilcoxon.{png,tif,svg} or Fig_GOchord_INGOR.{png,tif,svg}
#         (600 dpi) in --output-dir (default: ./output)
# Dependencies: python>=3.10, pandas, numpy, scipy, matplotlib, pycirclize
#   pip install pycirclize
#
# Usage:
#   python Fig3bd_GOchord.py --panel b   # Wilcoxon_GeneSet182.csv -> panel b
#   python Fig3bd_GOchord.py --panel d   # INGOR_GeneSet150.csv    -> panel d
# =====================================================================
import argparse
import os
import re
import textwrap

import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from pycirclize import Circos

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Nimbus Sans', 'TeX Gyre Heros', 'Arial', 'DejaVu Sans'],
    'svg.fonttype': 'path', 'figure.dpi': 300, 'savefig.dpi': 300,
})
PCMAP = LinearSegmentedColormap.from_list(
    'monago',
    ['#F39C12', '#E8C547', '#C9B037', '#8C8C3A', '#6E6E6E', '#3A4A6B', '#27408B', '#1F4FD8', '#3A6BFF'],
)

PANELS = {
    "b": dict(geneset="Wilcoxon_GeneSet182.csv", title="Wilcoxon gene set - GO BP (shared-gene chord)", outname="Fig_GOchord_Wilcoxon"),
    "d": dict(geneset="INGOR_GeneSet150.csv", title="INGOR gene set - GO BP (shared-gene chord)", outname="Fig_GOchord_INGOR"),
}


def wrap(s, width=24):
    return '\n'.join(textwrap.wrap(s, width=width)) or s


def clean(p):
    p = re.sub(r'\s*\(GO:\d+\)', '', p).strip()
    p = re.sub(r'\s*\(GO.*$', '', p).strip()
    return p


def parse_gmt(path):
    lib = {}
    for line in open(path, encoding='utf-8', errors='ignore'):
        x = line.rstrip("\n").split("\t")
        if len(x) >= 3:
            g = {t.split(',')[0].strip().upper() for t in x[2:] if t.strip()}
            if g:
                lib[x[0]] = g
    return lib


def ora(genes, lib, min_overlap=2, min_term=3):
    univ = set().union(*lib.values())
    N = len(univ)  # genome background = library gene universe
    q = {g.strip().upper() for g in genes} & univ
    n = len(q)
    rows = []
    for term, gs in lib.items():
        K = len(gs)
        k = len(q & gs)
        if K < min_term or k < min_overlap:
            continue
        p = hypergeom.sf(k - 1, N, K, n)
        rows.append((term, k, K, p, '; '.join(sorted(q & gs))))
    df = pd.DataFrame(rows, columns=['Pathway', 'Overlap', 'PathwaySize', 'p_value', 'Genes']).sort_values('p_value')
    m = len(df)
    fdr = df['p_value'].values * m / np.arange(1, m + 1)
    df['FDR'] = np.minimum(1, np.minimum.accumulate(fdr[::-1])[::-1])
    df['-log10FDR'] = -np.log10(df['FDR'].replace(0, np.nextafter(0, 1)))
    return df.reset_index(drop=True)


def chord(df, title, fname, output_dir, npool=22, nshow=14):
    d = df.head(npool).copy().reset_index(drop=True)
    d['lab'] = [clean(p) for p in d['Pathway']]
    d['gset'] = [set(str(g).split('; ')) for g in d['Genes']]
    n = len(d)
    S = np.zeros((n, n), int)
    for i in range(n):
        for j in range(i + 1, n):
            S[i, j] = S[j, i] = len(d['gset'][i] & d['gset'][j])
    deg = S.sum(1)
    keep = [i for i in range(n) if deg[i] > 0]
    keep = sorted(keep, key=lambda i: d['FDR'][i])[:nshow]
    keep = sorted(keep, key=lambda i: d['FDR'][i])
    lab2 = [wrap(d['lab'][i]) for i in keep]
    seen = {}
    for k, l in enumerate(lab2):
        if l in seen:
            lab2[k] = l + ' '
            seen[l] += 1
        else:
            seen[l] = 1
    M = pd.DataFrame(0, index=lab2, columns=lab2)
    for a, i in enumerate(keep):
        for b, j in enumerate(keep):
            if a != b:
                M.iloc[a, b] = S[i, j]
    M = M.loc[(M.sum(axis=1) > 0), (M.sum(axis=0) > 0)]
    kept = list(M.index)
    keep = [keep[lab2.index(l)] for l in kept]
    pv = np.array([d['p_value'][i] for i in keep])
    edges = np.linspace(pv.min(), pv.max(), 10)
    secol = {l: mpl.colors.to_hex(PCMAP((np.searchsorted(edges, p, side='right') - 1) / (len(edges) - 1)))
             for l, p in zip(kept, pv)}
    circos = Circos.chord_diagram(
        M, start=4, end=356, space=3, r_lim=(96, 100), cmap=secol,
        label_kws=dict(r=103, orientation='vertical', size=10.5, weight='bold'),
        link_kws=dict(color='#BFBFBF', alpha=0.55, ec='white', lw=0.3, direction=0),
    )
    fig = circos.plotfig(figsize=(13, 13))
    fig.suptitle(title, fontsize=16, fontweight='bold', y=1.05)
    handles = [Patch(facecolor=mpl.colors.to_hex(PCMAP(k / (len(edges) - 2))), edgecolor='none',
                      label='%.4f - %.4f' % (edges[k], edges[k + 1])) for k in range(len(edges) - 1)]
    fig.legend(handles=handles, title='p-value', loc='upper right', bbox_to_anchor=(1.0, 0.98),
               fontsize=9, title_fontproperties={'weight': 'bold', 'size': 10}, frameon=True,
               framealpha=0.95, edgecolor='#CCCCCC')
    for ext, kw in [('png', {'dpi': 600}),
                    ('tif', {'dpi': 600, 'pil_kwargs': {'compression': 'tiff_lzw', 'dpi': (600, 600)}}),
                    ('svg', {})]:
        fig.savefig(os.path.join(output_dir, f"{fname}.{ext}"), bbox_inches='tight', **kw)
    plt.close(fig)
    print("saved", fname, "| terms", len(kept))


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panel", choices=["b", "d"], required=True, help="b = Wilcoxon_GeneSet182 (panel b), d = INGOR_GeneSet150 (panel d)")
    p.add_argument("--input-dir", default="data")
    p.add_argument("--output-dir", default="output")
    p.add_argument("--go-gmt", default="GO_Biological_Process_2026.txt")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    cfg = PANELS[args.panel]
    genes = pd.read_csv(os.path.join(args.input_dir, cfg["geneset"])).iloc[:, 0].astype(str).tolist()
    GO = parse_gmt(os.path.join(args.input_dir, args.go_gmt))
    df = ora(genes, GO)
    chord(df, cfg["title"], cfg["outname"], args.output_dir)
    print("done")


if __name__ == "__main__":
    main()
