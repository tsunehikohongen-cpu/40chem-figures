#!/usr/bin/env python3
# =====================================================================
# Figure 3c  --  INGOR regulatory network of the 150 selected genes
#                (MCC-ranked nodes; activation/repression directed edges)
# Reproduction script (matches the published panel c).
#
# REQUIRED INPUT FILES (default: ./data/, override with --input-dir):
#   - INGOR_GeneSet150.csv   ('name' column = the 150 hub genes)
#   - 1000_ECv_edges.tsv     (source, target, edgeScore, gain ; INGOR ECv network)
# OUTPUT: Fig3C_Network_directed.{png,tif,svg} (600 dpi) in --output-dir (default: ./output)
# Dependencies: python>=3.10, pandas, numpy, networkx, scipy, matplotlib, adjustText
#   pip install adjustText
# =====================================================================
import argparse
import os
from math import factorial

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from adjustText import adjust_text

HUBN = 40  # number of top-MCC hub genes to label

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Nimbus Sans', 'TeX Gyre Heros', 'Arial', 'DejaVu Sans'],
    'svg.fonttype': 'path',
})


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", default="data")
    p.add_argument("--output-dir", default="output")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    inp = args.input_dir

    edges = pd.read_csv(os.path.join(inp, "1000_ECv_edges.tsv"), sep="\t")
    genes = [g.strip() for g in pd.read_csv(os.path.join(inp, "INGOR_GeneSet150.csv"))['name'].astype(str)]
    gset = set(genes)
    ind = edges[edges['source'].isin(gset) & edges['target'].isin(gset)]
    D = nx.DiGraph()
    D.add_nodes_from(genes)
    for _, r in ind.iterrows():
        D.add_edge(r['source'], r['target'], w=float(r['edgeScore']), gain=float(r['gain']))
    # bridge the isolated genes through one intermediate connector each
    U0 = D.to_undirected()
    iso = [n for n in genes if U0.degree(n) == 0]
    bridges = set()

    def addrow(r):
        D.add_edge(r['source'], r['target'], w=float(r['edgeScore']), gain=float(r['gain']))

    for n in iso:
        for _, r in edges[(edges['source'] == n) | (edges['target'] == n)].iterrows():
            o = r['target'] if r['source'] == n else r['source']
            lk = edges[((edges['source'] == o) & (edges['target'].isin(gset))) | ((edges['target'] == o) & (edges['source'].isin(gset)))]
            if len(lk):
                bridges.add(o)
                addrow(r)
                addrow(lk.iloc[0])
                break
    U = D.to_undirected()

    # Maximum Clique Centrality (CytoHubba MCC): sum over maximal cliques of (|C|-1)!
    mcc = {n: 0 for n in U.nodes()}
    for c in nx.find_cliques(U):
        f = factorial(len(c) - 1)
        for v in c:
            mcc[v] += f
    gmcc = sorted(genes, key=lambda x: mcc.get(x, 0), reverse=True)
    rank = {g: i + 1 for i, g in enumerate(gmcc)}

    pos = nx.kamada_kawai_layout(D)
    fig, ax = plt.subplots(figsize=(20, 18))
    cmap = mpl.colormaps['RdBu_r']
    ncolor = {g: cmap(1 - (rank[g] - 1) / (len(genes) - 1)) for g in genes}  # red = rank 1 (highest MCC)
    deg = dict(U.degree())
    nsize = {g: 300 + 70 * deg.get(g, 0) for g in genes}
    if bridges:
        nx.draw_networkx_nodes(D, pos, nodelist=list(bridges), node_shape='s', node_color='#E0E0E0',
                                edgecolors='#555555', linewidths=0.8, node_size=300, ax=ax)
    nx.draw_networkx_nodes(D, pos, nodelist=genes, node_color=[ncolor[g] for g in genes],
                            edgecolors='#333333', linewidths=0.7, node_size=[nsize[g] for g in genes], ax=ax)
    # curved directed arrows; colour by regulation sign (gain): red=activation, blue=repression
    for pos_sign, col in [(True, '#D6604D'), (False, '#4393C3')]:
        el = [(u, v) for u, v, dd in D.edges(data=True) if (dd.get('gain', 0) >= 0) == pos_sign]
        if not el:
            continue
        ws = [0.4 + 1.4 * D[u][v].get('w', 0.3) for u, v in el]
        nx.draw_networkx_edges(D, pos, edgelist=el, ax=ax, edge_color=col, width=ws, alpha=0.5,
                                arrows=True, arrowstyle='-|>', arrowsize=10, connectionstyle='arc3,rad=0.12',
                                node_size=600, min_source_margin=2, min_target_margin=5)
    # hub labels: black + white halo, de-cluttered
    hub = gmcc[:HUBN]
    texts = [ax.text(*pos[g], g, fontsize=12.5, fontweight='bold', ha='center', va='center', color='#111111',
                      zorder=7, path_effects=[pe.withStroke(linewidth=3.2, foreground='white')]) for g in hub]
    adjust_text(texts, ax=ax, expand=(1.4, 1.7), force_text=(0.4, 0.6),
                arrowprops=dict(arrowstyle='-', color='#7A7A7A', lw=0.7, shrinkA=4, shrinkB=4))
    ax.set_title('INGOR regulatory network of the %d selected genes (MCC-ranked)' % len(genes),
                 fontsize=21, fontweight='bold', pad=16)
    ax.axis('off')
    ax.margins(0.06)
    sm = mpl.cm.ScalarMappable(cmap='RdBu_r', norm=mpl.colors.Normalize(vmin=1, vmax=len(genes)))
    cb = fig.colorbar(sm, ax=ax, fraction=0.022, pad=0.01)
    cb.set_label('MCC rank  (1 = highest)', fontsize=15, fontweight='bold')
    cb.ax.invert_yaxis()
    cb.ax.tick_params(labelsize=12)
    leg = [Line2D([0], [0], color='#D6604D', lw=2.5, marker='>', markersize=8, label='Activation (gain >= 0)'),
           Line2D([0], [0], color='#4393C3', lw=2.5, marker='>', markersize=8, label='Repression (gain < 0)'),
           Line2D([0], [0], marker='s', ls='', markerfacecolor='#E0E0E0', markeredgecolor='#555555',
                  markersize=12, label='Bridge node (n=%d)' % len(bridges))]
    ax.legend(handles=leg, loc='lower center', fontsize=12, frameon=False, ncol=3)
    plt.tight_layout()
    outpath = os.path.join(args.output_dir, "Fig3C_Network_directed")
    for ext, kw in [('png', {'dpi': 600}),
                    ('tif', {'dpi': 600, 'pil_kwargs': {'compression': 'tiff_lzw', 'dpi': (600, 600)}}),
                    ('svg', {})]:
        fig.savefig(f"{outpath}.{ext}", bbox_inches='tight', **kw)
    plt.close(fig)
    print("nodes", D.number_of_nodes(), "edges", D.number_of_edges(), "bridges", bridges, "done")


if __name__ == "__main__":
    main()
