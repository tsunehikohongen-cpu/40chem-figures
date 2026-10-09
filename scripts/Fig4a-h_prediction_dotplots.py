"""
Dot-plot redesign of the 8 final-prediction figures (replacing bar charts).

Design:
  - Individual concentration replicates shown as small jittered gray dots
    (only where the source data actually contains per-concentration rows).
  - A thin horizontal line spans the replicate range (min-max).
  - The mean prediction is a large filled dot, colored by category
    (red = Category 1 / >=0.5, blue = Category 0 / <0.5).
  - Dashed vertical line at the 0.5 decision threshold.
  - Font: Helvetica (falls back to a metric-compatible clone / DejaVu Sans
    if unavailable on this system).

Inputs (default --input-dir data/fig4_predictions; run from the repository root).
Gene Set 1 = 182-gene Wilcoxon set, Gene Set 2 = 150-gene INGOR set:
  ML_multimodal_GeneSet182.csv        LR/SVM/RF probabilities per compound and concentration
  ML_multimodal_GeneSet150.csv
  ML_transcriptomics_GeneSet182.csv   LR/SVM/RF/Ensemble probabilities (index = compound_concentration)
  ML_transcriptomics_GeneSet150.csv
  NN_multimodal_GeneSet182.csv        NN output: 'Category' = predicted probability, 'MW' identifies the compound
  NN_multimodal_GeneSet150.csv
  NN_transcriptomics_GeneSet182.csv   NN output: 'Category' = predicted probability, 3 rows/compound
  NN_transcriptomics_GeneSet150.csv   (compound order: Deltamethrin, BDE-99, Bisphenol A, Chlorpyrifos, Imidacloprid)
"""
import argparse
import os

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

DPI = 600

plt.rcParams.update({
    'font.family':       ['Helvetica', 'Nimbus Sans', 'Arial', 'DejaVu Sans'],
    'font.size':         9,
    'axes.labelweight':  'normal',
    'axes.spines.top':   True,
    'axes.spines.right': True,
    'legend.frameon':    False,
})

CAT1_COLOR  = '#9e2a2b'   # deep maroon/red, matches reference figure
CAT0_COLOR  = '#2c3e66'   # navy blue, matches reference figure
REP_COLOR   = '#e69f00'   # Nature-style orange for individual concentrations
LINE_PALETTE = ['#1f77b4', '#2ca02c', '#9467bd', '#8c564b', '#17becf']  # one color per compound row
THRESHOLD   = 0.5

DISPLAY = {
    "2,2',4,4',5-Pentabromodiphenyl ether": "BDE-99",
    "imidacloprid": "Imidacloprid",
}
MW_MAP = {505.21: 'Deltamethrin', 564.69: 'BDE-99', 228.29: 'Bisphenol A',
          350.59: 'Chlorpyrifos', 255.66: 'Imidacloprid'}
COMPOUNDS = ['Deltamethrin', 'BDE-99', 'Bisphenol A', 'Chlorpyrifos', 'Imidacloprid']
LABELS = [c for c in COMPOUNDS for _ in range(3)]


def make_dotplot(rep_dict, title, stem, outdir):
    """
    rep_dict: {compound_display_name: list_of_replicate_values}
              (list may have length 1 if no replicate-level data exists)
    """
    # order by mean, descending
    means = {k: np.clip(np.mean(v), 0, None) for k, v in rep_dict.items()}
    order = sorted(means, key=lambda k: means[k], reverse=True)

    fig, ax = plt.subplots(figsize=(6.5, 3.0))

    for i, name in enumerate(order):
        vals = np.clip(np.array(rep_dict[name]), 0, None)
        mean = vals.mean()
        color = CAT1_COLOR if mean >= THRESHOLD else CAT0_COLOR

        if len(vals) > 1:
            # replicate range line (thin, one distinct color per compound)
            line_color = LINE_PALETTE[i % len(LINE_PALETTE)]
            ax.plot([vals.min(), vals.max()], [i, i],
                     color=line_color, linewidth=0.6, zorder=1)
            # jittered individual replicate dots (orange)
            jitter = np.linspace(-0.13, 0.13, len(vals))
            ax.scatter(vals, np.full(len(vals), i) + jitter,
                       s=16, color=REP_COLOR, alpha=0.9,
                       edgecolor='black', linewidth=0.3, zorder=2)

        # mean marker
        ax.scatter([mean], [i], s=75, color=color,
                   edgecolor='black', linewidth=0.5, zorder=3)
        ax.text(mean + 0.035, i, f'{mean:.3f}',
                va='center', ha='left', fontsize=9, fontweight='normal')

    ax.axvline(x=THRESHOLD, color='black', linestyle='--', linewidth=0.7, zorder=0)
    ax.set_xlim(-0.02, 1.18)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel('Prob (Category = 1)', fontweight='normal')

    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color('black')

    h1 = mlines.Line2D([], [], marker='o', linestyle='None', markersize=7,
                        markerfacecolor=CAT1_COLOR, markeredgecolor='black',
                        markeredgewidth=0.5, label='Category 1 (>= 0.5)')
    h0 = mlines.Line2D([], [], marker='o', linestyle='None', markersize=7,
                        markerfacecolor=CAT0_COLOR, markeredgecolor='black',
                        markeredgewidth=0.5, label='Category 0 (<0.5)')
    hrep = mlines.Line2D([], [], marker='o', linestyle='None', markersize=4,
                         markerfacecolor=REP_COLOR, markeredgecolor='black',
                         markeredgewidth=0.3, alpha=0.9, label='Individual concentration')
    ax.legend(handles=[h1, h0, hrep], loc='lower right', fontsize=7.5)

    ax.set_title(title, fontweight='normal', fontsize=10, pad=8)

    plt.tight_layout()
    for ext in ('png', 'tif'):
        path = os.path.join(outdir, f'{stem}.{ext}')
        plt.savefig(path, dpi=DPI, bbox_inches='tight')
        print(f'  Saved: {path}')
    plt.close()


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", default="data/fig4_predictions")
    p.add_argument("--output-dir", default="output/figures_dotplot")
    return p.parse_args()


def main():
    args = parse_args()
    base = args.input_dir
    outdir = args.output_dir
    os.makedirs(outdir, exist_ok=True)

    print('[1/8] ML x Multi-modal x GeneSet182')
    df = pd.read_csv(f'{base}/ML_multimodal_GeneSet182.csv')
    df['Chemical'] = df['Chemical'].replace(DISPLAY)
    df['row_mean'] = df[['LR_prob', 'SVM_prob', 'RF_prob']].mean(axis=1)
    reps = {name: g['row_mean'].tolist() for name, g in df.groupby('Chemical')}
    make_dotplot(reps, 'Unknown Chemical Prediction\nLabel A - GeneSet182 (Multi-modal)',
                 'ML_prediction_GS182_multimodal', outdir)

    print('[2/8] ML x Multi-modal x GeneSet150')
    df = pd.read_csv(f'{base}/ML_multimodal_GeneSet150.csv')
    df['Chemical'] = df['Chemical'].replace(DISPLAY)
    df['row_mean'] = df[['LR_prob', 'SVM_prob', 'RF_prob']].mean(axis=1)
    reps = {name: g['row_mean'].tolist() for name, g in df.groupby('Chemical')}
    make_dotplot(reps, 'Unknown Chemical Prediction\nLabel A - GeneSet150 (Multi-modal)',
                 'ML_prediction_GS150_multimodal', outdir)

    print('[3/8] ML x Transcriptomics x GeneSet182')
    df = pd.read_csv(f'{base}/ML_transcriptomics_GeneSet182.csv', index_col=0)
    df['Chemical'] = df['Chemical'].replace(DISPLAY)
    reps = {name: g['Ensemble'].tolist() for name, g in df.groupby('Chemical')}
    make_dotplot(reps, 'Unknown Chemical Prediction\nLabel A - GeneSet182 (Transcriptomics-based)',
                 'ML_prediction_GS182_transcriptomics', outdir)

    print('[4/8] ML x Transcriptomics x GeneSet150')
    df = pd.read_csv(f'{base}/ML_transcriptomics_GeneSet150.csv', index_col=0)
    df['Chemical'] = df['Chemical'].replace(DISPLAY)
    reps = {name: g['Ensemble'].tolist() for name, g in df.groupby('Chemical')}
    make_dotplot(reps, 'Unknown Chemical Prediction\nLabel A - GeneSet150 (Transcriptomics-based)',
                 'ML_prediction_GS150_transcriptomics', outdir)

    print('[5/8] NN x Multi-modal x GeneSet182')
    df = pd.read_csv(f'{base}/NN_multimodal_GeneSet182.csv')
    df['compound'] = df['MW'].round(2).map(MW_MAP)
    prob_col = df.columns[0]
    reps = {name: g[prob_col].tolist() for name, g in df.groupby('compound')}
    make_dotplot(reps, 'Unknown Chemical Prediction\nLabel A - GeneSet182 (Multi-modal)',
                 'NN_prediction_GS182_multimodal', outdir)

    print('[6/8] NN x Multi-modal x GeneSet150')
    df = pd.read_csv(f'{base}/NN_multimodal_GeneSet150.csv')
    df['compound'] = df['MW'].round(2).map(MW_MAP)
    prob_col = df.columns[0]
    reps = {name: g[prob_col].tolist() for name, g in df.groupby('compound')}
    make_dotplot(reps, 'Unknown Chemical Prediction\nLabel A - GeneSet150 (Multi-modal)',
                 'NN_prediction_GS150_multimodal', outdir)

    print('[7/8] NN x Transcriptomics x GeneSet182')
    df = pd.read_csv(f'{base}/NN_transcriptomics_GeneSet182.csv')
    df['compound'] = LABELS
    reps = {name: g['Category'].tolist() for name, g in df.groupby('compound')}
    make_dotplot(reps, 'Unknown Chemical Prediction\nLabel A - GeneSet182 (Transcriptomics-based)',
                 'NN_prediction_GS182_transcriptomics', outdir)

    print('[8/8] NN x Transcriptomics x GeneSet150')
    df = pd.read_csv(f'{base}/NN_transcriptomics_GeneSet150.csv')
    df['compound'] = LABELS
    reps = {name: g['Category'].tolist() for name, g in df.groupby('compound')}
    make_dotplot(reps, 'Unknown Chemical Prediction\nLabel A - GeneSet150 (Transcriptomics-based)',
                 'NN_prediction_GS150_transcriptomics', outdir)

    print(f'\nDone. All files saved to:\n  {outdir}')


if __name__ == "__main__":
    main()
