# 40chem-figures

Reproduction scripts for the figures in the 40-chemical transcriptomics
analysis (PCA, GO-pathway heatmaps, gene-expression heatmaps, GO chord
diagrams, an INGOR regulatory network, and 8 chemical-prediction dot plots).

This repository contains **code only**. The underlying expression/label
data are not included; see "Inputs" under each script below for the exact
file names and column formats expected, so the data can be supplied
separately (e.g. via the associated Zenodo record, or your own copy of the
dataset).

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Tested with Python 3.12. `Fig3b`/`Fig3d` additionally need `pycirclize`,
`Fig3c` needs `adjustText`, and `2b` needs `gseapy` (all listed in
`requirements.txt`).

Place your input files under `data/` (gitignored) and run each script from
the repository root; outputs are written to `output/` (gitignored) unless
you pass `--output-dir`.

## Scripts

### `2a_reproduce_pca_high_concentration_600dpi.py` — Figure 2A
PCA of the 40 chemicals' high-concentration samples (Label A groups; DMSO
controls excluded), rendered at 600 dpi with PIL (no matplotlib).

- Inputs: `data/40_Chem_DESeq_normalized.csv` (DESeq-normalized expression
  matrix; first column = gene id, remaining columns = sample names like
  `Acitretin_1`, `DMSO-1 0.1`), `data/Table_Labels_A.xlsx` (sheet
  `Labels A_B`, header row 3, columns `Chemical` / `Label A`).
- Run: `python scripts/2a_reproduce_pca_high_concentration_600dpi.py`
- Outputs: `pca_high_concentration_600dpi.{png,tiff}`,
  `pca_high_concentration_scores.csv`, `pca_high_concentration_gene_loadings.csv`

### `2b_Pathway_Heatmap_GO_BP_Development_reproduce.py` — Figure 2B
Per-gene Delta-rlog (chemical vs. DMSO mean) → ssGSEA against
`GO_Biological_Process_2023` → NES matrix, then the top-30
development-related pathways by cross-chemical variance, row-ordered by
Ward hierarchical clustering.

- Inputs: `data/40_Chem_DESeq_normalized.csv`.
- Run: `python scripts/2b_Pathway_Heatmap_GO_BP_Development_reproduce.py`
- Note: the ssGSEA step (`gseapy`) downloads the `GO_Biological_Process_2023`
  gene set from Enrichr on first run (internet required), unless
  `output/ssgsea_NES_GO_Biological_Process_2023.csv` already exists, in
  which case Part 1 is skipped and the cached NES matrix is reused.
- Outputs: `ssgsea_NES_GO_Biological_Process_2023.csv`,
  `Pathway_Heatmap_GO_BP_Development.{png,pdf}`

### `Fig3a_heatmap.py` — Figure 3a
Expression heatmap (row Z-score) of the 182-gene Wilcoxon set across the 35
labeled chemicals, split into UP/DOWN-regulated blocks and Cat.1/Cat.0 rows.

- Inputs: `data/Top1000_data.csv` (columns: `symbol` + per-sample expression,
  sample names like `Acitretin_1`), `data/Wilcoxon_GeneSet182.csv` (column
  `Gene`), `data/1000_ECv_wilcoxon_nodes_LabelA_complete.tsv` (column `name`,
  `direction` = UP/DOWN/NS).
- Run: `python scripts/Fig3a_heatmap.py`
- Outputs: `Fig3A_Heatmap_bigfont.{png,tif,svg}`

### `Fig3bd_GOchord.py` — Figures 3b and 3d
GO Biological Process over-representation analysis (hypergeometric test +
BH-FDR) followed by a shared-gene chord diagram, for either gene set.

- Inputs: `data/Wilcoxon_GeneSet182.csv` (panel b) or
  `data/INGOR_GeneSet150.csv` (panel d) — single gene-name column;
  `data/GO_Biological_Process_2026.txt` (Enrichr GMT format, tab-separated).
- Run:
  `python scripts/Fig3bd_GOchord.py --panel b`
  `python scripts/Fig3bd_GOchord.py --panel d`
- Outputs: `Fig_GOchord_Wilcoxon.{png,tif,svg}` or `Fig_GOchord_INGOR.{png,tif,svg}`

### `Fig3c_network.py` — Figure 3c
INGOR regulatory network of the 150 selected genes; nodes ranked by Maximum
Clique Centrality (CytoHubba MCC), edges colored by regulation sign
(activation/repression).

- Inputs: `data/INGOR_GeneSet150.csv` (column `name`),
  `data/1000_ECv_edges.tsv` (columns `source`, `target`, `edgeScore`, `gain`).
- Run: `python scripts/Fig3c_network.py`
- Outputs: `Fig3C_Network_directed.{png,tif,svg}`

### `4a-h_reproduce_8_figures.py` — Figures 4a–4h
Dot plots (replacing bar charts) of chemical-category predictions from 4
model/gene-set combinations, each shown for both a multi-modal and a
transcriptomics-only feature set (8 panels total).

- Inputs: see the expected `--input-dir` layout documented in the script's
  module docstring (`LabelA/GeneSet182/Final_prediction.csv`,
  `ML_analysis/LabelA/GeneSet182/Final_prediction/predictions.csv`,
  `NN_Final_182.csv`, `NN_Final_182_transcriptomics.csv`, etc., and their
  GeneSet150 counterparts).
- Run: `python "scripts/4a-h_reproduce_8_figures.py"`
- Outputs: 8 dot-plot PNG/TIFF pairs under `output/figures_dotplot/`.

## License

See `LICENSE` (MIT). The code may be reused and adapted with attribution;
this license does not cover any dataset distributed separately.

## Citing

See `CITATION.cff`. If you use this code, please cite the associated
Zenodo record: [10.5281/zenodo.22108794](https://doi.org/10.5281/zenodo.22108794)
and/or the manuscript.
