# 40chem-figures

Code and input data for the main-text figures of the 40-chemical iPSC
transcriptomics manuscript (mechanism-of-action structure of human iPSC
responses and developmental-toxicity prediction): PCA and NES-profile
correlation (Fig. 2), feature-gene selection (Fig. 3), predictions for the
five Unknown chemicals (Fig. 4) and read-across of the Unknown chemicals
(Fig. 6).

Companion repository (Supplementary Information code and data):
[40chem-supplement](https://github.com/tsunehikohongen-cpu/40chem-supplement)

## Which script makes which figure?

| Figure | Script | Main inputs (in `data/`) |
|---|---|---|
| Fig. 2a, 2b — PCA; within- vs between-category NES correlation | `scripts/Fig2_PCA_and_NES_correlation.py` | `40_Chem_DESeq_normalized.csv`, `NES_matrix_GO_BP_2023_40chemicals.csv` |
| Fig. 3a — expression heatmap of Gene Set 1 | `scripts/Fig3a_heatmap.py` | see below |
| Fig. 3b, 3d — GO chord diagrams | `scripts/Fig3bd_GOchord.py` | see below |
| Fig. 3c — INGOR network of Gene Set 2 | `scripts/Fig3c_network.py` | see below |
| Fig. 4a-h — predictions for the Unknown chemicals | `scripts/Fig4a-h_prediction_dotplots.py` | see below |
| Fig. 6 — read-across of the Unknown chemicals | `scripts/Fig6_read_across.py` | `NES_matrix_GO_BP_2023_40chemicals.csv` |
| NES matrix used by Fig. 2b and Fig. 6 | `scripts/compute_ssGSEA_NES_matrix.py` | `40_Chem_DESeq_normalized.csv` |

`scripts/fig_style.py` holds the shared matplotlib style and colour palette.
Figures of the Supplementary Information are made by the code in
[40chem-supplement](https://github.com/tsunehikohongen-cpu/40chem-supplement).

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Tested with Python 3.12. `Fig3b`/`Fig3d` additionally need `pycirclize`,
`Fig3c` needs `adjustText`, and `compute_ssGSEA_NES_matrix.py` needs `gseapy`
(all listed in `requirements.txt`).

Run each script from the repository root; outputs are written to `output/`
(git-ignored) unless the script says otherwise.

## Data (`data/`)

- **`40_Chem_DESeq_normalized.csv`** — DESeq-normalized (rlog-scale) expression
  matrix, 57,923 genes x 126 samples (40 chemicals at three concentrations plus six
  DMSO vehicle controls). The raw RNA-seq data are deposited in GEO (GSE345109).
  Sample columns are named like `Acitretin_1` or `DMSO-1 0.1`.
- **`NES_matrix_GO_BP_2023_40chemicals.csv`** — ssGSEA normalized enrichment
  scores (2,692 GO Biological Process terms x 40 chemicals; columns in the fixed
  order Category 1 (20) | Category 0 (15) | Unknown (5)). Derived from the matrix
  above by `scripts/compute_ssGSEA_NES_matrix.py`; re-running that script
  reproduces the shipped file.
- The gene lists, INGOR network tables and model-prediction files used by Fig. 3
  and Fig. 4 are not included; their expected names and formats are described under
  each script below.

## Scripts

### `Fig2_PCA_and_NES_correlation.py` — Figure 2
(a) PCA of the 40 high-concentration samples (DMSO controls excluded; normalized
values used as they are, genes with zero variance removed, gene-wise centring;
21,435 genes, PC1 23.09%, PC2 11.29%). (b) Pairwise Pearson correlation of the
ssGSEA NES profiles of the 35 labelled chemicals, within- vs between-category
(Mann-Whitney U test).

- Run: `python scripts/Fig2_PCA_and_NES_correlation.py`
- Output: `output/Fig2_PCA_and_NES_correlation.{png,pdf}`

### `compute_ssGSEA_NES_matrix.py` — ssGSEA NES matrix
Delta-rlog (chemical at its highest concentration minus the mean of the six DMSO
controls) -> ssGSEA (`gseapy`, library `GO_Biological_Process_2023` from Enrichr,
internet required) -> NES matrix (2,692 terms x 40 chemicals). Genes with
|delta-rlog| < 0.01 in all 40 chemicals are excluded (21,494 genes).

- Run: `python scripts/compute_ssGSEA_NES_matrix.py`
- Output: `data/NES_matrix_GO_BP_2023_40chemicals.csv`

### `Fig6_read_across.py` — Figure 6
Cosine similarity of the ssGSEA NES profile (2,692 GO BP terms) of each of the
five Unknown chemicals to an in-group anchor set (valproic acid, tretinoin,
isotretinoin, methotrexate, thalidomide; Category 1) and an out-group anchor set
(L-ascorbic acid, saccharin sodium salt hydrate, D-glucitol, acetylsalicylic
acid; Category 0 non-toxic references): (a) mean similarity to each set,
(b) 5 x 9 heatmap of the individual similarities.

- Run: `python scripts/Fig6_read_across.py`
- Output: `output/Fig6_read_across.{png,pdf}`

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

### `Fig4a-h_prediction_dotplots.py` — Figures 4a–4h
Dot plots (replacing bar charts) of chemical-category predictions from 4
model/gene-set combinations, each shown for both a multi-modal and a
transcriptomics-only feature set (8 panels total).

- Inputs: see the expected `--input-dir` layout documented in the script's
  module docstring (`LabelA/GeneSet182/Final_prediction.csv`,
  `ML_analysis/LabelA/GeneSet182/Final_prediction/predictions.csv`,
  `NN_Final_182.csv`, `NN_Final_182_transcriptomics.csv`, etc., and their
  GeneSet150 counterparts).
- Run: `python scripts/Fig4a-h_prediction_dotplots.py`
- Outputs: 8 dot-plot PNG/TIFF pairs under `output/figures_dotplot/`.

## License

See `LICENSE` (MIT). The code may be reused and adapted with attribution;
this license does not cover any dataset distributed separately.

## Citing

See `CITATION.cff`. If you use this code, please cite the associated
Zenodo record: [10.5281/zenodo.22108793](https://doi.org/10.5281/zenodo.22108793)
(concept DOI; always resolves to the latest version) and/or the manuscript.
