# 40chem-figures

Code and input data for the main-text figures of the 40-chemical iPSC
transcriptomics manuscript (mechanism-of-action structure of human iPSC
responses and developmental-toxicity prediction): PCA and NES-profile
correlation (Fig. 2), feature-gene selection (Fig. 3), predictions for the
five Unknown chemicals (Fig. 4), mechanism-of-action (MoA) structure of the
transcriptional response (Fig. 5) and read-across of the Unknown chemicals
(Fig. 6).

Companion repository (Supplementary Information code and data):
[40chem-supplement](https://github.com/tsunehikohongen-cpu/40chem-supplement)

## Which script makes which figure?

| Figure | Script | Main inputs (in `data/`) |
|---|---|---|
| Fig. 2a, 2b — PCA; within- vs between-category NES correlation | `scripts/Fig2_PCA_and_NES_correlation.py` | `40_Chem_DESeq_normalized.csv`, `NES_matrix_GO_BP_2023_40chemicals.csv` |
| Fig. 3a — expression heatmap of Gene Set 1 | `scripts/Fig3a_heatmap.py` | `Top1000_data.csv`, `Wilcoxon_GeneSet182.csv`, `1000_ECv_wilcoxon_nodes_LabelA_complete.tsv` |
| Fig. 3b, 3d — GO chord diagrams | `scripts/Fig3bd_GOchord.py` | `Wilcoxon_GeneSet182.csv` / `INGOR_GeneSet150.csv`, `GO_Biological_Process_2026.txt` |
| Fig. 3c — INGOR network of Gene Set 2 | `scripts/Fig3c_network.py` | `INGOR_GeneSet150.csv`, `1000_ECv_edges.tsv` |
| Fig. 4a-h — predictions for the Unknown chemicals | `scripts/Fig4a-h_prediction_dotplots.py` | the 8 files in `fig4_predictions/` |
| Fig. 5 — MoA structure (clustering, MoA vs category, per-MoA coherence) | `scripts/Fig5_MoA_similarity.py` | `NES_matrix_GO_BP_2023_40chemicals.csv` |
| Fig. 5 / Supplementary Table 7 — chemical-level permutation test | `scripts/Fig5_MoA_permutation_test.py` | `NES_matrix_GO_BP_2023_40chemicals.csv` |
| Fig. 6 — read-across of the Unknown chemicals | `scripts/Fig6_read_across.py` | `NES_matrix_GO_BP_2023_40chemicals.csv` |
| Per-feature contributions / modality attribution (Results) | `scripts/run_shap_contribution_export.py` | `results/factor_analysis/*.csv` (see below) |
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
- **Fig. 3 inputs**
  - `Top1000_data.csv` — expression of the top-1,000 variable genes (iDEP; column
    `symbol` plus one column per sample, e.g. `Acitretin_1`).
  - `Wilcoxon_GeneSet182.csv` — Gene Set 1, the 182 genes selected by INGOR edge
    score >= 0.5 and Wilcoxon signed-rank test FDR <= 0.01 (column `Gene`).
  - `INGOR_GeneSet150.csv` — Gene Set 2, the top-150 hub genes by Maximum Clique
    Centrality (column `name`).
  - `1000_ECv_edges.tsv`, `1000_ECv_wilcoxon_nodes_LabelA_complete.tsv` — INGOR
    network edges (`source, target, edgeScore, gain`) and per-gene node statistics
    (`name`, `direction` = UP/DOWN/NS).
  - `GO_Biological_Process_2026.txt` — Enrichr `GO_Biological_Process_2026` library
    (GMT format), used for the over-representation analysis.
- **Fig. 4 inputs (`fig4_predictions/`)** — predicted probability of Category 1
  for the five Unknown chemicals at three concentrations each, for the four
  model / gene-set combinations in the multi-modal and the transcriptomics-only
  setting (Gene Set 1 = 182 genes, Gene Set 2 = 150 genes):
  `ML_multimodal_GeneSet{182,150}.csv`, `ML_transcriptomics_GeneSet{182,150}.csv`
  (LR / SVM / RF probabilities; the ML value plotted is their mean) and
  `NN_multimodal_GeneSet{182,150}.csv`, `NN_transcriptomics_GeneSet{182,150}.csv`
  (neural-network output: `Category` = predicted probability followed by the input
  features; rows are ordered Deltamethrin, BDE-99, Bisphenol A, Chlorpyrifos,
  Imidacloprid, three concentrations each).

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

### `Fig5_MoA_similarity.py` — Figure 5
Cosine similarity of the ssGSEA NES profiles of the 35 labelled chemicals:
(a) hierarchically clustered similarity matrix with category and MoA-subgroup
bars (16 MoA subgroups); (b) same MoA vs same category / different MoA vs
different category; (c, d) same-MoA vs different-MoA pairs within Category 1 and
Category 0; (e) per-MoA coherence (Cliff's delta). Permutation P values shown in
the figure come from `Fig5_MoA_permutation_test.py`.

- Run: `python scripts/Fig5_MoA_similarity.py`
- Output: `output/Fig5_MoA_similarity.{png,pdf}`

### `Fig5_MoA_permutation_test.py` — chemical-level permutation test (Fig. 5, Supplementary Table 7)
MoA labels are shuffled among chemicals within each category (10,000
permutations, fixed seed 42) and the test statistic (Cliff's delta or mean
difference) recomputed; P_perm is the fraction of permutations with a statistic
at least as large as the observed one. Results are printed to the console.

- Run: `python scripts/Fig5_MoA_permutation_test.py` (a few minutes)

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

- Inputs: the 8 files in `data/fig4_predictions/` (see Data above and the
  script's module docstring).
- Run: `python scripts/Fig4a-h_prediction_dotplots.py`
- Outputs: 8 dot-plot PNG/TIFF pairs under `output/figures_dotplot/`.

### `results/factor_analysis/` — per-feature contribution tables
Per-feature contribution percentages (`Contribution [%]`, `Contribution_positive [%]`,
`Contribution_negative [%]`; each file sums to 100%) for the four combinations
CACTUS (multi-modal) / NonCACTUS (transcriptomics-only) x Gene Set 1 (182 genes) /
Gene Set 2 (150 genes), for both model families:

- `{CACTUS,NonCACTUS}_ML_GeneSet{182,150}.csv` — ML (LR / SVM / RF ensemble), SHAP-based,
  produced by `scripts/run_shap_contribution_export.py`.
- `{CACTUS,NonCACTUS}_NN_GeneSet{182,150}.csv` — neural network, factor analysis of the
  first-layer activations (produced by a separate NN pipeline that is not included).

The "modality attribution" in the Results (chemical descriptors = 11-19% of the
discriminative weight in the ML classifiers, 5-6% in the neural network) is the summed
contribution of the 10 descriptor rows (MW, LogP, TPSA, HBD, HBA, QED, BBB_penetrant,
BBB_score, Brenk_alerts, PAINS_alerts) in the CACTUS files: ML 11.2% (GeneSet182) and
19.0% (GeneSet150); NN 5.6% (GeneSet182) and 5.4% (GeneSet150). The NonCACTUS files
contain no descriptors.

Reproducing the ML contributions (needs a clone of
[40chem-supplement](https://github.com/tsunehikohongen-cpu/40chem-supplement) next to this repository):

```bash
python scripts/build_setup_all_pkl.py --supplement ../40chem-supplement --out setup_all.pkl
python scripts/run_shap_contribution_export.py --setup-pkl setup_all.pkl --time-budget 120   # repeat until "ALL COMPLETE"
```

`build_setup_all_pkl.py` fits LR / SVM / RF on the 105-sample final training matrices stored in the supplement repository.
The regenerated tables agree with the shipped ones to within about 0.7 percentage points per feature (SHAP sampling and random
forest are stochastic) and give the same modality attribution (descriptors 11.8 % / 19.5 % for Gene Set 1 / 2 versus 11.8 % /
19.3 % shipped). The neural-network tables were produced with Multi-Sigma and are provided as results only.

## License

See `LICENSE` (MIT). The code may be reused and adapted with attribution;
this license does not cover any dataset distributed separately.

## Citing

See `CITATION.cff`. If you use this code, please cite the associated
Zenodo record: [10.5281/zenodo.22108793](https://doi.org/10.5281/zenodo.22108793)
(concept DOI; always resolves to the latest version) and/or the manuscript.
