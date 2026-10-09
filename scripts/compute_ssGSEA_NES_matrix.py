#!/usr/bin/env python3
"""
Compute the ssGSEA normalized-enrichment-score (NES) matrix used by Fig. 2b, Fig. 5,
Fig. 6 and Supplementary Fig. 1 (2,692 GO Biological Process terms x 40 chemicals).

Per-gene Delta-rlog (chemical at its highest concentration vs. the mean of the six DMSO
vehicle controls) -> ssGSEA (gseapy) against the Enrichr library GO_Biological_Process_2023
-> NES matrix. Columns are the 40 chemicals in the fixed order
Category 1 (20) | Category 0 (15) | Unknown (5) (No. 1-40 of Table 1).

Input (default: ./data/, override with --input-dir / --input-csv):
  - 40_Chem_DESeq_normalized.csv  (DESeq-normalized, rlog-scale expression matrix;
     first column = gene, remaining columns = sample names, e.g. 'Acitretin_1', 'DMSO-1 0.1')

Output (default: data/NES_matrix_GO_BP_2023_40chemicals.csv, override with --output-csv):
  - NES matrix, 2,692 terms x 40 chemicals. A copy of the resulting file is shipped
    in this repository (data/NES_matrix_GO_BP_2023_40chemicals.csv).

Run (repository root):
  python scripts/compute_ssGSEA_NES_matrix.py

Notes
  - Genes with |Delta-rlog| < 0.01 in all 40 chemicals are excluded (21,494 genes retained).
  - gseapy downloads the 'GO_Biological_Process_2023' library from Enrichr on first use
    (internet access required).

Dependencies: python>=3.10, pandas, numpy, gseapy
"""
import argparse
import os
import warnings

import pandas as pd

warnings.filterwarnings('ignore')

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

DMSO_COLS = ['DMSO-1 0.1', 'DMSO-1 0.01', 'DMSO-1 0.001',
             'DMSO-2 0.1', 'DMSO-2 0.01', 'DMSO-2 0.001']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", default="data")
    p.add_argument("--input-csv", default="40_Chem_DESeq_normalized.csv")
    p.add_argument("--output-csv", default="data/NES_matrix_GO_BP_2023_40chemicals.csv")
    return p.parse_args()


def main():
    args = parse_args()
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
    print(f'ssGSEA input: {fc_mat.shape[0]} genes x {fc_mat.shape[1]} chemicals')

    ss = gp.ssgsea(
        data=fc_mat,
        gene_sets='GO_Biological_Process_2023',
        outdir=None,
        sample_norm_method='rank',
        no_plot=True,
        threads=4,
    )
    nes_df = ss.res2d.pivot(index='Term', columns='Name', values='NES')
    nes_df = nes_df.reindex(columns=CHEM_NAMES).astype(float)
    os.makedirs(os.path.dirname(args.output_csv) or '.', exist_ok=True)
    nes_df.to_csv(args.output_csv)
    print(f'saved {args.output_csv}  ({nes_df.shape[0]} terms x {nes_df.shape[1]} chemicals)')


if __name__ == '__main__':
    main()
