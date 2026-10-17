#!/usr/bin/env python3
"""Build setup_all.pkl, the input of run_shap_contribution_export.py.

For each of the four combinations {CACTUS (multi-modal), NonCACTUS (transcriptomics-only)} x
{GeneSet182, GeneSet150}, the LR / SVM / RF models of the manuscript are fitted on the complete
feature matrix of all 35 labelled chemicals (105 samples) after standardization, and stored together
with the scaled matrix and the feature names:

    {tag: {'X_scaled': DataFrame, 'models': {'LR','RF','SVM'}, 'feat_names': list}}

Inputs (from the 40chem-supplement repository; the files are the "Final_prediction" training inputs):
    <supplement>/data/multisigma_inputs/CACTUS_GeneSet{182,150}_Final_prediction/input_data1.csv
    (+ output_data1.csv = category label; the 10 descriptor columns are dropped for NonCACTUS)

Run:
    python scripts/build_setup_all_pkl.py --supplement ../40chem-supplement --out setup_all.pkl
"""
import argparse
import pickle

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

DESCRIPTORS = ['MW', 'LogP', 'TPSA', 'HBD', 'HBA', 'QED', 'BBB_penetrant', 'BBB_score',
               'Brenk_alerts', 'PAINS_alerts']


def make_models():
    return {
        'LR': LogisticRegression(C=0.1, max_iter=2000, class_weight='balanced', random_state=42),
        'RF': RandomForestClassifier(n_estimators=200, class_weight='balanced_subsample', random_state=42),
        'SVM': SVC(kernel='rbf', C=1.0, probability=True, class_weight='balanced', random_state=42),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--supplement', default='../40chem-supplement', help='path of the 40chem-supplement repository')
    ap.add_argument('--out', default='setup_all.pkl')
    args = ap.parse_args()

    out = {}
    for gs in ('182', '150'):
        d = f'{args.supplement}/data/multisigma_inputs/CACTUS_GeneSet{gs}_Final_prediction/'
        X = pd.read_csv(d + 'input_data1.csv')
        X.index = pd.read_csv(d + 'input_data1_rowkey.csv')['sample'].values
        y = pd.read_csv(d + 'output_data1.csv').iloc[:, 0].values
        for pre in ('CACTUS', 'NonCACTUS'):
            Xp = X if pre == 'CACTUS' else X.drop(columns=DESCRIPTORS)
            Xs = pd.DataFrame(StandardScaler().fit_transform(Xp.values), index=Xp.index, columns=Xp.columns)
            models = make_models()
            for m in models.values():
                m.fit(Xs.values, y)
            out[f'{pre}_GeneSet{gs}'] = {'X_scaled': Xs, 'models': models, 'feat_names': list(Xp.columns)}
            print(pre, f'GeneSet{gs}', Xs.shape)
    with open(args.out, 'wb') as f:
        pickle.dump(out, f)
    print('saved', args.out)


if __name__ == '__main__':
    main()
