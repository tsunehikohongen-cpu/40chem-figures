import sys,pickle
v=sys.argv[1]
src=open('/home/claude/40chem-supplement/scripts/ml_pipeline/run_gs1_gs2_analysis.py').read()
head=src.split('def run_fold_eval')[0]
g={'INPUT_DIR':f'/tmp/claude-0/ml/inputs_{v}','OUTPUT_DIR':f'/tmp/claude-0/ml/shap/tmp_{v}'}
exec(head,g)
import numpy as np,pandas as pd
from sklearn.preprocessing import StandardScaler
nos=g['LABELED_A']; y=np.array([g['label_lookup'][n] for n in nos for _ in range(len(g['get_logfc_cols'](n)))])
out={}
for pre,fn in (('CACTUS',g['build_cactus_X']),('NonCACTUS',g['build_logfc_X'])):
    for gs,genes in (('GeneSet182',g['genes_182']),('GeneSet150',g['genes_150'])):
        X=fn(nos,genes); sc=StandardScaler(); Xs=pd.DataFrame(sc.fit_transform(X.values),index=X.index,columns=X.columns)
        ms=g['make_models']()
        for m in ms.values(): m.fit(Xs.values,y)
        out[f'{pre}_{gs}']={'X_scaled':Xs,'models':ms,'feat_names':list(X.columns)}
        print(pre,gs,Xs.shape)
pickle.dump(out,open(f'setup_{v}.pkl','wb'))
