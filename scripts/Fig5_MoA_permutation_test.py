#!/usr/bin/env python3
"""Chemical-level permutation test for the MoA structure analyses in Fig. 5
(results reported in Supplementary Table 7).
Null: MoA labels randomly shuffled among chemicals within each category.
Observed statistic preserved; P_perm = Pr(shuffled stat >= observed) under null.

Run from the repository root (results are printed to the console; fixed seed 42):
  python scripts/Fig5_MoA_permutation_test.py
Input: data/NES_matrix_GO_BP_2023_40chemicals.csv"""
import numpy as np, pandas as pd
from itertools import combinations
from sklearn.metrics.pairwise import cosine_similarity

nes=pd.read_csv('data/NES_matrix_GO_BP_2023_40chemicals.csv', index_col=0)
MOA={'Acitretin':'Retinoid signalling','Isotretinoin':'Retinoid signalling','Tretinoin':'Retinoid signalling',
'Busulfan':'Alkylating DNA damage','Cyclophosphamide':'Alkylating DNA damage',
'5-Fluorouracil':'Antimetabolite DNA synthesis','Hydroxyurea':'Antimetabolite DNA synthesis',
'Methotrexate':'Antimetabolite DNA synthesis','Cytarabine':'Antimetabolite DNA synthesis',
'Carbamazepine':'Antiepileptic ion channel','Phenytoin':'Antiepileptic ion channel',
'Topiramate':'Antiepileptic ion channel','Trimethadione':'Antiepileptic ion channel',
'Valproic acid':'HDAC inhibition antiepileptic','Thalidomide':'Cereblon IMiD','Pomalidomide':'Cereblon IMiD',
'Aspirin':'COX inhibition','Ibuprofen':'COX inhibition','Fluconazole':'Azole antifungal CYP',
'Vismodegib':'Hedgehog pathway inhibition',
'Ibrutinib':'Kinase inhibitor','Pazopanib':'Kinase inhibitor','Dabrafenib':'Kinase inhibitor',
'Dasatinib':'Kinase inhibitor','Imatinib':'Kinase inhibitor',
'Ribavirin':'Antiviral nucleoside analog','Tacrolimus':'Calcineurin immunosuppressant',
'Bosentan':'Endothelin receptor antagonist','Cisplatin':'Platinum DNA crosslinking',
'Rotenone':'Mitochondrial complex I inhibition',
'Acetaminophen':'Non toxic reference','Acetylsalicylic acid':'Non toxic reference',
'D-Glucitol':'Non toxic reference','L-Ascorbic acid':'Non toxic reference',
'Saccharin Sodium Salt hydrate':'Non toxic reference'}
CAT1=['Acitretin','Busulfan','Carbamazepine','Cyclophosphamide','Fluconazole','5-Fluorouracil','Hydroxyurea','Isotretinoin','Methotrexate','Phenytoin','Thalidomide','Topiramate','Tretinoin','Valproic acid','Pomalidomide','Aspirin','Cytarabine','Ibuprofen','Trimethadione','Vismodegib']
CAT0=['Ibrutinib','Pazopanib','Ribavirin','Tacrolimus','Bosentan','Cisplatin','Dabrafenib','Dasatinib','Imatinib','Rotenone','Acetaminophen','Acetylsalicylic acid','D-Glucitol','L-Ascorbic acid','Saccharin Sodium Salt hydrate']
NON_TOX='Non toxic reference'

cols=[c for c in nes.columns if c in CAT1+CAT0]
sim=pd.DataFrame(cosine_similarity(nes[cols].T.values), index=cols, columns=cols)

def cliffs(a,b):
    if len(a)==0 or len(b)==0: return np.nan
    a=np.asarray(a); b=np.asarray(b)
    g=np.sum(a[:,None]>b[None,:]); l=np.sum(a[:,None]<b[None,:])
    return (g-l)/(len(a)*len(b))

def within_cat_sm_dm(members, moa):
    sm,dm=[],[]
    for a,b in combinations(members,2):
        if moa.get(a)==NON_TOX or moa.get(b)==NON_TOX: continue
        r=sim.loc[a,b]
        if moa[a]==moa[b]: sm.append(r)
        else: dm.append(r)
    return np.array(sm), np.array(dm)

rng=np.random.default_rng(42)
N_PERM=10000

# === Analysis 1: Within Cat 1, same-MoA vs diff-MoA ===
sm_obs, dm_obs = within_cat_sm_dm(CAT1, MOA)
delta_obs = cliffs(sm_obs, dm_obs)
mean_diff_obs = sm_obs.mean()-dm_obs.mean()
print(f'=== Within Cat1 ===')
print(f'Observed: same n={len(sm_obs)}, diff n={len(dm_obs)}, Cliff\'s δ={delta_obs:+.3f}, mean_diff={mean_diff_obs:+.4f}')

# Permutation: shuffle MoA labels among Cat1 chemicals
cat1_moas=[MOA[c] for c in CAT1]
delta_null=[]; md_null=[]
for _ in range(N_PERM):
    perm=list(cat1_moas); rng.shuffle(perm)
    moa_shuf={c:m for c,m in zip(CAT1,perm)}
    for c in CAT0: moa_shuf[c]=MOA[c]
    sm,dm=within_cat_sm_dm(CAT1, moa_shuf)
    if len(sm)>0 and len(dm)>0:
        delta_null.append(cliffs(sm,dm))
        md_null.append(sm.mean()-dm.mean())
delta_null=np.array(delta_null); md_null=np.array(md_null)
p_perm_delta=(delta_null>=delta_obs).mean()
p_perm_md=(md_null>=mean_diff_obs).mean()
print(f'Permutation (N={len(delta_null)}): P_perm(δ)={p_perm_delta:.4f}, P_perm(mean_diff)={p_perm_md:.4f}')
print(f'Null δ: mean={delta_null.mean():+.3f}, 95% CI=[{np.percentile(delta_null,2.5):+.3f}, {np.percentile(delta_null,97.5):+.3f}]')

# === Analysis 2: Within Cat 0 ===
print(f'\n=== Within Cat0 ===')
sm0_obs, dm0_obs = within_cat_sm_dm(CAT0, MOA)
delta0_obs = cliffs(sm0_obs, dm0_obs)
print(f'Observed: same n={len(sm0_obs)}, diff n={len(dm0_obs)}, Cliff\'s δ={delta0_obs:+.3f}')
cat0_moas=[MOA[c] for c in CAT0]
delta0_null=[]
for _ in range(N_PERM):
    perm=list(cat0_moas); rng.shuffle(perm)
    moa_shuf={c:m for c,m in zip(CAT0,perm)}
    for c in CAT1: moa_shuf[c]=MOA[c]
    sm,dm=within_cat_sm_dm(CAT0, moa_shuf)
    if len(sm)>0 and len(dm)>0:
        delta0_null.append(cliffs(sm,dm))
delta0_null=np.array(delta0_null)
p0=(delta0_null>=delta0_obs).mean()
print(f'Permutation: P_perm(δ)={p0:.4f}')

# === Analysis 3: Per-MoA coherence (Fig 5e) ===
print(f'\n=== Per-MoA coherence (within-MoA vs between-MoA) ===')
labeled=[c for c in CAT1+CAT0 if MOA.get(c)!=NON_TOX]
tox_moas=sorted(set(MOA[c] for c in labeled))
for m in tox_moas:
    members=[c for c in labeled if MOA[c]==m]
    others=[c for c in labeled if MOA[c]!=m]
    if len(members)<2: continue
    within=[sim.loc[a,b] for a,b in combinations(members,2)]
    between=[sim.loc[a,b] for a in members for b in others]
    delta=cliffs(within, between)
    # permutation: shuffle membership labels (which chemicals belong to this MoA)
    n_in=len(members)
    null_deltas=[]
    for _ in range(N_PERM):
        picked=rng.choice(labeled, size=n_in, replace=False)
        picked_set=set(picked)
        memb_p=[c for c in labeled if c in picked_set]
        oth_p=[c for c in labeled if c not in picked_set]
        wi=[sim.loc[a,b] for a,b in combinations(memb_p,2)]
        bt=[sim.loc[a,b] for a in memb_p for b in oth_p]
        if len(wi)>0: null_deltas.append(cliffs(wi,bt))
    null_deltas=np.array(null_deltas)
    pp=(null_deltas>=delta).mean()
    print(f'  {m:40s} n={n_in}  δ={delta:+.3f}  P_perm={pp:.4f}')

# === Analysis 4: Fig 5b three-way comparison ===
print(f'\n=== Fig 5b: same-MoA vs same-cat/diff-MoA vs diff-cat ===')
cat={c:'1' for c in CAT1}; cat.update({c:'0' for c in CAT0})
groups={'same_moa':[],'same_cat_diff_moa':[],'diff_cat':[]}
for a,b in combinations(labeled,2):
    r=sim.loc[a,b]
    if MOA[a]==MOA[b]: groups['same_moa'].append(r)
    elif cat[a]==cat[b]: groups['same_cat_diff_moa'].append(r)
    else: groups['diff_cat'].append(r)
print(f'  same-MoA n={len(groups["same_moa"])}, same-cat/diff-MoA n={len(groups["same_cat_diff_moa"])}, diff-cat n={len(groups["diff_cat"])}')
d1=cliffs(groups['same_moa'], groups['diff_cat'])
d2=cliffs(groups['same_cat_diff_moa'], groups['diff_cat'])
d3=cliffs(groups['same_moa'], groups['same_cat_diff_moa'])
print(f'  same-MoA vs diff-cat: δ={d1:+.3f}')
print(f'  same-cat/diff-MoA vs diff-cat: δ={d2:+.3f}')
print(f'  same-MoA vs same-cat/diff-MoA: δ={d3:+.3f}')

# Permutation: shuffle MoA labels within each category (preserve category structure)
null_d1=[]; null_d2=[]; null_d3=[]
cat1_m=[MOA[c] for c in CAT1 if MOA[c]!=NON_TOX]
cat1_lab=[c for c in CAT1 if MOA[c]!=NON_TOX]
cat0_m=[MOA[c] for c in CAT0 if MOA[c]!=NON_TOX]
cat0_lab=[c for c in CAT0 if MOA[c]!=NON_TOX]
for _ in range(N_PERM):
    p1=list(cat1_m); rng.shuffle(p1)
    p0=list(cat0_m); rng.shuffle(p0)
    ms={c:m for c,m in zip(cat1_lab,p1)}
    ms.update({c:m for c,m in zip(cat0_lab,p0)})
    g={'same_moa':[],'same_cat_diff_moa':[],'diff_cat':[]}
    for a,b in combinations(cat1_lab+cat0_lab,2):
        r=sim.loc[a,b]
        if ms[a]==ms[b]: g['same_moa'].append(r)
        elif cat[a]==cat[b]: g['same_cat_diff_moa'].append(r)
        else: g['diff_cat'].append(r)
    if len(g['same_moa'])>0:
        null_d1.append(cliffs(g['same_moa'], g['diff_cat']))
        null_d2.append(cliffs(g['same_cat_diff_moa'], g['diff_cat']))
        null_d3.append(cliffs(g['same_moa'], g['same_cat_diff_moa']))
null_d1=np.array(null_d1); null_d2=np.array(null_d2); null_d3=np.array(null_d3)
p1=(null_d1>=d1).mean(); p2=(null_d2>=d2).mean(); p3=(null_d3>=d3).mean()
print(f'  P_perm(same-MoA > diff-cat)={p1:.4f}')
print(f'  P_perm(same-cat/diff-MoA > diff-cat)={p2:.4f}  [under null, labels shuffle within cat, so category-level signal preserved]')
print(f'  P_perm(same-MoA > same-cat/diff-MoA)={p3:.4f}')
