"""Sensitivity analyses for the MoA structure result (Fig. 5c,e; Supplementary Table 7).

(1) Leave-out: repeat the within-Category 1 same-MoA vs different-MoA comparison after
    removing all chemicals of the antimetabolite and retinoid subgroups (the two subgroups
    with significant per-MoA coherence). Reported in the text as Cliff's delta = -0.08.
(2) Holm correction across the seven per-MoA permutation P values (Supplementary Table 7).

Run from the repository root:  python scripts/Fig5_MoA_leave_out_sensitivity.py
Re-uses the data loading / MoA dictionary / Cliff's delta of Fig5_MoA_permutation_test.py.
"""
import sys, numpy as np
from itertools import combinations
sys.argv = ['x']
exec(open('scripts/Fig5_MoA_permutation_test.py').read().split('rng=np.random')[0])

EXCL = ('Retinoid signalling', 'Antimetabolite DNA synthesis')

def pairs(members, moa):
    same, diff = [], []
    for a, b in combinations(members, 2):
        if moa[a] == NON_TOX or moa[b] == NON_TOX:
            continue
        (same if moa[a] == moa[b] else diff).append(sim.loc[a, b])
    return np.array(same), np.array(diff)

s_all, d_all = pairs(CAT1, MOA)
keep = [c for c in CAT1 if MOA[c] not in EXCL]
s_ex, d_ex = pairs(keep, MOA)
print(f'Within Category 1, all chemicals     : n_same={len(s_all)} n_diff={len(d_all)} delta={cliffs(s_all, d_all):+.3f}')
print(f'Within Category 1, antimetabolite+retinoid chemicals removed ({len(CAT1)-len(keep)} chemicals): '
      f'n_same={len(s_ex)} n_diff={len(d_ex)} delta={cliffs(s_ex, d_ex):+.3f}')

# Holm across the seven per-MoA permutation P values (Supplementary Table 7)
P = {'Antimetabolite DNA synthesis': 0.023, 'Retinoid signalling': 0.025, 'COX inhibition': 0.077,
     'Cereblon IMiD': 0.212, 'Alkylating DNA damage': 0.347, 'Antiepileptic ion channel': 0.766,
     'Kinase inhibitor': 0.453}
items = sorted(P.items(), key=lambda kv: kv[1]); m = len(items); run = 0
print('Holm-adjusted P:')
for i, (k, p) in enumerate(items):
    run = max(run, min(1, (m - i) * p)); print(f'  {k:32s} raw {p:.3f}  Holm {run:.3f}')
