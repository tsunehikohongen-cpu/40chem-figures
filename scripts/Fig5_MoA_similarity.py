#!/usr/bin/env python3
"""Build Fig. 5 — mechanism-of-action (MoA) structure of the transcriptional response
of the 35 Category 1 / Category 0 chemicals (ssGSEA, GO BP 2023).

  a  Hierarchically clustered cosine-similarity matrix of the 35 labelled chemicals
     (colour bars: category and MoA subgroup, hatch patterns distinguish MoA subgroups)
  b  Three-way comparison: same MoA / same category-different MoA / different category
     (one-sided MWU P and chemical-level permutation P)
  c  Within Category 1: same-MoA vs different-MoA pairs
  d  Within Category 0: same-MoA vs different-MoA pairs
  e  Per-MoA transcriptional coherence (Cliff's delta, permutation P)

Run from the repository root.
Input : data/NES_matrix_GO_BP_2023_40chemicals.csv  (2,692 GO BP terms x 40 chemicals; the
        columns are reordered internally to the Table 1 order Cat1 | Cat0 | Unknown)
Output: output/Fig5_MoA_similarity.png (600 dpi) and .pdf

The permutation P values (PERM below) are produced by Fig5_MoA_permutation_test.py
(same values as Supplementary Table 7).
"""
import os
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.metrics.pairwise import cosine_similarity
from itertools import combinations
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Patch

# ---- Table 1 chemical order (fixed; this is the paper's No.1-40 order) ----
CAT1_CHEMS = ['Acitretin','Busulfan','Carbamazepine','Cyclophosphamide','Fluconazole',
    '5-Fluorouracil','Hydroxyurea','Isotretinoin','Methotrexate','Phenytoin',
    'Thalidomide','Topiramate','Tretinoin','Valproic acid','Pomalidomide',
    'Aspirin','Cytarabine','Ibuprofen','Trimethadione','Vismodegib']
CAT0_CHEMS = ['Ibrutinib','Pazopanib','Ribavirin','Tacrolimus','Bosentan',
    'Cisplatin','Dabrafenib','Dasatinib','Imatinib','Rotenone',
    'Acetaminophen','Acetylsalicylic acid','D-Glucitol','L-Ascorbic acid',
    'Saccharin Sodium Salt hydrate']
UNK_CHEMS  = ['Deltamethrin',"2,2',4,4',5-Pentabromodiphenyl ether",
    'Bisphenol A','Chlorpyrifos','imidacloprid']
CATEGORY   = {c:'Cat1' for c in CAT1_CHEMS}
CATEGORY.update({c:'Cat0' for c in CAT0_CHEMS})
CATEGORY.update({c:'Unknown' for c in UNK_CHEMS})
ORDERED    = CAT1_CHEMS + CAT0_CHEMS + UNK_CHEMS   # 40

# ---- Load NES matrix (any column order accepted) ----
np.random.seed(42)   # fixed seed: jitter of the points in panels b-d
nes_raw = pd.read_csv('data/NES_matrix_GO_BP_2023_40chemicals.csv', index_col=0)
# Normalise column names (strip quotes/whitespace) and reorder to Table 1
nes_raw.columns = [c.strip().strip('"') for c in nes_raw.columns]
missing = set(ORDERED) - set(nes_raw.columns)
if missing:
    raise ValueError(f'NES matrix missing chemicals: {missing}')
nes = nes_raw[ORDERED]                                   # terms × 40, fixed order
chems = list(nes.columns)
CAT   = [CATEGORY[c] for c in chems]                     # ['Cat1']*20+['Cat0']*15+['Unknown']*5

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
cat_map=dict(zip(chems, CAT))
NON_TOX='Non toxic reference'

labelled=[c for c in chems if cat_map[c]!='Unknown']
sim_full=pd.DataFrame(cosine_similarity(nes.T.values), index=chems, columns=chems)
sim=sim_full.loc[labelled, labelled]

def cliffs(a,b):
    a=np.asarray(a); b=np.asarray(b)
    if len(a)==0 or len(b)==0: return np.nan
    g=np.sum(a[:,None]>b[None,:]); l=np.sum(a[:,None]<b[None,:])
    return (g-l)/(len(a)*len(b))

# ======== Permutation P values (Fig5_MoA_permutation_test.py; Supplementary Table 7) ========
PERM={'Within Cat1': 0.030, 'Within Cat0': 0.949,
      'same_moa_vs_diff_cat': 0.196, 'same_cat_diff_moa_vs_diff_cat': 0.805, 'same_moa_vs_same_cat_diff_moa': 0.191,
      'Antimetabolite DNA synthesis': 0.023, 'Retinoid signalling': 0.025,
      'COX inhibition': 0.077, 'Cereblon IMiD': 0.212,
      'Alkylating DNA damage': 0.347, 'Kinase inhibitor': 0.453,
      'Antiepileptic ion channel': 0.766}

# ---- Panel a: hierarchical clustering on labelled 35 ----
Z=linkage(sim.values, method='average', metric='euclidean')
order=dendrogram(Z, no_plot=True)['leaves']
chem_order=[labelled[i] for i in order]
sim_ord=sim.loc[chem_order, chem_order]

# ---- Panel b: three-way comparison ----
groups={'sm':[],'scdm':[],'dc':[]}
for a,b in combinations(labelled,2):
    ma,mb=MOA[a],MOA[b]
    if ma==NON_TOX or mb==NON_TOX: continue
    r=sim.loc[a,b]
    if ma==mb: groups['sm'].append(r)
    elif cat_map[a]==cat_map[b]: groups['scdm'].append(r)
    else: groups['dc'].append(r)
# one-sided Mann-Whitney U tests (first group greater than second); computed here, not hard-coded
p_smVdc=mannwhitneyu(groups['sm'], groups['dc'], alternative='greater')[1]
p_scdmVdc=mannwhitneyu(groups['scdm'], groups['dc'], alternative='greater')[1]
p_smVscdm=mannwhitneyu(groups['sm'], groups['scdm'], alternative='greater')[1]

# ---- Panels c, d: within Cat 1 and Cat 0 ----
def within_cat(cat_name):
    members=[c for c in labelled if cat_map[c]==cat_name]
    sm,dm=[],[]
    for a,b in combinations(members,2):
        ma,mb=MOA[a],MOA[b]
        if ma==NON_TOX or mb==NON_TOX: continue
        r=sim.loc[a,b]
        if ma==mb: sm.append(r)
        else: dm.append(r)
    return np.array(sm), np.array(dm)
sm1,dm1=within_cat('Cat1'); sm0,dm0=within_cat('Cat0')
p1=mannwhitneyu(sm1,dm1,alternative='greater')[1]
p0=mannwhitneyu(sm0,dm0,alternative='greater')[1]
d1,d0=cliffs(sm1,dm1), cliffs(sm0,dm0)

# ---- Panel e: per-MoA coherence (use permutation P from PERM dict) ----
toxic_moas=[m for m in sorted(set(MOA.values())) if m!=NON_TOX]
rows=[]
for m in toxic_moas:
    members=[c for c in labelled if MOA[c]==m]
    others=[c for c in labelled if MOA[c]!=m and MOA[c]!=NON_TOX]
    if len(members)<2: continue
    within=[sim.loc[a,b] for a,b in combinations(members,2)]
    between=[sim.loc[a,b] for a in members for b in others]
    rows.append({'MoA':m,'n':len(members),'d':cliffs(within,between), 'p_perm':PERM.get(m,np.nan)})
df=pd.DataFrame(rows).sort_values('d', ascending=False)

# ======== 16-MoA colour + hatch map ========
MOA_LIST=['Retinoid signalling','Alkylating DNA damage','Antimetabolite DNA synthesis','Antiepileptic ion channel',
          'HDAC inhibition antiepileptic','Cereblon IMiD','COX inhibition','Azole antifungal CYP',
          'Hedgehog pathway inhibition','Kinase inhibitor','Antiviral nucleoside analog','Calcineurin immunosuppressant',
          'Endothelin receptor antagonist','Platinum DNA crosslinking','Mitochondrial complex I inhibition','Non toxic reference']
COLS=['#C0392B','#8E44AD','#2E86DE','#5D4037','#E91E63','#27AE60','#616161','#CDDC39',
      '#F39C12','#1ABC9C','#AED581','#B0BEC5','#FFB74D','#C62828','#FF6B6B','#BDBDBD']
HATCH=['','//','xx','..','\\\\','||','++','oo','**','--','OO','//','\\\\','xx','..','']
MOA_COL=dict(zip(MOA_LIST, COLS))
MOA_HATCH=dict(zip(MOA_LIST, HATCH))

CAT_COL={'Cat1':'#C0392B','Cat0':'#2166AC'}

# ======== Build composite figure ========
plt.rcParams.update({'font.family':'sans-serif',
    'font.sans-serif':['Helvetica','Nimbus Sans','Arial'],
    'axes.spines.top':False,'axes.spines.right':False})

fig=plt.figure(figsize=(20, 22))
gs_outer=gridspec.GridSpec(3, 1, height_ratios=[2.2, 1.0, 0.95], hspace=0.60, left=0.07, right=0.97, top=0.96, bottom=0.05)

# ---- Panel a: heatmap + dendrogram + category/MoA bars with more spacing ----
# 7 columns: dendrogram | heatmap | MoA | Cat | chem-label area | SPACER | colorbar
# — panel a shifted left; colorbar stays at the far right with a gap after labels
gs_a=gridspec.GridSpecFromSubplotSpec(1,7, subplot_spec=gs_outer[0],
    width_ratios=[4.0, 28, 0.9, 0.9, 7.0, 1.5, 0.9], wspace=0.05)
ax_dend=fig.add_subplot(gs_a[0]); ax_hm=fig.add_subplot(gs_a[1])
ax_moa=fig.add_subplot(gs_a[2]); ax_cat=fig.add_subplot(gs_a[3])
ax_labs=fig.add_subplot(gs_a[4]); cax=fig.add_subplot(gs_a[6])
# dendrogram
dendrogram(Z, orientation='left', labels=labelled, no_labels=True, ax=ax_dend, color_threshold=0, above_threshold_color='grey', link_color_func=lambda x:'grey')
ax_dend.axis('off')
# heatmap — y-tick labels on RIGHT so dendrogram lines never collide with text
im=ax_hm.imshow(sim_ord.values[::-1], cmap='RdBu_r', aspect='auto', vmin=0.3, vmax=0.9, interpolation='nearest')
ax_hm.set_xticks(np.arange(len(chem_order)))
ax_hm.set_xticklabels(chem_order, rotation=55, ha='right', fontsize=10, fontweight='bold')
ax_hm.set_yticks(np.arange(len(chem_order)))
ax_hm.set_yticklabels([], fontsize=10)  # suppress on left
ax_hm.tick_params(axis='x', pad=4); ax_hm.tick_params(axis='y', length=0)
# MoA bar
ax_moa.set_xlim(0,1); ax_moa.set_ylim(0,len(chem_order))
for i,c in enumerate(chem_order[::-1]):
    m=MOA[c]
    ax_moa.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=MOA_COL.get(m,'#aaa'), edgecolor='white',
                                    hatch=MOA_HATCH.get(m,''), linewidth=0.3))
ax_moa.axis('off')
ax_moa.text(0.5, len(chem_order)+0.8, 'MoA', ha='center', fontsize=12, fontweight='bold')
# Category bar — close to MoA (wspace minimal)
ax_cat.set_xlim(0,1); ax_cat.set_ylim(0,len(chem_order))
for i,c in enumerate(chem_order[::-1]):
    ax_cat.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=CAT_COL[cat_map[c]], edgecolor='white', linewidth=0.3))
ax_cat.axis('off')
ax_cat.text(0.5, len(chem_order)+0.8, 'Cat', ha='center', fontsize=12, fontweight='bold')
# Chemical labels — their own axis to the RIGHT of Cat, LEFT of colorbar
ax_labs.set_xlim(0,1); ax_labs.set_ylim(0,len(chem_order)); ax_labs.axis('off')
for i,c in enumerate(chem_order[::-1]):
    ax_labs.text(0.02, i+0.5, c, va='center', ha='left', fontsize=11, fontweight='bold')
# colorbar (narrow strip at right)
cax.set_aspect(20)
cbar=fig.colorbar(im, cax=cax)
cbar.set_label('Cosine similarity', fontsize=12, fontweight='bold')
cbar.ax.tick_params(labelsize=11)
ax_hm.text(-0.16, 1.03, 'a', transform=ax_hm.transAxes, fontsize=24, fontweight='bold', va='bottom', ha='left')
ax_hm.set_title('Transcriptomic (GSEA-NES) similarity among the 35 labelled chemicals', fontsize=14, fontweight='bold', pad=14)

# ---- Panels b, c, d, e in middle+bottom ----
gs_mid=gridspec.GridSpecFromSubplotSpec(1,3, subplot_spec=gs_outer[1], width_ratios=[1.4, 1.0, 1.0], wspace=0.35)
ax_b=fig.add_subplot(gs_mid[0])
ax_c=fig.add_subplot(gs_mid[1])
ax_d=fig.add_subplot(gs_mid[2])
# panel b: three-way
parts=ax_b.violinplot([groups['sm'], groups['scdm'], groups['dc']], positions=[0,1,2], showextrema=False, widths=0.7)
for pc,c in zip(parts['bodies'], ['#2ECC71','#3498DB','#E74C3C']):
    pc.set_facecolor(c); pc.set_alpha(0.5); pc.set_edgecolor('none')
for i,data in enumerate([groups['sm'], groups['scdm'], groups['dc']]):
    xs=np.full(len(data), i) + np.random.uniform(-0.08, 0.08, len(data))
    ax_b.scatter(xs, data, s=8, color='black', alpha=0.5, zorder=3)
    ax_b.plot([i-0.3, i+0.3], [np.median(data)]*2, 'k-', lw=2)
ax_b.set_xticks([0,1,2])
ax_b.set_xticklabels([f'Same MoA\n(n={len(groups["sm"])})',
                     f'Same category,\ndifferent MoA\n(n={len(groups["scdm"])})',
                     f'Different\ncategory\n(n={len(groups["dc"])})'], fontsize=13)
ax_b.set_ylabel('Cosine similarity (GSEA-NES)', fontsize=14, fontweight='bold')
ax_b.tick_params(axis='y', labelsize=13)
# bracket + P values with both MWU (computed) and perm
def bracket(ax, x1, x2, y, text):
    ax.plot([x1,x1,x2,x2],[y, y+0.015, y+0.015, y], color='black', lw=1.0)
    ax.text((x1+x2)/2, y+0.028, text, ha='center', fontsize=12)
ymax=max(max(groups['sm']), max(groups['scdm']), max(groups['dc']))
bracket(ax_b, 0, 2, ymax+0.03, f'MWU P={p_smVdc:.3f}; Pperm={PERM["same_moa_vs_diff_cat"]:.3f}')
bracket(ax_b, 1, 2, ymax+0.12, f'MWU P={p_scdmVdc:.3f}; Pperm={PERM["same_cat_diff_moa_vs_diff_cat"]:.3f}')
bracket(ax_b, 0, 1, ymax+0.21, f'MWU P={p_smVscdm:.3f}; Pperm={PERM["same_moa_vs_same_cat_diff_moa"]:.3f}')
ax_b.set_ylim(min(min(groups['sm']), min(groups['scdm']), min(groups['dc']))-0.05, ymax+0.33)
ax_b.set_title('b  Category vs MoA (three-way)', fontsize=16, fontweight='bold', loc='left', pad=10)

# panel c: within Cat1
parts=ax_c.violinplot([sm1, dm1], positions=[0,1], showextrema=False, widths=0.7)
for pc,c in zip(parts['bodies'], ['#C0392B','#9aa0a8']):
    pc.set_facecolor(c); pc.set_alpha(0.6); pc.set_edgecolor('none')
for i,data in enumerate([sm1,dm1]):
    xs=np.full(len(data), i) + np.random.uniform(-0.08, 0.08, len(data))
    ax_c.scatter(xs, data, s=12, color='black', zorder=3)
    ax_c.plot([i-0.3,i+0.3],[np.median(data)]*2,'k-',lw=2)
ax_c.set_xticks([0,1]); ax_c.set_xticklabels(['Same MoA','Different MoA'], fontsize=13)
ax_c.set_ylabel('Cosine similarity', fontsize=14, fontweight='bold')
ax_c.tick_params(axis='y', labelsize=13)
ax_c.set_title(f'c  Within Category 1\nδ={d1:+.3f}  MWU P={p1:.3f}  Pperm={PERM["Within Cat1"]:.3f}', fontsize=15, fontweight='bold', loc='left', pad=10)

# panel d: within Cat0
parts=ax_d.violinplot([sm0, dm0], positions=[0,1], showextrema=False, widths=0.7)
for pc,c in zip(parts['bodies'], ['#2166AC','#9aa0a8']):
    pc.set_facecolor(c); pc.set_alpha(0.6); pc.set_edgecolor('none')
for i,data in enumerate([sm0,dm0]):
    xs=np.full(len(data), i) + np.random.uniform(-0.08, 0.08, len(data))
    ax_d.scatter(xs, data, s=12, color='black', zorder=3)
    ax_d.plot([i-0.3,i+0.3],[np.median(data)]*2,'k-',lw=2)
ax_d.set_xticks([0,1]); ax_d.set_xticklabels(['Same MoA','Different MoA'], fontsize=13)
ax_d.set_ylabel('Cosine similarity', fontsize=14, fontweight='bold')
ax_d.tick_params(axis='y', labelsize=13)
ax_d.set_title(f'd  Within Category 0\nδ={d0:+.3f}  MWU P={p0:.3f}  Pperm={PERM["Within Cat0"]:.3f}', fontsize=15, fontweight='bold', loc='left', pad=10)

# ---- Panel e: per-MoA coherence with permutation P ----
gs_bot=gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec=gs_outer[2], width_ratios=[2.5, 1.5], wspace=0.15)
ax_e=fig.add_subplot(gs_bot[0])
colors=['#C0392B' if (dv>0.3 and pv<0.1) else '#f39c12' if dv>0.3 else '#95a5a6' if abs(dv)<0.1 else '#2166AC' for dv,pv in zip(df['d'], df['p_perm'])]
y=np.arange(len(df))[::-1]
bars=ax_e.barh(y, df['d'], color=colors, edgecolor='black', linewidth=0.6)
ax_e.set_yticks(y); ax_e.set_yticklabels([f'{m} (n={n})' for m,n in zip(df['MoA'], df['n'])], fontsize=13)
ax_e.set_xlabel("Cliff's δ (within-MoA vs between-MoA)", fontsize=14, fontweight='bold')
ax_e.tick_params(axis='x', labelsize=13)
ax_e.axvline(0, color='black', lw=0.8)
for bar,p in zip(bars, df['p_perm']):
    xpos=bar.get_width()
    # Negative bars: put label on the POSITIVE side of 0 to avoid y-axis tick-label collision
    if xpos>=0:
        ax_e.text(xpos+0.025, bar.get_y()+bar.get_height()/2, f'Pperm={p:.3f}', fontsize=13, va='center', ha='left')
    else:
        ax_e.text(0.025, bar.get_y()+bar.get_height()/2, f'Pperm={p:.3f}', fontsize=13, va='center', ha='left')
ax_e.set_xlim(min(df['d'].min()-0.15, -0.35), max(df['d'].max()+0.30, 1.05))
ax_e.set_title('e  Per-MoA transcriptional coherence (chemical-level permutation)', fontsize=16, fontweight='bold', loc='left', pad=10)

# ---- Legend for MoA (right side of bottom row) ----
ax_leg=fig.add_subplot(gs_bot[1]); ax_leg.axis('off')
handles=[Patch(facecolor=MOA_COL[m], hatch=MOA_HATCH[m], edgecolor='black', linewidth=0.3, label=m) for m in MOA_LIST]
ax_leg.legend(handles=handles, loc='center', fontsize=12, ncol=2, title='MoA subgroup', title_fontsize=14, frameon=False, handlelength=2.0, handleheight=1.3)

os.makedirs('output', exist_ok=True)
for ext in ('png','pdf'):
    fig.savefig(f'output/Fig5_MoA_similarity.{ext}', dpi=600 if ext=='png' else None, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('saved output/Fig5_MoA_similarity.png / .pdf')
print(f'3-way: smVdc MWU={p_smVdc:.4f} Pperm={PERM["same_moa_vs_diff_cat"]:.3f}; scdmVdc {p_scdmVdc:.4f}/{PERM["same_cat_diff_moa_vs_diff_cat"]:.3f}; smVscdm {p_smVscdm:.4f}/{PERM["same_moa_vs_same_cat_diff_moa"]:.3f}')
