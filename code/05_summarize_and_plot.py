import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr, mannwhitneyu

p=argparse.ArgumentParser()
p.add_argument('root',type=Path)
a=p.parse_args()
root=a.root
r=root/'analysis/results'
figdir=root/'submission/figures'
figdir.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.titlesize':10,
                     'axes.labelsize':9,'pdf.fonttype':42,'ps.fonttype':42,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'savefig.facecolor':'white'})
blue='#247F9C';red='#CC654C';green='#41856B';ink='#25292C'
cmap=LinearSegmentedColormap.from_list('effects',[blue,'#FFFFFF',red])
cmap.set_bad('#E6E6E6')
def save(fig,name):
    fig.savefig(figdir/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(figdir/f'{name}.png',dpi=600,bbox_inches='tight')
    plt.close(fig)

joint=pd.read_csv(r/'cross_species_common_genes.csv')
candidates=joint[joint.candidate].copy()
bulk=pd.read_csv(r/'candidate_all_bulk_effects.csv')
sen=pd.read_csv(r/'discovery_sensitivity_summary.csv')
mods=pd.read_csv(r/'module_sensitivity_summary.csv')
membership=pd.read_csv(r/'all_scenario_membership.csv')
camera=pd.read_csv(r/'CAMERA_primary.csv')
camall=pd.read_csv(r/'CAMERA_all_scenarios.csv')
rows=[]
paired_names=[x for x in membership.scenario.unique() if x.startswith('paired_')]
for _,gene in candidates.iterrows():
    g=bulk[bulk.chick_gene_id.eq(gene.chick_gene_id)]
    onset=g[g.dataset.eq('GSE227724')]
    progression=g[g.dataset.eq('GSE261232')]
    direction=np.sign(gene.chick_logFC)
    onset_opposite=np.sign(onset.logFC).eq(-direction)
    paired=membership[membership.chick_gene_id.eq(gene.chick_gene_id)&membership.scenario.isin(paired_names)]
    rows.append({**gene.to_dict(),
      'onset_timepoints_tested':len(onset),'onset_opposite_times':int(onset_opposite.sum()),
      'onset_opposite_nominal_times':int((onset_opposite & onset.PValue.lt(.05)).sum()),
      'onset_opposite_FDR_times':int((onset_opposite & onset.FDR.lt(.05)).sum()),
      'progression_FDR_times':int(progression.FDR.lt(.05).sum()),
      'paired_schemes_retaining_candidate':int(paired.candidate.sum()),
      'strict_QC_retained':bool(membership.loc[membership.chick_gene_id.eq(gene.chick_gene_id)&membership.scenario.eq('fibro_strict'),'candidate'].any())})
ranked=pd.DataFrame(rows).sort_values(['onset_opposite_FDR_times','paired_schemes_retaining_candidate',
    'onset_opposite_nominal_times','chick_FDR'],ascending=[False,False,False,True])
ranked['descriptive_rank']=np.arange(1,len(ranked)+1)
ranked.to_csv(r/'candidate_prioritization.csv',index=False)

summary={'common_genes':len(joint),'up':int(joint.module.eq('Recovery_up').sum()),
 'down':int(joint.module.eq('Recovery_down').sum()),'candidates':len(candidates),
 'candidate_up':int(candidates.module.eq('Recovery_up').sum()),
 'candidate_down':int(candidates.module.eq('Recovery_down').sum()),
 'top10':ranked.head(10).symbol.tolist(),
 'rho':float(spearmanr(joint.chick_logFC,joint.mouse_logFC).statistic),
 'spearman_p':float(spearmanr(joint.chick_logFC,joint.mouse_logFC).pvalue)}
for mod,alt in [('Recovery_up','less'),('Recovery_down','greater')]:
    ix=joint.module.eq(mod)
    summary[mod]={'median_mouse_logFC':float(joint.loc[ix,'mouse_logFC'].median()),
      'MWU_p':float(mannwhitneyu(joint.loc[ix,'mouse_logFC'],joint.loc[~ix,'mouse_logFC'],alternative=alt).pvalue)}
key=[]
for scenario in camall.scenario.unique():
    z=camall[camall.scenario.eq(scenario)]
    for dataset,tp in [('GSE227724','12h'),('GSE261232','Z20')]:
        block=z[z.dataset.eq(dataset)&z.contrast.eq(tp)]
        for _,v in block.iterrows():
            key.append({'scenario':scenario,'dataset':dataset,'contrast':tp,'module':v.module,
                        'Direction':v.Direction,'FDR':v.FDR_across_chick_stages})
pd.DataFrame(key).to_csv(r/'key_timepoint_sensitivity.csv',index=False)
(r/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))

# Figure 1: sample-model uncertainty beside descriptive cross-species context.
fig,axs=plt.subplots(1,2,figsize=(9.2,3.8),gridspec_kw={'width_ratios':[1.15,1]},layout='constrained')
z=sen[sen.model.eq('unpaired')|sen.model.str.startswith('paired_')]
labels=['Unpaired (primary)']+[f'Pairing {v[-3:]}' for v in z.model.iloc[1:]]
axs[0].barh(np.arange(len(z)),z.FDR05,color=[green]+['#BBC3C6']*(len(z)-1),height=.6)
axs[0].set_yticks(np.arange(len(z)),labels)
axs[0].invert_yaxis();axs[0].set_xlim(0,z.FDR05.max()*1.17)
for i,n in enumerate(z.FDR05):axs[0].text(n+9,i,str(n),va='center',fontsize=8)
axs[0].set_xlabel('Genes at FDR < 0.05');axs[0].set_title('A  Discovery depends on pool model',loc='left',fontweight='bold')
for i,mod in enumerate(['Recovery_down','Recovery_up']):
    med=summary[mod]['median_mouse_logFC'];n=int(joint.module.eq(mod).sum())
    axs[1].plot([0,med],[i,i],color=red if i==0 else blue,lw=2)
    axs[1].scatter([med],[i],color=red if i==0 else blue,s=55,zorder=3)
    axs[1].text(med+.008,i,f'{med:.3f}',va='center',fontsize=8)
axs[1].set_yticks([0,1],[f'Recovery-down\n{summary["down"]} genes',f'Recovery-up\n{summary["up"]} genes'])
axs[1].set_ylim(1.6,-.6);axs[1].axvline(0,color='#90989D',lw=.8)
axs[1].set_xlim(-.015,max(summary['Recovery_down']['median_mouse_logFC'],summary['Recovery_up']['median_mouse_logFC'])+.075)
axs[1].set_xlabel('Median mouse log2 fold change\nSaline relative to tropicamide')
axs[1].set_title('B  Cross-species direction',loc='left',fontweight='bold')
save(fig,'Figure_1')

ct=pd.read_csv(r/'celltype_directional_tests.csv')
ct=ct[ct.module.eq('Opposite_candidates')].sort_values('median_expected_effect',ascending=False)
fig,ax=plt.subplots(figsize=(7.1,3.7),layout='constrained')
for i,row in enumerate(ct.itertuples()):
    color=green if row.celltype=='Fibroblast' else blue
    ax.plot([0,row.median_expected_effect],[i,i],color=color,lw=1.6)
    ax.scatter(row.median_expected_effect,i,s=54,facecolors=color if row.FDR<.05 else 'white',edgecolors=color,zorder=3)
    ax.text(.99,i,f'{row.n_genes} genes; q = {row.FDR:.2g}',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=8)
ax.set_yticks(np.arange(len(ct)),[s.replace('_',' ') for s in ct.celltype]);ax.invert_yaxis()
ax.axvline(0,color='#A8AFB2',lw=.8);ax.set_xlim(min(-.06,ct.median_expected_effect.min()-.03),ct.median_expected_effect.max()+.23)
ax.set_xlabel('Median direction-aligned log2 fold change')
ax.set_title('Candidate context across mouse cell types',loc='left',fontweight='bold')
save(fig,'Figure_2')

modules=['Recovery_down','Recovery_up','Opposite_candidates']
panels=[]
for scenario in ['unpaired','estimated_gene_correlation']:
    for dataset,times in [('GSE227724',['00h','04h','08h','12h','16h','20h']),('GSE261232',['Z00','Z04','Z08','Z12','Z16','Z20'])]:
        sub=camall[camall.dataset.eq(dataset)&camall.scenario.eq(scenario)]
        values=np.full((3,6),np.nan);qs=values.copy()
        for i,mod in enumerate(modules):
            for j,tp in enumerate(times):
                row=sub[sub.module.eq(mod)&sub.contrast.eq(tp)].iloc[0]
                qs[i,j]=row.FDR_across_chick_stages
                values[i,j]=(-np.log10(max(row.FDR_across_chick_stages,1e-300)))*(1 if row.Direction=='Up' else -1)
        panels.append((values,qs))
limit=float(np.ceil(max(np.nanmax(np.abs(x[0])) for x in panels)))
fig,axs=plt.subplots(2,2,figsize=(9.8,6.3),sharey=True,layout='constrained')
titles=['A  Onset; fixed correlation 0.01','B  Progression; fixed correlation 0.01',
        'C  Onset; estimated correlation','D  Progression; estimated correlation']
for ax,(values,qs),title in zip(axs.flat,panels,titles):
    im=ax.imshow(values,cmap=cmap,vmin=-limit,vmax=limit,aspect='auto')
    ax.set_xticks(range(6),['0','4','8','12','16','20']);ax.set_xlabel('Zeitgeber time (h)')
    ax.set_yticks(range(3),['Recovery-down','Recovery-up',f'{len(candidates)}-gene candidates'])
    ax.set_title(title,loc='left',fontweight='bold')
    for i in range(3):
        for j in range(6):
            ax.text(j,i,f'{values[i,j]:.1f}',ha='center',va='center',fontsize=8,color='white' if abs(values[i,j])>.6*limit else ink)
            if qs[i,j]<.05:
                ax.add_patch(plt.Rectangle((j-.47,i-.47),.94,.94,fill=False,color=ink,lw=1))
    ax.set_xticks(np.arange(-.5,6,1),minor=True);ax.set_yticks(np.arange(-.5,3,1),minor=True)
    ax.grid(which='minor',color='white',lw=1);ax.tick_params(which='minor',bottom=False,left=False)
fig.colorbar(im,ax=axs.ravel().tolist(),shrink=.85,pad=.03,label='Signed -log10(q), 36-test correction')
save(fig,'Figure_3')

top=ranked.head(10)
columns=[('GSE227724',v) for v in ['00h','04h','08h','12h','16h','20h']]+[('GSE261232',v) for v in ['Z00','Z04','Z08','Z12','Z16','Z20']]+[('GSE256200',v) for v in ['Choroid','Sclera']]
values=np.full((len(top),len(columns)),np.nan)
for i,g in enumerate(top.chick_gene_id):
    for j,(dataset,tp) in enumerate(columns):
        sub=bulk[bulk.chick_gene_id.eq(g)&bulk.dataset.eq(dataset)&bulk.contrast.eq(tp)]
        if len(sub): values[i,j]=sub.iloc[0].logFC
effect=pd.DataFrame(values,index=top.symbol,columns=[f'{g}:{t}' for g,t in columns])
effect.to_csv(r/'figure4_selected_effects.csv',index_label='symbol')
limit=max(.2,float(np.ceil(np.nanmax(np.abs(values))*10)/10))
fig,ax=plt.subplots(figsize=(10.1,5.0),layout='constrained')
im=ax.imshow(values,cmap=cmap,vmin=-limit,vmax=limit,aspect='auto')
ax.set_yticks(range(len(top)),top.symbol)
ax.set_xticks(range(14),['0','4','8','12','16','20','0','4','8','12','16','20','Choroid','Sclera'],rotation=45,ha='right')
ax.axvline(5.5,color=ink,lw=1.4);ax.axvline(11.5,color=ink,lw=1.4)
for x,label in [(2.5,'Onset'),(8.5,'Progression'),(12.5,'Guinea pig')]:
    ax.text(x,1.04,label,transform=ax.get_xaxis_transform(),ha='center',fontweight='bold',fontsize=9)
for i in range(values.shape[0]):
    for j in range(values.shape[1]):
        if np.isfinite(values[i,j]):ax.text(j,i,f'{values[i,j]:.2f}',ha='center',va='center',fontsize=6.6,color='white' if abs(values[i,j])>.65*limit else ink)
ax.set_xticks(np.arange(-.5,14,1),minor=True);ax.set_yticks(np.arange(-.5,len(top),1),minor=True)
ax.grid(which='minor',color='white',lw=.7);ax.tick_params(which='minor',bottom=False,left=False)
ax.set_xlabel('Zeitgeber time (h), or tissue');fig.colorbar(im,ax=ax,shrink=.85,pad=.02,label='log2 fold change')
save(fig,'Figure_4')

print('Figures and summary complete',flush=True)
