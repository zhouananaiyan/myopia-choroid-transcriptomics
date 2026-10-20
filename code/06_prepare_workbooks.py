import json
import sys
from pathlib import Path
import pandas as pd

root=Path(sys.argv[1]);r=root/'analysis/results';inp=root/'analysis/inputs'
out=root/'analysis/workbook_data';out.mkdir(parents=True,exist_ok=True)
sets={i:[] for i in range(1,8)}
def add(n,name,file,frame=None):
    df=pd.read_csv(file) if frame is None else frame
    df=df.drop(columns=[c for c in df.columns if c.startswith('Unnamed:')],errors='ignore')
    record={'name':name[:31],'source':str(Path(file).relative_to(root)).replace('\\','/'),
            'columns':df.columns.tolist(),
            'rows':df.astype(object).where(pd.notna(df),None).values.tolist()}
    sets[n].append(record)

add(1,'Pseudobulk_summary',r/'pseudobulk_summary.csv')
add(1,'Discovery_sensitivity',r/'discovery_sensitivity_summary.csv')
add(1,'Module_sensitivity',r/'module_sensitivity_summary.csv')
for species in ['chick','mouse']:
    add(1,f'{species}_sample_cell_counts',inp/f'{species}_alltypes_metadata.csv')

add(2,'Candidate_prioritization',r/'candidate_prioritization.csv')
add(2,'Recovery_modules',r/'recovery_modules.csv')
add(2,'Common_tested_genes',r/'cross_species_common_genes.csv')
add(2,'Candidate_bulk_effects',r/'candidate_all_bulk_effects.csv')
add(2,'Chicken_mouse_one2one',inp/'chick_mouse_one2one.csv')
add(2,'Chicken_GP_one2one',inp/'chick_gp_one2one.csv')
# Keep the full scenario-membership table in analysis/results, outside the workbook.

cam=pd.read_csv(r/'CAMERA_all_scenarios.csv')
for n,gse in [(3,'GSE227724'),(4,'GSE261232'),(5,'GSE256200')]:
    add(n,'CAMERA_primary',r/'CAMERA_all_scenarios.csv',cam[cam.dataset.eq(gse)&cam.scenario.eq('unpaired')])
    add(n,'CAMERA_sensitivity',r/'CAMERA_all_scenarios.csv',cam[cam.dataset.eq(gse)])
    add(n,'Sample_metadata',inp/f'{gse}_metadata.csv')
    for f in sorted(r.glob(f'{gse}_*_DE.csv')):add(n,f.stem.replace(gse+'_',''),f)
en=pd.read_csv(r/'enrichment_complete.csv')
cols=['query','source','native','name','p_value','significant','term_size','query_size',
      'intersection_size','effective_domain_size','precision','recall']
add(7,'Functional_enrichment',r/'enrichment_complete.csv',en[[c for c in cols if c in en]])
add(6,'Celltype_tests',r/'celltype_directional_tests.csv')
add(6,'Celltype_gene_effects',r/'celltype_gene_effects.csv')
add(6,'Subcluster_summary',r/'mouse_subcluster_summary.csv')
add(6,'Composition_tests',r/'mouse_subcluster_composition.csv')
add(6,'Subcluster_cell_counts',inp/'mouse_fibro_clusters_metadata.csv')
for f in sorted(r.glob('mouse_cluster_*_DE.csv')):add(6,f.stem.replace('mouse_',''),f)

notes=[['Field or rule','Definition'],
 ['Version','Reanalysis 2026-10-14. Numerical cells retain source precision within Excel numeric precision. Authoritative CSV and JSON results are in Additional file 8. No hidden rows.'],
 ['Experimental unit','Biological sample pool for single-cell pseudobulk; paired bird within each chick bulk time point. Cells are not independent replicates.'],
 ['Discovery primary','Condition-only model; uncertainty of pool matching is assessed with all six numbered pairing assignments.'],
 ['logFC / log2FC','Recovery vs control (chick discovery); saline vs tropicamide (mouse); occluded vs open eye (chick bulk); albino myopic vs albino hyperopic (guinea pig).'],
 ['PValue','Unadjusted probability from the named test.'],
 ['FDR in DE tables','BH correction over all tested genes in that cell type or tissue/time contrast.'],
 ['FDR_across_chick_stages','BH over 3 sets x 6 time points x 2 chick datasets, separately within each sensitivity scenario. Primary manuscript uses this 36-test correction.'],
 ['FDR_across_time','BH over 18 tests per chick dataset, or six tests across the two guinea-pig tissues.'],
 ['FDR_within_contrast','BH over the three gene sets at one tissue/time contrast, retained for transparency.'],
 ['Correlation','CAMERA inter-gene correlation: fixed 0.01 except estimated_gene_correlation, where estimated separately for each set.'],
 ['Candidate set','Exploratory selection using chick FDR<0.05, abs(logFC)>=0.25; mouse nominal P<0.05, abs(logFC)>=0.10 and opposite signs.'],
 ['Descriptive rank','Post hoc ordering by onset FDR-supported opposite-direction counts, six-pairing retention, onset nominal counts, then discovery FDR. Not a validated clinical score.'],
 ['Blank cell','Missing/not tested/not applicable; not zero.'],
 ['ortholog_one2one','Direct Ensembl Compara release 105 gene-ID mapping. Main analysis requires high confidence; all confidence records are retained for audit.'],
 ['Enrichment p_value','Already adjusted using the requested g:Profiler FDR method, within each query; exploratory mouse ortholog enrichment with a custom common-tested background.'],
 ['Celltype directional tests','Exploratory gene-distribution tests; genes are correlated and the fibroblast candidate set was selected in the same mouse comparison.'],
 ['Cell labels','Inherited from archived annotations; QC/subsampling changes retain these labels. Not a new independent annotation validation.'],
 ['GEO sources','https://www.ncbi.nlm.nih.gov/geo/ ; accessions GSE235684, GSE293607, GSE227724, GSE261232, GSE256200.'],
 ['Provenance','Source paths refer to the analysis archive. Orthology tables are in Additional file 2; full source URLs, g:Profiler request/response JSON, and derivation code are in Additional file 8.']]
for n,sheets in sets.items():
    names=[s['name'] for s in sheets]
    assert len(set(names))==len(names),names
    sheets.append({'name':'README','source':'Analysis protocol and data dictionary','columns':notes[0],
                   'rows':notes[1:]+[[s['name'],s['source']] for s in sheets]})
    (out/f'Additional_file_{n}.json').write_text(json.dumps(sheets,ensure_ascii=False),encoding='utf-8')
    print(n,len(sheets),'sheets',sum(len(s['rows']) for s in sheets),'rows')
