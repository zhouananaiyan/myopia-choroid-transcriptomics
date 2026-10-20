from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import csv
import re

from docx import Document
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_clean.docx"
OUT = ROOT / "submission" / "BMC_Ophthalmology_final_revised_20261014.docx"
FIG_DIR = ROOT / "submission" / "figures"
RESULTS = ROOT / "analysis" / "results"


def set_text(paragraph, text):
    for run in list(paragraph.runs):
        paragraph._p.remove(run._r)
    paragraph.add_run(text)


def top_candidates():
    with (RESULTS / "candidate_prioritization.csv").open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))[:10]
    result = [[
        "Gene", "Chick recovery log2FC", "Chick FDR", "Mouse log2FC",
        "Mouse FDR", "Onset opposite (n/6)", "Rank"
    ]]
    for row in rows:
        result.append([
            row["symbol"],
            f'{float(row["chick_logFC"]):.2f}',
            f'{float(row["chick_FDR"]):.3g}',
            f'{float(row["mouse_logFC"]):.2f}',
            f'{float(row["mouse_FDR"]):.3g}',
            f'{int(row["onset_opposite_times"])}/6',
            row["descriptive_rank"],
        ])
    return result


def update_table(table, rows):
    for i, values in enumerate(rows):
        for cell, value in zip(table.rows[i].cells, values):
            cell.text = str(value)


def replace_media(docx_path):
    temp = docx_path.with_suffix(".tmp.docx")
    with ZipFile(docx_path, "r") as source, ZipFile(temp, "w", ZIP_DEFLATED) as target:
        for item in source.infolist():
            payload = source.read(item.filename)
            if item.filename == "word/media/image1.png":
                payload = (FIG_DIR / "Figure_1.png").read_bytes()
            elif item.filename == "word/media/image2.png":
                payload = (FIG_DIR / "Figure_2.png").read_bytes()
            elif item.filename == "word/media/image3.png":
                payload = (FIG_DIR / "Figure_3.png").read_bytes()
            elif item.filename == "word/media/image4.png":
                payload = (FIG_DIR / "Figure_4.png").read_bytes()
            target.writestr(item, payload)
    temp.replace(docx_path)


def fit_figures_to_source_aspect_ratio(document):
    """Replace figures at their existing width, preserving the new aspect ratio."""
    for index, shape in enumerate(document.inline_shapes, start=1):
        figure = FIG_DIR / f"Figure_{index}.png"
        if not figure.exists():
            continue
        with Image.open(figure) as img:
            width, height = img.size
        shape.height = round(shape.width * height / width)


def clean_final_ooxml(docx_path):
    """Set figure descriptions and remove orphaned review metadata."""
    descriptions = {
        "1": ("Study design and evidence hierarchy across five public transcriptomic datasets.",
              "Two-panel figure: gene counts at FDR < 0.05 for the primary condition-only model and six pool-pairing sensitivities; median mouse log2 fold changes for recovery-down and recovery-up genes. Counts are genes, not animals or samples."),
        "2": ("Exploratory mouse posterior eye-cup cell-type context for the 31-gene opposite-direction candidate set.",
              "Seven-panel cell-type summary of direction-aligned mouse effects for 31 candidate genes. Labels give tested gene counts; the fibroblast comparison is selection-dependent."),
        "3": ("Time-specific CAMERA enrichment during chick myopia onset and established progression.",
              "Four-panel heatmap of signed minus-log10 CAMERA FDR values for fixed and estimated inter-gene correlations across onset and progression time points. The displayed scale spans plus or minus 17."),
        "4": ("Cross-dataset effect-size matrix for the ten highest-ranked descriptive candidates.",
              "Heatmap of log2 fold changes for the ten highest-ranked descriptive candidates across six onset, six progression, and two guinea-pig tissue contrasts; gray indicates missing data."),
    }
    with ZipFile(docx_path, "r") as source:
        members = {item.filename: source.read(item.filename) for item in source.infolist()}
        infos = {item.filename: item for item in source.infolist()}

    document_name = "word/document.xml"
    document = members[document_name].decode("utf-8")
    def replace_docpr(match):
        tag = match.group(0)
        ident_match = re.search(r'\bid="([^"]+)"', tag)
        if ident_match is None or ident_match.group(1) not in descriptions:
            return tag
        title, descr = descriptions[ident_match.group(1)]
        escaped_title = title.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        escaped_descr = descr.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        tag = re.sub(r'\stitle="[^"]*"', "", tag)
        tag = re.sub(r'\sdescr="[^"]*"', "", tag)
        return tag[:-2] + f' descr="{escaped_descr}" title="{escaped_title}"/>'
    document, count = re.subn(r'<wp:docPr\b[^>]*/>', replace_docpr, document)
    if count != 4:
        raise RuntimeError(f"Expected four figure docPr elements, found {count}")
    members[document_name] = document.encode("utf-8")

    core_name = "docProps/core.xml"
    if core_name in members:
        core = members[core_name].decode("utf-8")
        core = re.sub(r"<dc:subject>.*?</dc:subject>", "<dc:subject>BMC Ophthalmology manuscript</dc:subject>", core, count=1)
        core = re.sub(r"<dc:creator>.*?</dc:creator>", "<dc:creator></dc:creator>", core, count=1)
        core = re.sub(r"<dc:description>.*?</dc:description>", "<dc:description>Final revised manuscript</dc:description>", core, count=1)
        core = re.sub(r"<cp:lastModifiedBy>.*?</cp:lastModifiedBy>", "<cp:lastModifiedBy></cp:lastModifiedBy>", core, count=1)
        members[core_name] = core.encode("utf-8")

    rels_name = "word/_rels/document.xml.rels"
    if rels_name in members:
        rels = members[rels_name].decode("utf-8")
        rels = re.sub(r"<Relationship\b[^>]*/>", lambda m: "" if re.search(r'Target="(?:comments(?:Ids|Extensible)?|people)\.xml"', m.group(0)) else m.group(0), rels)
        members[rels_name] = rels.encode("utf-8")

    ct_name = "[Content_Types].xml"
    if ct_name in members:
        content_types = members[ct_name].decode("utf-8")
        content_types = re.sub(r"<Override\b[^>]*/>", lambda m: "" if re.search(r'PartName="/word/(?:comments(?:Ids|Extensible)?|people)\.xml"', m.group(0)) else m.group(0), content_types)
        members[ct_name] = content_types.encode("utf-8")

    for orphan in ("word/comments.xml", "word/commentsIds.xml", "word/commentsExtensible.xml", "word/people.xml"):
        members.pop(orphan, None)
        infos.pop(orphan, None)
    temp = docx_path.with_suffix(".clean.tmp.docx")
    with ZipFile(temp, "w", ZIP_DEFLATED) as target:
        for name, payload in members.items():
            target.writestr(infos.get(name, name), payload)
    temp.replace(docx_path)


doc = Document(SRC)
p = doc.paragraphs
original_references = {i: p[96+i].text.split('. ', 1)[1] for i in range(1,32)}

paragraph_updates = {
    1: "[AUTHORS: insert full names in publication order]",
    2: "[AFFILIATIONS: department, institution, city, country]",
    3: "Corresponding author: [insert full name]; Email: [insert email address]",
    7: "Background: Recovery experiments provide a reference for transcriptional changes during restoration of ocular growth. Whether a choroidal program defined during recovery shows the same direction during myopia onset and established progression remains uncertain.",
    8: "Methods: We reanalysed five public GEO datasets while preserving each study's biological-sample structure. In GSE235684, cell-level counts were aggregated within cell type and pool; a condition-only model was primary because the public files did not verify the control-recovery pool pairing, and all six possible pool assignments were sensitivity analyses. High-confidence one-to-one Ensembl Compara release 105 orthologs were used for chick-mouse matching. GSE293607 provided mouse cell-type context, and recovery modules were evaluated in independent chick onset (GSE227724), established progression (GSE261232), and guinea-pig choroid/sclera (GSE256200) datasets using pseudobulk or bulk models and CAMERA. Fixed and dataset-estimated inter-gene correlations were both assessed.",
    9: "Results: Strict orthology yielded 7,335 common genes, including 92 recovery-up genes, 123 recovery-down genes, and 31 opposite-direction candidates. The six pool assignments produced different numbers of FDR-positive genes, confirming uncertainty in the paired analysis. With a fixed CAMERA correlation of 0.01, selected onset modules showed inverse behavior at 12 h in GSE227724, whereas at ZT20 in GSE261232 the recovery-down module decreased and the recovery-up module increased, matching the recovery directions. After resolving two inconsistent GSE256200 tissue labels, no guinea-pig module retained significance after the six-test correction across choroid and sclera; sclera showed only nominal within-contrast evidence. With estimated correlations, the corresponding stage-level signals did not remain significant after multiple-testing correction. The top descriptive candidates were CLEC3B, PLXDC2, PPFIBP2, SLMAP, and RCAN2.",
    10: "Conclusions: Recovery-defined choroidal programs show time- and stage-dependent directionality, but the strength of statistical evidence depends on pool pairing and gene-correlation assumptions. The candidate genes and module patterns are hypothesis-generating and require prospective, tissue-resolved validation before causal or biomarker claims can be made.",
    14: "Myopia is a common and increasingly prevalent refractive disorder [1]. Large genetic studies implicate multiple biological pathways in susceptibility to myopia [2]. Visual experience is translated into changes in ocular growth through coordinated signaling among the retina, retinal pigment epithelium, choroid, and sclera [3-6]. Understanding how this signaling axis changes across myopia stages is important for interpreting molecular studies of ocular growth.",
    18: "We tested whether a program defined during recovery shows a consistent direction during myopia onset and established progression. The workflow assigned complementary roles to the five datasets: sample-level pseudobulk discovery in chick recovery, cell-type context in mouse, independent stage/time evaluation in chick, and exploratory tissue-level comparison in guinea pig. Because public metadata did not establish the individual control-recovery pool pairing in GSE235684, the primary analysis was condition-only and the six possible pairings were treated as sensitivity analyses.",
    21: "This study was a secondary analysis of de-identified, publicly available transcriptomic datasets from the NCBI Gene Expression Omnibus (GEO) [18]. No new human participants were recruited and no new animal procedures were performed. The five datasets were GSE235684, GSE293607, GSE227724, GSE261232, and GSE256200. Accession-level sample descriptions, processing decisions, and provenance are supplied in Additional files 1-8.",
    24: "Primary single-cell inputs were archived Seurat objects containing retained cells, count matrices, and cell-type annotations [19,20]. This reanalysis inherited the existing annotations and did not independently repeat clustering or validate all upstream quality-control decisions. In GSE235684, fibroblast sensitivity subsets required 500-6,000 detected genes and a mitochondrial proportion <10%, or an annotation prediction score >=0.9. Cell counts are in Additional file 1; the export rules and metadata are in Additional file 8.",
    25: "GSE293607 comprised posterior eye-cup samples containing posterior sclera, choroid, and retina. The archived object supplied primary cell labels and original Seurat cluster assignments. A stricter mouse fibroblast sensitivity subset required 500-6,000 detected genes and a mitochondrial proportion <10%. Broad cell types were tested only when all six samples contributed at least 20 cells; Uncertain and Cycling labels were excluded. Cluster-level analyses required all six samples and at least 10 cells per sample, with clusters below 30 cells per sample treated as exploratory.",
    27: "Raw UMI counts were aggregated within each cell type and biological sample. The biological sample, rather than the cell, served as the experimental unit for differential expression analysis. For GSE235684, edgeR quasi-likelihood models with TMM normalization [21-24] were fitted using a condition-only design as the primary analysis because the public record did not verify the C1-C3 to R1-R3 pool membership. Six permutations of the three control and three recovery labels were fitted as paired sensitivity analyses. Fibroblast annotation filtering and 80% cell downsampling with seeds 1, 42, and 20261014 were also evaluated. Cell-level tests were not used as evidence of independent biological replication.",
    29: "Chick and mouse results were matched through direct high-confidence one-to-one orthologs from Ensembl Compara release 105 [25], using stable Ensembl gene identifiers and retaining ortholog_one2one records with is_high_confidence = 1. Recovery-up and recovery-down modules required chick FDR < 0.05 and absolute log2FC >=0.25 in the common tested universe. Opposite-direction candidates additionally required mouse nominal P < 0.05, absolute log2FC >=0.10, and opposite effect signs; mouse FDR significance was not required. The mapping and common tested genes are in Additional file 2. All-confidence one-to-one mapping and stricter mouse QC were sensitivity analyses.",
    31: "Raw count matrices and GEO metadata for GSE227724 and GSE261232 were obtained from GEO. Analyses were restricted to choroid samples and complete occluded/open-eye pairs at each Zeitgeber time. Models included bird and eye-condition terms within each time point. Discovery-defined modules were tested with limma-voom CAMERA [26,27], using fixed inter-gene correlation 0.01 and, in sensitivity analysis, module-specific estimated correlations. BH correction covered 36 tests (3 modules x 6 times x 2 chick datasets) separately within each scenario. Guinea-pig module tests used a separate six-test correction across the two tissues.",
    33: "For GSE256200, Gallus gallus and Cavia porcellus gene identifiers were mapped directly using high-confidence one-to-one records from the Ensembl Compara release 105 orthology table for the primary analysis; all one-to-one records were a sensitivity analysis. Two GEO sample titles and count descriptions conflicted with their source/tissue fields, so GSM8087498 (retina, A995_RR) was resolved as Retina and GSM8087499 (choroid, A995LC) as Choroid; the curation rule is recorded in Additional file 8. Four albino-myopic and four albino-hyperopic samples were compared separately in choroid and sclera with edgeR. The analysis was exploratory because tissue dissection, species, and biological design differed from the chick discovery dataset. Module-level results were reported under the fixed-correlation and estimated-correlation CAMERA settings; no causal tissue attribution was made.",
    35: "Mouse fibroblasts were stratified by the original Seurat clusters. Counts were aggregated by cluster and biological sample, and cluster composition was compared at the sample level. Cluster-specific differential-expression and composition results were treated as cell-state context and are supplied in Additional file 6. The analysis did not define the primary recovery modules or provide independent biological replication.",
    37: "Mouse ortholog identifiers were submitted to g:Profiler for GO biological process, molecular function and cellular component, Reactome, KEGG, and WikiPathways enrichment [28], using the common tested mouse orthologs as background and the requested g:Profiler FDR adjustment. Additional file 7 contains selected result fields; the complete CSV, request and response JSON are in Additional file 8. Gene-level FDR used Benjamini-Hochberg adjustment [24]. Cross-species Spearman correlation and one-sided Mann-Whitney comparisons of module effects against other common genes were descriptive. Cell-type summaries aligned mouse effects to the direction opposite chick recovery and used exploratory gene-distribution tests with BH adjustment across 21 cell-type/module comparisons; these tests do not account for inter-gene dependence.",
    39: "Sample-level modeling defines a recovery-associated program",
    40: "The archived GSE235684 object contained 70,860 retained cells, including 14,374 cells in the broad fibroblast class. Annotation-score filtering retained 8,531 fibroblasts, and each 80% downsampling analysis retained 11,498 cells. The condition-only primary model identified 347 fibroblast genes at FDR < 0.05, all with absolute log2FC >=0.25. The numbered C1/R1, C2/R2, C3/R3 pairing yielded 665 FDR-positive genes, whereas all six hypothetical pairings ranged from 164 to 665. These sensitivity analyses do not establish the true pool pairing. The 7,335-gene common ortholog universe contained 92 recovery-up and 123 recovery-down genes, with 31 opposite-direction candidates (6 recovery-up and 25 recovery-down).",
    41: "Mouse data provide partial cross-species context",
    42: "After quality control, 65,989 mouse cells remained, and seven broad cell types were represented in all six samples used for cell-type-stratified testing. Genome-wide concordance between chick recovery and mouse saline-versus-tropicamide effects was weak (Spearman rho = -0.0569, P = 1.10 x 10-6). In the mouse contrast, recovery-down genes had a positive median log2 fold change of 0.072 (Mann-Whitney P = 0.0032), whereas recovery-up genes had a median of 0.041 (P = 0.889). This comparison provides cell-type context, but it does not constitute a matched disease-stage replication because the intervention and tissue preparation differ from chick recovery.",
    45: "Cell-type comparisons provide context for candidate genes",
    46: "The 31-gene candidate set was summarized across seven analyzed mouse cell types. The direction-aligned median effect was 0.434 in fibroblasts (exploratory FDR = 2.54 x 10-17) and 0.301 in Muller glia (FDR = 0.0186); the remaining candidate-set comparisons did not pass FDR < 0.05. Because candidates were selected using this same mouse fibroblast contrast, the fibroblast result is selection-dependent, not independent validation. No mouse fibroblast gene passed gene-level FDR < 0.05. The full cell-type/module summaries and gene effects are in Additional file 6.",
    44: "Figure 1. Sensitivity of discovery to pool model and cross-species directional context. Panel A shows genes at FDR < 0.05 for the unpaired primary model and six possible pairings; Panel B shows median mouse log2 fold change for 123 recovery-down and 92 recovery-up genes. Numbers are gene counts, not animal or sample counts.",
    48: "Figure 2. Cell-type context of the 31-gene opposite-direction candidate set in GSE293607. Mouse log2FC is multiplied by the negative sign of chick recovery log2FC, so positive values indicate opposition to recovery. Filled points have exploratory FDR < 0.05 across 21 cell-type/module comparisons. Labels give the tested gene count, not animal or sample counts. Selection used the same fibroblast contrast, so its small FDR does not demonstrate independent replication or cell-type specificity.",
    54: "Figure 3. Competitive gene-set results across chick myopia-onset and established-progression choroidal time series. Cells show signed -log10 CAMERA FDR values. Positive values denote higher expression in occluded eyes and negative values denote lower expression; outlined cells have FDR < 0.05. The displayed color scale is symmetric and spans +/-17 signed -log10(q) in this rendering (values are not clipped at +/-8).",
    50: "Using the fixed recovery modules rather than stage-specific gene lists, we tested whether disease-state changes were opposite to or in the same direction as recovery. GSE227724 included 5-6 complete occluded/open-eye pairs per time point. Under the fixed CAMERA correlation of 0.01, the strongest inverse-direction evidence occurred at 12 h, where recovery-down, recovery-up, and the 31-gene opposite-direction candidate module had stage-level FDR values of 5.57e-6, 2.08e-3, and 3.29e-3, respectively. These values are descriptive because the estimated-correlation sensitivity analysis did not retain significance after stage-level correction.",
    52: "In GSE261232, each time point contained six complete pairs. Under the fixed correlation of 0.01, ZT20 showed strong module-level signals for recovery-down, recovery-up, and opposite-direction candidates, with stage-level FDR values of 2.70e-17, 1.27e-6, and 4.27e-6, respectively. At this time point, recovery-down decreased and recovery-up increased in occluded eyes, matching the recovery directions. The estimated-correlation sensitivity analysis produced a common stage-level FDR of 0.254 for the three modules at ZT20. Therefore, ZT20 is reported as a time-specific pattern under the fixed-correlation model, not as correlation-robust confirmation.",
    56: "Cross-dataset comparison prioritized a descriptive list of 31 opposite-direction candidates after strict orthology and recovery-module filtering. The top ten by the descriptive ranking were CLEC3B, PLXDC2, PPFIBP2, SLMAP, RCAN2, PLCB4, EFEMP1, MT4, DGKH, and SHE. These ten genes were used to illustrate the ranking in Fig. 4; they should not be interpreted as having equivalent evidence across datasets.",
    57: "Table 2. Descriptively ranked candidate genes and cross-dataset evidence.",
    58: "Note: Post hoc ranking orders candidates by the number of onset time points with opposite direction and gene-level FDR < 0.05, retention across six hypothetical pairings, the number of nominally significant opposite-direction onset contrasts, and then chick FDR. The n/6 column counts opposite-direction time points regardless of significance. Mouse selection used nominal P values; none of the ten genes passed mouse gene-level FDR < 0.05. Complete component statistics for all 31 candidates are in Additional file 2.",
    59: "The top-ten effect-size matrix is shown in Fig. 4; the complete 31-gene/component statistics and sensitivity fields are supplied in Additional file 2. The highest-ranked genes are proposed only for future tissue-specific qPCR, RNAscope, or protein-localization studies after independent experimental confirmation.",
    61: "Figure 4. Cross-dataset effect-size matrix for the ten highest-ranked descriptive candidates. The 14 columns show six chick onset time points, six chick progression time points, and guinea-pig choroid and sclera contrasts. Mouse context is summarized in Fig. 2 and Additional file 6 rather than in this matrix. The matrix is a visualization of cross-dataset evidence, not a claim of biological equivalence or causal validation.",
    63: "Five mouse fibroblast clusters were represented in all six samples. Cluster-level differential-expression and composition results were summarized only when all samples were represented and were reported as cell-state context in Additional file 6. Because cluster sizes and annotations are dataset-specific, these analyses were not used to define the primary recovery modules or to support independent replication.",
    65: "After resolving the two inconsistent GSE256200 tissue labels, fixed-correlation CAMERA found no six-test-corrected module-level signal in choroid (q = 0.406-0.476). Sclera showed nominal within-contrast evidence for all three modules (q = 0.031), but none retained significance after the six-test correction across tissues (q = 0.0625). Estimated-correlation results were not significant after correction. The guinea-pig dataset therefore provides, at most, a weak tissue-context hypothesis and does not confirm a choroid- or sclera-specific recovery program.",
    67: "The methodological contribution is a workflow that uses sample-level pseudobulk to define the recovery program, fixes its direction before stage evaluation, maps species through stable-ID one-to-one orthologs, and treats uncertain pool pairing and gene-correlation assumptions as explicit sensitivity dimensions. This design reduces pseudoreplication, but it cannot recover experimental information absent from public metadata.",
    68: "During onset, opposite-direction changes were most apparent at selected time points under the fixed CAMERA correlation, with the clearest module-level pattern at 12 h. Because estimated correlations weakened the corrected evidence, the timing result should be interpreted as a time-window hypothesis rather than a definitive temporal mechanism.",
    69: "At ZT20 during established progression, the fixed-correlation model showed recovery-down decreases and recovery-up increases in occluded eyes. The estimated-correlation analysis attenuated these signals after correction. The result therefore illustrates how the biological interpretation depends on the statistical treatment of correlated genes.",
    70: "The mouse analysis provides cell-type context but requires cautious interpretation. The mouse contrast captures a tropicamide-versus-saline and posterior-eye-cup comparison rather than a direct matched recovery experiment. The strict one-to-one orthology map and sample-level pseudobulk model make the comparison auditable, but the weak or module-specific concordance cannot establish conservation of mechanism.",
    71: "Cross-dataset integration prioritized CLEC3B, PLXDC2, PPFIBP2, SLMAP, and RCAN2 as the highest-ranked descriptive candidates. Their ranking reflects multiple evidence fields and sensitivity retention, not independent replication. Prospective tissue-resolved experiments are required before any of these genes can be proposed as a mechanistic or clinical biomarker.",
    72: "After metadata resolution, the guinea-pig data did not provide six-test-corrected module significance in either tissue. Sclera retained nominal within-contrast evidence, whereas choroid did not; estimated-correlation sensitivity analysis was also non-significant. The guinea-pig dataset is therefore best used to assess tissue-context uncertainty rather than to confirm the chick-derived program.",
    74: "The datasets span species, tissue dissections, disease stages, interventions, and sampling designs, so inference remains cross-study and context-specific. Three pool-level samples per condition in the discovery and mouse datasets limit power. The public GSE235684 metadata did not provide a verifiable C1-C3 to R1-R3 membership map; accordingly, the unpaired model was primary and all six pool assignments were sensitivity analyses.",
    75: "The mouse contrast captures a combined response to tropicamide treatment and posterior-eye-cup preparation. Candidate ranking is therefore exploratory. In addition, CAMERA results depend on the inter-gene correlation assumption: fixed 0.01 produced stronger time-specific signals than dataset-estimated correlations. This difference is reported explicitly and prevents the stage results from being described as correlation-robust validation.",
    76: "Bulk choroidal data do not resolve the cellular origin of the program, and time points sampled from different animals describe cross-sectional states rather than individual trajectories. Guinea-pig tissue comparisons also have different dissection and species context. No new wet-lab experiment was performed, and the study does not support causal, therapeutic, or clinical biomarker claims.",
    78: "This study provides a sample-level framework for testing whether a program defined during recovery changes in the same direction across myopia stages. Under a fixed gene-correlation model, the strict orthology modules showed time-specific directional patterns in chick onset and established progression, while estimated-correlation sensitivity analyses showed that statistical strength was not fully robust. The resulting candidate genes and enrichment signals are hypothesis-generating resources for future tissue-resolved experiments, not evidence of a confirmed pathogenic mechanism.",
    85: "The datasets analysed during the current study are publicly available in the NCBI Gene Expression Omnibus (GEO) under accession numbers GSE235684 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235684), GSE293607 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293607), GSE227724 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227724), GSE261232 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261232), and GSE256200 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256200). Sample-level summaries and supplementary result workbooks are supplied as Additional files 1-7. Full differential-expression outputs, session information, stable-ID mapping tables, and the reproducible analysis scripts are supplied in the Additional file 8 ZIP package.",
    92: "[AUTHOR CHECK REQUIRED: provide author initials and individual contributions. All authors should read and approve the final manuscript.]",
    95: "AI-assisted tools were used for code review, reproducibility checks, and language organization. The authors remain responsible for the analysis, interpretation, and final manuscript text.",
    128: "Additional files",
    129: "Additional file 1: Supplementary Table S1 (.xlsx). Sample-level pseudobulk summaries, cell-type counts, and discovery and module sensitivity summaries. Dataset inventory is in Table 1; export filters are documented in Additional file 8.",
    130: "Additional file 2: Supplementary Table S2 (.xlsx). Complete 31-gene opposite-direction candidate list, component effect sizes, ranking fields, and cross-dataset evidence.",
    131: "Additional file 3: Supplementary Table S3 (.xlsx). GSE227724 time-point differential-expression and CAMERA results under fixed and estimated gene-correlation settings.",
    132: "Additional file 4: Supplementary Table S4 (.xlsx). GSE261232 time-point differential-expression and CAMERA results under fixed and estimated gene-correlation settings.",
    133: "Additional file 5: Supplementary Table S5 (.xlsx). Guinea-pig choroid and sclera differential-expression and module-level results.",
    134: "Additional file 8: Supplementary Methods and analysis code (.pdf and .zip). The PDF gives the reproducibility methods; the ZIP contains the R/Python workflow, software information, session information, mapping rules, full result CSVs, and execution notes.",
}
for index, text in paragraph_updates.items():
    set_text(p[index], text)
set_text(p[28], "Cross-species context using recovery-defined modules")

# Number citations in order of first appearance and omit uncited bibliography entries.
citation_pattern = re.compile(r"\[([0-9,\- ]+)\]")
def citation_numbers(match):
    values = []
    for token in match.group(1).replace(' ', '').split(','):
        if '-' in token:
            first, last = map(int, token.split('-'))
            values.extend(range(first, last+1))
        else:
            values.append(int(token))
    return values
ordered_refs = []
for paragraph in p[:96]:
    for match in citation_pattern.finditer(paragraph.text):
        for ref in citation_numbers(match):
            if ref not in ordered_refs:
                ordered_refs.append(ref)
renumber = {old: new for new, old in enumerate(ordered_refs, 1)}
for paragraph in p[:96]:
    if citation_pattern.search(paragraph.text):
        set_text(paragraph, citation_pattern.sub(
            lambda m: '[' + ','.join(str(n) for n in sorted(renumber[v] for v in citation_numbers(m))) + ']',
            paragraph.text))
for i, old in enumerate(ordered_refs, 97):
    set_text(p[i], f'{renumber[old]}. {original_references[old]}')
for paragraph in p[97+len(ordered_refs):128]:
    paragraph._element.getparent().remove(paragraph._element)
assert set(ordered_refs) <= set(original_references)

# Preserve author placeholders and restore the contributions heading.
set_text(p[92], "Authors' contributions")
for index in (1, 2, 3, 93):
    assert p[index].text == Document(SRC).paragraphs[index].text
for index in (39, 41, 45):
    p[index].style = "Heading 2"
for index in (40, 42, 46):
    p[index].style = "Normal"

inserted_s6 = p[134].insert_paragraph_before("Additional file 6: Supplementary Table S6 (.xlsx). Mouse fibroblast subcluster pseudobulk, composition, and cell-type context results.")
inserted_s6.style = p[133].style
inserted_s7 = p[134].insert_paragraph_before("Additional file 7: Supplementary Table S7 (.xlsx). Functional enrichment results with 12 selected fields for all 9,442 returned records. The complete CSV and g:Profiler request/response JSON are in Additional file 8.")
inserted_s7.style = p[133].style

dataset_rows = [
    ["GSE235684", "Gallus gallus", "Choroid scRNA-seq", "Control versus 24-h recovery", "3 pools/condition; public pairing not verified", "Discovery"],
    ["GSE293607", "Mus musculus", "Posterior eye-cup scRNA-seq", "Tropicamide versus saline", "3 pools/condition", "Cross-species cell-type context"],
    ["GSE227724", "Gallus gallus", "Choroid bulk RNA-seq", "Onset: occluded vs open eye", "5-6 paired birds/time", "Independent evaluation"],
    ["GSE261232", "Gallus gallus", "Choroid bulk RNA-seq", "Progression: occluded vs open eye", "6 paired birds/time", "Independent evaluation"],
    ["GSE256200", "Cavia porcellus", "Choroid and sclera bulk RNA-seq", "Albino myopic vs albino hyperopic", "4 samples/group/tissue", "Exploratory tissue context"],
]
for row, values in zip(doc.tables[0].rows[1:], dataset_rows):
    for cell, value in zip(row.cells, values):
        cell.text = value

update_table(doc.tables[1], top_candidates())
fit_figures_to_source_aspect_ratio(doc)
OUT.parent.mkdir(parents=True, exist_ok=True)
doc.save(OUT)
replace_media(OUT)
clean_final_ooxml(OUT)
print(OUT)
