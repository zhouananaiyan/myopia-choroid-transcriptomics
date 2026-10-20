ADDITIONAL FILE 8: ANALYSIS INPUTS, CODE AND RESULTS
Reanalysis date: 2026-10-14

1. Scope and experimental unit

This archive preserves the current analysis code, exported inputs and result
tables. Biological sample pools, not individual cells, are the experimental
units for single-cell pseudobulk differential expression. The GSE235684 pool
pairing could not be verified from public metadata: the condition-only model
is primary, and all six pairing assignments are sensitivity analyses.

This is a rerunnable exported-input archive, not a complete raw-single-cell
processing pipeline. In a fresh copy with the required R packages installed,
03_reanalyse.R can run using the included analysis/inputs files without the
original machine's D: or C: paths. Recreating those exported inputs from raw
single cells is a separate task and requires omitted annotated RDS checkpoints.

2. Contents and precise software-information paths

- code/: current Python (.py), R (.R), JavaScript (.mjs), configuration and
  documentation files; dependency folders and caches are not distributed.
- analysis/inputs/: nine sample-level pseudobulk count matrices; corresponding
  metadata; chick/mouse feature keys; bulk count tables and metadata for all
  three bulk datasets; Compara mapping tables/provenance; archived cell metadata
  and feature information; the available GSE235684 GEO family metadata file.
- analysis/results/: full differential-expression, CAMERA, enrichment,
  candidate, sensitivity, composition and figure-data CSV/JSON outputs.
- analysis/inputs/export_sessionInfo.txt: R environment used for the Seurat
  checkpoint-to-pseudobulk export (Seurat 5.4.0, Matrix 1.7-4).
- analysis/results/sessionInfo.txt: R environment used for the statistical
  reanalysis (R 4.5.2, edgeR 4.8.2, limma 3.66.0 and supporting packages).
- analysis/run_log.txt: archived R run log, ending with REANALYSIS COMPLETE.
- packaging_build_environment.json: the Python environment that assembled
  this ZIP/PDF, NOT the original Python statistical-analysis environment.
- MANIFEST.csv: SHA-256 and byte size for each archived file except the manifest
  itself; this supports verification of the exact current code snapshot.

The nine compressed pseudobulk matrices total about 4.9 MB, the two feature
keys about 2.2 MB, and the three compressed bulk count tables about 9.0 MB.
These are included so the principal R analysis no longer depends on omitted
counts. Large single-cell raw matrices/reads and annotated RDS files are not
included. No user library installation, credentials or dependency directory is
included.

The GSE256200 metadata contains two GEO records whose sample titles/count
descriptions disagree with their source/tissue fields. The supplied metadata
resolves GSM8087498 (retina, A995_RR) as Retina and GSM8087499 (choroid,
A995LC) as Choroid. The curation is recorded in mapping_provenance.json and
is enforced by 01_prepare_mapping.py. The GSE256200 count table is copied as
one validated gzip member; trailing non-gzip bytes are removed and recorded
in the same provenance file.

3. Rerun from the supplied exported inputs (recommended reproducibility path)

Work in a NEW extracted copy: reruns overwrite analysis/results and later
packaging steps overwrite submission artifacts. Preserve the original archive.
Commands below are shell-neutral examples; ROOT means the absolute directory
containing code/ and analysis/. Replace <ROOT> before running the commands.

1) Install R dependencies matching analysis/results/sessionInfo.txt. The main
   analysis requires edgeR, limma and their dependencies (including statmod).
   Seurat is NOT needed when starting from the supplied exported inputs.
2) Run:
   Rscript code/03_reanalyse.R "<ROOT>"
   Optional custom R library: Rscript code/03_reanalyse.R "<ROOT>" "<R_LIBRARY>"
   Alternative logged runner:
   python code/run_analysis.py "<ROOT>" --rscript "<PATH_TO_RSCRIPT>"
   An optional --r-library argument is supported by that runner.
3) To regenerate rankings, summary and figures:
   python code/05_summarize_and_plot.py "<ROOT>"
   Install dependencies listed in code/requirements-analysis.txt first.
   Arial is the plotting preference; font substitution may change appearance,
   but does not alter numeric results.
4) Enrichment is already archived in enrichment_request.json,
   enrichment_response.json and enrichment_complete.csv. Leave these files
   unchanged to retain the reported database snapshot. Only if intentionally
   refreshing the external service, run:
   python code/04_enrichment.py "<ROOT>"
   That command requires network access and may return different results after
   a g:Profiler database update. It overwrites the archived enrichment files.
5) To recreate workbook input specifications from the current result tables:
   python code/06_prepare_workbooks.py "<ROOT>"
   To recreate each Excel workbook with the original writer:
   node code/07_build_workbook.mjs "<ROOT>" 1
   Repeat for workbook numbers 2 through 7. This requires Node.js and the
   @oai/artifact-tool package used by the originating Codex runtime. That package
   is not bundled or claimed to be publicly installable; without it, the CSV
   results and existing XLSX remain available, but exact XLSX regeneration is
   not guaranteed. This packaging dependency is not needed for the R analysis.
6) To rebuild the methods PDF and ZIP from the current directory state:
   python code/08_build_submission_package.py
   This script resolves ROOT relative to its own code/ directory and uses
   reportlab. Build the ZIP LAST, after all code and final manuscript-related
   code edits, so its snapshot includes the final files. The builder deliberately
   fails if required rerun inputs or session records are missing.

update_manuscript.py is an auxiliary manuscript editor, not a statistical
analysis step. It requires the source_clean.docx template and publication
figures; the template, author information and final manuscript are not included
in this analysis ZIP. Do not run it merely to reproduce statistical results.

4. Optional recreation of exported inputs (not self-contained)

The local working tree may contain code/config.json for the original machine,
but that file is deliberately excluded from the submission ZIP. It is NOT a
portable configuration. Start from config.example.json and replace every path.
Do not run steps 01/02 in an extracted archive unless
the external data dependencies below have actually been obtained.

01_prepare_mapping.py --config <YOUR_CONFIG> needs the original feature files,
the source archive's external_validation subdirectories, bulk count files,
sample_metadata.csv files and cached chick Compara file. It downloads the mouse
Compara file and GSE235684 GEO family metadata, but is not a general raw-data
downloader. Typical invocation:
   python code/01_prepare_mapping.py --config "<YOUR_CONFIG>"

02_export_pseudobulk.R needs these original annotated Seurat checkpoints under
source_root (not included and not claimed to be available from GEO):
   GSE235684/GSE235684_allcells_annotated.rds
   GSE293607/GSE293607_annotation_checkpoint.rds
It also needs the original chick/mouse feature files and the recorded Seurat
environment. Typical invocation:
   Rscript code/02_export_pseudobulk.R "<YOUR_CONFIG>"
GEO accessions alone do not reconstruct the original annotation/QC decisions.
The supplied cell metadata, feature keys and pseudobulk matrices preserve the
inputs used here, but are not a substitute for a fully documented raw-cell
annotation workflow. Steps 01/02 overwrite files in analysis/inputs.

5. Public data and reference URLs

GEO series pages (study descriptions and supplementary downloads):
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235684
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293607
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227724
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261232
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256200

Reference sources recorded in analysis/inputs/mapping_provenance.json:
https://ftp.ensembl.org/pub/release-105/tsv/ensembl-compara/homologies/mus_musculus/Compara.105.protein_default.homologies.tsv.gz
https://ftp.ensembl.org/pub/release-105/tsv/ensembl-compara/homologies/gallus_gallus/Compara.105.protein_default.homologies.tsv.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235684/soft/GSE235684_family.soft.gz

g:Profiler endpoint used by 04_enrichment.py:
https://biit.cs.ut.ee/gprofiler/api/gost/profile/

6. Supplementary table and provenance conventions

S1: sample-level cell counts, pseudobulk and discovery/module sensitivity.
S2: all 31 candidates, 215 module genes, 7,335 common-tested genes, stable-ID
    mapping tables and all candidate bulk effects.
S3: GSE227724 onset bulk DE and CAMERA results.
S4: GSE261232 progression bulk DE and CAMERA results.
S5: GSE256200 guinea-pig choroid/sclera bulk DE and CAMERA results.
S6: mouse cell-type context, fibroblast subclusters and composition.
S7: all 9,442 enrichment rows in a 12-column reader table. The unabridged
    18-column enrichment_complete.csv and request/response JSON are in S8 ZIP,
    not inside S7. Query sources were GO:BP, GO:MF, GO:CC, REAC, KEGG and WP;
    Human Protein Atlas was not queried.
S8: supplementary methods PDF and this analysis-input/code/result ZIP.

Workbook README source paths are relative to the archive root. Mapping source
URLs, file checksums and release evidence are in
analysis/inputs/mapping_provenance.json, not the enrichment workbook S7.
The submission ZIP contains no machine-specific config.json. Mapping provenance
uses archive-relative filenames together with public source URLs and checksums;
use the included relative files for the rerun.
Blank values mean missing/not tested/not applicable, not zero. Guinea-pig
FDR_across_chick_stages is not applicable. Candidate rankings are exploratory,
not a validated prediction score; fixed-correlation CAMERA findings did not
retain equivalent corrected support under estimated inter-gene correlation.

7. Environment limitations

requirements-analysis.txt and requirements-packaging.txt describe dependency
roles. They are not retrospective lockfiles. The original Python statistical
environment was not separately frozen, so exact bit-for-bit Python rerun output
is not asserted. The recorded R sessionInfo files apply only to their stated
export and analysis stages. The generated packaging_build_environment.json must
not be cited as the original analysis software environment.
