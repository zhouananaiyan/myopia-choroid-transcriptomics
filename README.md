# Choroidal transcriptional programs in experimental myopia

Analysis code, exported inputs and archived results for **Recovery-defined choroidal transcriptional programs show stage- and time-specific directionality in experimental myopia**.

## Scope

The analysis compares recovery-defined choroidal gene modules across experimental myopia datasets. Biological samples or sample pools, not individual cells, are the units of differential-expression inference. The condition-only discovery model is primary because control-recovery pool pairing could not be verified. All six possible pairings are evaluated as sensitivity analyses.

The repository supports rerunning the statistical analysis from the included exported inputs. It is not a complete raw-single-cell processing workflow: the annotated Seurat checkpoints used for export are not included. See [the reproducibility guide](code/README_reproducibility.txt) for inputs, provenance, commands and limitations.

## Contents

- `code/`: R analysis, Python preprocessing and summaries, and optional publication-packaging scripts.
- `analysis/inputs/`: exported pseudobulk matrices, sample metadata, bulk counts, orthology maps and provenance.
- `analysis/results/`: archived differential-expression, CAMERA, enrichment and sensitivity results.
- `analysis/inputs/export_sessionInfo.txt`: recorded R environment for single-cell export.
- `analysis/results/sessionInfo.txt`: recorded R environment for statistical analysis.
- `MANIFEST.csv`: SHA-256 checksums and file sizes for this release, excluding the manifest itself.
- `packaging_build_environment.json`: historical archive-packaging environment, not the statistical-analysis environment.

## Rerun the analysis

Use a separate working copy: the analysis overwrites result files. Install R with edgeR, limma and statmod; consult the recorded session information for the original versions.

From the repository root:

```sh
Rscript code/03_reanalyse.R "/absolute/path/to/repository"
```

To regenerate candidate rankings, summary statistics and figures:

```sh
python -m pip install -r code/requirements-analysis.txt
python code/05_summarize_and_plot.py "/absolute/path/to/repository"
```

`04_enrichment.py` queries the external g:Profiler service and overwrites archived enrichment results. Run it only when intentionally updating that analysis; current database responses may differ from the archived snapshot.

Steps 01 and 02 require external source files and annotated Seurat checkpoints. Do not run them to reproduce results from the supplied exported inputs. Copy `code/config.example.json` to a local configuration only if recreating the export stage.

## Optional packaging

Steps 06-08 prepare supplementary workbooks and archive/PDF outputs. The original Excel writer in step 07 depends on `@oai/artifact-tool`, which is not bundled and is not guaranteed to be publicly installable. This dependency is not required for the statistical analysis. `update_manuscript.py` is an auxiliary document editor requiring a separate manuscript template; it is not an analysis step.

## Data sources

Public GEO datasets: [GSE235684](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235684), [GSE293607](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293607), [GSE227724](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227724), [GSE261232](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261232), and [GSE256200](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256200).

Orthology uses Ensembl Compara release 105. Source URLs, checksums and the GSE256200 tissue-label curation are recorded in `analysis/inputs/mapping_provenance.json`.

## Interpretation

The cross-study comparisons are exploratory. Fixed-correlation CAMERA signals do not retain equivalent corrected support with estimated correlations. Mouse candidate selection and subsequent fibroblast summaries are not independent validation. The guinea-pig modules do not pass the six-test correction across choroid and sclera. Candidate genes are proposed for experimental follow-up, not as validated biomarkers or causal mechanisms.

## Release preparation

This copy was prepared from the submission analysis archive on 20 October 2026. Comment and docstring wording was simplified without changing analysis logic, parameters, inputs or archived results. AI-assistance disclosure and software-dependency provenance were retained. No author identities or software license have been inferred; the authors should choose a license before enabling reuse beyond applicable default terms.
