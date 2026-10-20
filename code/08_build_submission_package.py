from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from datetime import datetime, timezone
import csv
import hashlib
import io
import json
import platform
import sys
from importlib.metadata import PackageNotFoundError, version

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, KeepTogether


ROOT = Path(__file__).resolve().parents[1]
SUBMISSION = ROOT / "submission"
PDF_OUT = SUBMISSION / "Additional_file_8_Supplementary_Methods.pdf"
ZIP_OUT = SUBMISSION / "Additional_file_8_Analysis_Code.zip"


def build_pdf():
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(
        name="MethodsTitle", parent=styles["Title"], fontName="Helvetica-Bold",
        fontSize=16, leading=20, alignment=TA_CENTER, textColor=colors.HexColor("#1F2937"),
        spaceAfter=8 * mm,
    ))
    styles.add(ParagraphStyle(
        name="MethodsHeading", parent=styles["Heading2"], fontName="Helvetica-Bold",
        fontSize=11.5, leading=14, textColor=colors.HexColor("#1F4E78"),
        spaceBefore=4 * mm, spaceAfter=2 * mm,
    ))
    styles.add(ParagraphStyle(
        name="MethodsBody", parent=styles["BodyText"], fontName="Helvetica",
        fontSize=9.2, leading=12.5, textColor=colors.HexColor("#222222"),
        spaceAfter=2.5 * mm,
    ))
    styles.add(ParagraphStyle(
        name="MethodsSmall", parent=styles["BodyText"], fontName="Helvetica",
        fontSize=8, leading=10, textColor=colors.HexColor("#333333"),
        spaceAfter=1.5 * mm,
    ))
    styles.add(ParagraphStyle(
        name="MethodsTableHeader", parent=styles["BodyText"], fontName="Helvetica-Bold",
        fontSize=8, leading=10, textColor=colors.white,
    ))

    def table_cells(rows):
        converted = []
        for row_index, row in enumerate(rows):
            style = styles["MethodsTableHeader"] if row_index == 0 else styles["MethodsSmall"]
            converted.append([Paragraph(str(value), style) for value in row])
        return converted

    story = []
    story.append(Paragraph("Supplementary Methods and Reproducibility Information", styles["MethodsTitle"]))
    story.append(Paragraph(
        "BMC Ophthalmology submission revision package. Reanalysis dated 14 October 2026. "
        "This document describes the reanalysis represented by Additional files 1-7 and the accompanying code package.",
        styles["MethodsBody"],
    ))

    story.append(Paragraph("1. Study design and datasets", styles["MethodsHeading"]))
    story.append(Paragraph(
        "This study is a secondary analysis of five publicly available transcriptomic datasets from the NCBI Gene Expression Omnibus (GEO). "
        "The workflow uses sample-level pseudobulk discovery in chick recovery, mouse cell-type context, independent chick onset and progression contrasts, and exploratory guinea-pig tissue comparisons. "
        "No new participants or animals were recruited for this analysis.", styles["MethodsBody"]))
    dataset_rows = [
        ["Dataset", "Species / material", "Contrast", "Role"],
        ["GSE235684", "Chick choroid scRNA-seq", "Control vs 24 h recovery", "Primary discovery"],
        ["GSE293607", "Mouse posterior eye-cup scRNA-seq", "Saline vs tropicamide", "Cell-type context"],
        ["GSE227724", "Chick choroid bulk RNA-seq", "Occluded vs open eye, onset", "Independent evaluation"],
        ["GSE261232", "Chick choroid bulk RNA-seq", "Occluded vs open eye, progression", "Independent evaluation"],
        ["GSE256200", "Guinea-pig choroid and sclera bulk RNA-seq", "Albino myopic vs hyperopic", "Exploratory tissue context"],
    ]
    table = Table(table_cells(dataset_rows), colWidths=[25*mm, 48*mm, 55*mm, 45*mm], repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1F4E78")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("LEADING", (0, 0), (-1, -1), 10),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#D1D5DB")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F3F6F8")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 5), ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4), ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(table)

    story.append(Paragraph("2. Preprocessing and experimental unit", styles["MethodsHeading"]))
    story.append(Paragraph(
        "Archived single-cell annotations and UMI count matrices were re-read from the local analysis archive. "
        "Primary cell labels and upstream QC were inherited from those checkpoints, not independently reconstructed or re-clustered in this reanalysis. "
        "Stricter fibroblast QC required 500-6,000 detected genes and mitochondrial transcripts below 10%; the chick annotation-score subset required prediction.score.max &gt;=0.9. "
        "Counts were aggregated within cell type and biological sample pool. The biological sample pool, rather than an individual cell, was the experimental unit for differential expression. Cell-level tests were not treated as independent replication.", styles["MethodsBody"]))
    story.append(Paragraph(
        "The primary GSE235684 model was condition-only because the public metadata did not verify the membership pairing between the three control pools and three recovery pools. "
        "All six possible assignments were fitted as sensitivity analyses. Fibroblast annotation filtering, 80% downsampling with seeds 1, 42, and 20261014, and all-confidence one-to-one mapping were retained in the audit outputs. "
        "Broad cell-type testing required six samples with at least 20 cells each, excluding Cycling and Uncertain labels. Genes were filtered with filterByExpr, TMM normalized, and tested with robust edgeR quasi-likelihood models; gene-level FDR used Benjamini-Hochberg correction. "
        "The chick object contained 70,860 cells including 14,374 fibroblasts, and the mouse object contained 65,989 cells including 5,083 fibroblasts. Seven mouse broad types met testing eligibility.", styles["MethodsBody"]))

    story.append(Paragraph("3. Cross-species mapping and recovery modules", styles["MethodsHeading"]))
    story.append(Paragraph(
        "Chick and mouse results were matched using stable Ensembl gene identifiers and Ensembl Compara release 105 ortholog_one2one records. "
        "The primary map required the high-confidence flag. The resulting common tested universe contained 7,335 genes. "
        "At the stated analysis thresholds, the condition-only discovery analysis defined 92 recovery-up genes and 123 recovery-down genes; 31 genes were retained as exploratory opposite-direction candidates. "
        "Chicken-to-guinea-pig mapping used direct one-to-one stable-ID records from the same Compara release, with high-confidence records used for the primary gene-set tests.", styles["MethodsBody"]))

    story.append(Paragraph("4. Bulk validation and gene-set testing", styles["MethodsHeading"]))
    story.append(Paragraph(
        "For GSE227724 and GSE261232, only complete occluded/open-eye pairs were retained at each time point, with bird and eye-condition terms in each model. "
        "The fixed CAMERA inter-gene correlation of 0.01 was used for the primary descriptive analysis, and dataset- and module-specific correlations estimated from the data were used as sensitivity analyses. "
        "Within each sensitivity scenario, Benjamini-Hochberg correction covered three modules, six time points and two chick datasets (36 tests). "
        "For GSE256200, edgeR contrasts used four albino-myopic and four albino-hyperopic samples per tissue. The two internally inconsistent GEO records were curated using the sample title and count-column description: GSM8087498 was assigned to retina/Retina (A995_RR) and GSM8087499 to choroid/Choroid (A995LC); the curation is recorded in mapping_provenance.json and the reproducibility README. "
        "The supplied GSE256200 count archive was validated as a single gzip member; 16,384 trailing non-gzip bytes in the local source copy were removed from the analysis copy and documented in mapping_provenance.json. Module FDR was corrected over six tests across choroid and sclera. These results are exploratory because species, dissection and biological design differ from discovery.", styles["MethodsBody"]))

    story.append(KeepTogether([
        Paragraph("5. Candidate ranking and enrichment", styles["MethodsHeading"]),
        Paragraph(
            "Candidate selection used chick FDR &lt;0.05 and absolute log2 fold change &gt;=0.25, with mouse nominal P &lt;0.05, absolute log2 fold change &gt;=0.10 and opposite signs. "
            "Post hoc ranking used onset opposite-direction FDR counts, retention across the six pairings, onset nominal counts, and then chick discovery FDR. "
            "It is not a validated score and does not establish causality or clinical biomarker performance. "
            "The candidate set was selected in the same mouse fibroblast contrast used in Fig. 2; its gene-distribution FDR is selection-dependent and does not establish independent replication or cell-type specificity. Those exploratory tests use BH correction across 21 cell-type/module comparisons without modeling inter-gene dependence. "
            "Mouse ortholog identifiers were submitted to g:Profiler for GO biological process, molecular function and cellular component, Reactome, KEGG and WikiPathways enrichment using the 7,335 mapped tested genes as the custom background. "
            "The request used the FDR method, a 0.05 threshold and all_results=true. Enrichment was interpreted as functional context rather than mechanistic confirmation; Human Protein Atlas was not queried.", styles["MethodsBody"]),
    ]))

    story.append(Paragraph("6. Reproducibility package and attachment map", styles["MethodsHeading"]))
    attachment_rows = [
        ["Attachment", "Contents"],
        ["Additional file 1.xlsx", "Sample-level cell counts, pseudobulk and discovery/module sensitivity summaries."],
        ["Additional file 2.xlsx", "Candidate list, modules, cross-species tested genes, mapping and bulk effects."],
        ["Additional file 3.xlsx", "GSE227724 time-point DE and CAMERA results."],
        ["Additional file 4.xlsx", "GSE261232 time-point DE and CAMERA results."],
        ["Additional file 5.xlsx", "GSE256200 choroid/sclera DE and CAMERA results."],
        ["Additional file 6.xlsx", "Mouse cell-type, fibroblast subcluster and composition context."],
        ["Additional file 7.xlsx", "All enrichment result rows in a 12-column reader table; full CSV fields and request/response JSON are in the Additional file 8 ZIP."],
        ["Additional file 8 PDF", "This supplementary methods and reproducibility document."],
        ["Additional file 8 ZIP", "Current Python/R/JavaScript code, exported input counts and feature keys, bulk count tables, mappings, full result CSV/JSON, software records and rerun instructions."],
    ]
    attachment_table = Table(table_cells(attachment_rows), colWidths=[42*mm, 131*mm], repeatRows=1)
    attachment_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1F4E78")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#D1D5DB")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F3F6F8")]),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5), ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4), ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(attachment_table)
    story.append(Spacer(1, 3 * mm))
    story.append(Paragraph(
        "The ZIP includes the exported pseudobulk inputs and bulk count tables required by 03_reanalyse.R, so that step can be rerun after installing the recorded R dependencies. "
        "It does not include large single-cell raw matrices or the original annotated Seurat RDS checkpoints. Rebuilding the pseudobulk inputs from raw cells therefore requires those checkpoints, their feature files and the original annotation workflow; it is not an end-to-end raw-data archive. "
        "Public source URLs, run order and limitations are provided in README.txt and code/README_reproducibility.txt.", styles["MethodsSmall"]))
    story.append(Paragraph(
        "The exported-input environment is recorded in analysis/inputs/export_sessionInfo.txt, and the edgeR/limma analysis environment in analysis/results/sessionInfo.txt. "
        "packaging_build_environment.json records only the Python runtime used to assemble the ZIP/PDF, not proof of the original statistical-analysis Python environment. "
        "The archived g:Profiler request/response preserves the reported enrichment snapshot; a new API call may use a newer database release.", styles["MethodsSmall"]))

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor("#6B7280"))
        canvas.drawString(18 * mm, 12 * mm, "Supplementary Methods and Reproducibility Information")
        canvas.drawRightString(192 * mm, 12 * mm, f"Page {doc.page}")
        canvas.restoreState()

    PDF_OUT.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(
        str(PDF_OUT), pagesize=A4, leftMargin=18*mm, rightMargin=18*mm,
        topMargin=17*mm, bottomMargin=18*mm, title="Supplementary Methods and Reproducibility Information",
        author="Study authors",
    )
    doc.build(story, onFirstPage=footer, onLaterPages=footer)


def package_files():
    """Collect code, exported inputs and results, excluding local config and caches."""
    code_suffixes = {".py", ".r", ".mjs", ".js", ".json", ".txt", ".md"}
    include = [p for p in (ROOT / "code").iterdir()
               if p.is_file() and p.suffix.lower() in code_suffixes
               and p.name not in {"config.json"}]
    include.extend(p for p in (ROOT / "analysis" / "results").iterdir()
                   if p.is_file() and p.suffix.lower() in {".csv", ".json", ".txt"})
    # Include the pseudobulk matrices, feature keys and bulk count tables for reruns.
    include.extend(p for p in (ROOT / "analysis" / "inputs").iterdir()
                   if p.is_file() and p.suffix.lower() in {".csv", ".json", ".txt", ".gz"})
    include.append(ROOT / "analysis" / "run_log.txt")
    required = ["code/07_build_workbook.mjs", "code/update_manuscript.py",
                "code/README_reproducibility.txt", "code/config.example.json",
                "code/requirements-analysis.txt", "code/requirements-packaging.txt",
                "analysis/inputs/export_sessionInfo.txt", "analysis/results/sessionInfo.txt",
                "analysis/inputs/chick_feature_key.csv", "analysis/inputs/mouse_feature_key.csv",
                "analysis/inputs/chick_mouse_one2one.csv", "analysis/inputs/chick_gp_one2one.csv",
                "analysis/inputs/GSE227724_rawCounts_dnrs_142samp.txt.gz",
                "analysis/inputs/GSE261232_raw_counts_DNRS_144samp.txt.gz",
                "analysis/inputs/GSE256200_count_table.txt.gz"]
    prefixes = ["chick_alltypes", "chick_fibro_strict", "chick_fibro_confident",
                "chick_fibro_downsample_1", "chick_fibro_downsample_42",
                "chick_fibro_downsample_20261014", "mouse_alltypes",
                "mouse_fibro_strict", "mouse_fibro_clusters"]
    required.extend(f"analysis/inputs/{prefix}_{kind}" for prefix in prefixes
                    for kind in ["counts.csv.gz", "metadata.csv"])
    required.extend(f"analysis/inputs/{gse}_metadata.csv"
                    for gse in ["GSE227724", "GSE261232", "GSE256200"])
    missing = [name for name in required if not (ROOT / name).is_file()]
    if missing:
        raise FileNotFoundError("Required reproducibility inputs missing: " + ", ".join(missing))
    return sorted(set(include))


def build_environment():
    packages = {}
    for name in ["reportlab", "python-docx", "Pillow", "pandas", "numpy", "scipy", "matplotlib", "requests"]:
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = "not installed in packaging runtime"
    return {
        "role": "ZIP/PDF packaging runtime only; not the original statistical-analysis environment",
        "built_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "platform": platform.platform(),
        "installed_packages_at_package_build": packages,
        "R_export_environment": "analysis/inputs/export_sessionInfo.txt",
        "R_analysis_environment": "analysis/results/sessionInfo.txt",
        "limitation": "The original Python analysis environment was not separately frozen. Requirements files list dependencies, not a retrospective version lock.",
    }


def build_zip():
    include = package_files()
    readme = (ROOT / "code" / "README_reproducibility.txt").read_text(encoding="utf-8")
    generated = {
        "README.txt": readme.encode("utf-8"),
        "packaging_build_environment.json": json.dumps(build_environment(), indent=2).encode("utf-8"),
    }
    manifest = io.StringIO(newline="")
    writer = csv.writer(manifest)
    writer.writerow(["archive_path", "bytes", "sha256"])
    for name, content in generated.items():
        writer.writerow([name, len(content), hashlib.sha256(content).hexdigest()])
    for path in include:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        writer.writerow([path.relative_to(ROOT).as_posix(), path.stat().st_size, digest.hexdigest()])
    generated["MANIFEST.csv"] = manifest.getvalue().encode("utf-8")
    ZIP_OUT.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(ZIP_OUT, "w", ZIP_DEFLATED, compresslevel=6) as zf:
        for name, content in generated.items():
            zf.writestr(name, content)
        for path in include:
            rel = path.relative_to(ROOT).as_posix()
            zf.write(path, rel)
    with ZipFile(ZIP_OUT) as zf:
        bad = zf.testzip()
        if bad is not None:
            raise RuntimeError(f"Archive CRC check failed: {bad}")


if __name__ == "__main__":
    build_pdf()
    build_zip()
    print(PDF_OUT)
    print(ZIP_OUT)
