args <- commandArgs(trailingOnly=TRUE)
cfg <- jsonlite::fromJSON(args[1])
if (!is.null(cfg$r_library)) .libPaths(c(cfg$r_library, .libPaths()))
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
out <- file.path(cfg$output_root, "analysis", "inputs")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
set.seed(cfg$seed)

write_gz <- function(x, name) {
  con <- gzfile(file.path(out, name), "wt")
  on.exit(close(con))
  write.csv(x, con, row.names=FALSE)
}

export_set <- function(object, feature_path, prefix) {
  counts <- GetAssayData(object, assay="RNA", layer="counts")
  meta <- object[[]]
  features <- read.delim(feature_path, header=FALSE, stringsAsFactors=FALSE)
  feature_names <- gsub("_","-",make.unique(features[[2]]),fixed=TRUE)
  hit <- match(rownames(counts), feature_names)
  if(anyNA(hit)) {
    print(head(rownames(counts)[is.na(hit)],30))
    stop("Features in checkpoint could not be mapped to the submitted 10X feature IDs")
  }
  stopifnot(!anyNA(hit))
  key <- data.frame(feature=rownames(counts), gene_id=sub("\\..*$", "", features[[1]][hit]),
                    symbol=features[[2]][hit])
  write.csv(key, file.path(out,paste0(prefix,"_feature_key.csv")),row.names=FALSE)
  if (prefix == "mouse") {
    meta$group <- ifelse(meta$group == "Contracted", "Saline", "Tropicamide")
  }
  meta$barcode <- rownames(meta)
  write_gz(meta,paste0(prefix,"_cell_metadata.csv.gz"))
  aggregate_counts <- function(keep, labels, name) {
    labs <- factor(paste(labels[keep], meta$sample[keep], sep="__"))
    design <- sparse.model.matrix(~0+labs)
    pb <- as.matrix(counts[,keep,drop=FALSE] %*% design)
    colnames(pb) <- levels(labs)
    pb <- rowsum(pb, key$gene_id, reorder=FALSE)
    tab <- data.frame(gene_id=rownames(pb), pb, check.names=FALSE)
    write_gz(tab,paste0(prefix,"_",name,"_counts.csv.gz"))
    parts <- do.call(rbind,strsplit(colnames(pb),"__",fixed=TRUE))
    pmeta <- data.frame(column=colnames(pb), celltype=parts[,1], sample=parts[,2],
                       n_cells=as.integer(table(labs)),stringsAsFactors=FALSE)
    sm <- unique(meta[,c("sample","group")])
    pmeta$group <- sm$group[match(pmeta$sample,sm$sample)]
    write.csv(pmeta,file.path(out,paste0(prefix,"_",name,"_metadata.csv")),row.names=FALSE)
    message(prefix, " ", name, " ", nrow(pb), " genes; ", sum(keep), " cells")
  }
  aggregate_counts(rep(TRUE,nrow(meta)), meta$celltype, "alltypes")
  fibro <- meta$celltype == "Fibroblast"
  stringent <- fibro & meta$nFeature_RNA >= 500 & meta$nFeature_RNA <= 6000 & meta$percent.mt < 10
  aggregate_counts(stringent, meta$celltype, "fibro_strict")
  if (prefix == "chick") {
    confident <- fibro & meta$prediction.score.max >= .9
    aggregate_counts(confident,meta$celltype,"fibro_confident")
    for (s in c(1,42,20261014)) {
      set.seed(s)
      keep <- rep(FALSE,nrow(meta))
      for (sample_id in unique(meta$sample)) {
        ix <- which(fibro & meta$sample == sample_id)
        keep[sample(ix, floor(.8*length(ix)))] <- TRUE
      }
      aggregate_counts(keep,meta$celltype,paste0("fibro_downsample_",s))
    }
  } else {
    aggregate_counts(fibro,as.character(meta$seurat_clusters),"fibro_clusters")
  }
}

message("Loading chick annotated checkpoint")
x <- readRDS(file.path(cfg$source_root,"GSE235684","GSE235684_allcells_annotated.rds"))
export_set(x,cfg$chick_features,"chick")
rm(x); gc()
message("Loading mouse annotated checkpoint")
x <- readRDS(file.path(cfg$source_root,"GSE293607","GSE293607_annotation_checkpoint.rds"))
export_set(x,cfg$mouse_features,"mouse")
writeLines(capture.output(sessionInfo()),file.path(out,"export_sessionInfo.txt"))
message("EXPORT COMPLETE")
