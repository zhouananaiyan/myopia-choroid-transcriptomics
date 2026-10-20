args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1],winslash="/",mustWork=TRUE)
if(length(args)>=2) .libPaths(c(args[2],.libPaths()))
if (nzchar(Sys.getenv("REVISION_R_LIBRARY"))) .libPaths(c(Sys.getenv("REVISION_R_LIBRARY"),.libPaths()))
suppressPackageStartupMessages({library(edgeR);library(limma)})
inputs <- file.path(root,"analysis","inputs")
out <- file.path(root,"analysis","results")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(20261014)
write_out <- function(x,name) write.csv(x,file.path(out,name),row.names=FALSE)
read_counts <- function(name) {
  z <- read.csv(gzfile(file.path(inputs,name)),check.names=FALSE)
  x <- as.matrix(z[,-1,drop=FALSE]);storage.mode(x)<-"double"
  rownames(x)<-sub("\\..*$","",z[[1]])
  if (anyDuplicated(rownames(x))) x<-rowsum(x,rownames(x),reorder=FALSE)
  stopifnot(all(is.finite(x)),all(x>=0),all(x==floor(x)))
  x
}
fit_de <- function(x,design,coef) {
  stopifnot(qr(design)$rank==ncol(design))
  y <- DGEList(x)
  y <- y[filterByExpr(y,design=design),,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y,design,robust=TRUE)
  fit <- glmQLFit(y,design,robust=TRUE)
  de <- topTags(glmQLFTest(fit,coef=coef),n=Inf,sort.by="PValue")$table
  de$gene_id <- rownames(de)
  list(y=y,de=de,design=design,coef=match(coef,colnames(design)))
}

keys <- lapply(c("chick","mouse"),function(p) read.csv(file.path(inputs,paste0(p,"_feature_key.csv"))))
names(keys)<-c("chick","mouse")
add_symbols <- function(de,species) {
  de$symbol <- keys[[species]]$symbol[match(de$gene_id,keys[[species]]$gene_id)]
  de
}
chick_x<-read_counts("chick_alltypes_counts.csv.gz")
chick_m<-read.csv(file.path(inputs,"chick_alltypes_metadata.csv"))
mouse_x<-read_counts("mouse_alltypes_counts.csv.gz")
mouse_m<-read.csv(file.path(inputs,"mouse_alltypes_metadata.csv"))
pb_models<-list();pb_summary<-list();mouse_fits<-list()
for (species in c("chick","mouse")) {
  xall<-if(species=="chick") chick_x else mouse_x
  mall<-if(species=="chick") chick_m else mouse_m
  for (ct in setdiff(unique(mall$celltype),c("Uncertain","Cycling"))) {
    m<-mall[mall$celltype==ct,]
    if(nrow(m)!=6||any(m$n_cells<20)) next
    m$group<-relevel(factor(m$group),ref=if(species=="chick") "Control" else "Tropicamide")
    design<-model.matrix(~group,m)
    coef<-colnames(design)[2]
    fit<-fit_de(xall[,m$column,drop=FALSE],design,coef)
    fit$de<-add_symbols(fit$de,species)
    write_out(fit$de,paste0(species,"_",ct,"_unpaired_DE.csv"))
    pb_summary[[paste(species,ct)]]<-data.frame(species=species,celltype=ct,
      samples=nrow(m),min_cells=min(m$n_cells),genes=nrow(fit$de),FDR05=sum(fit$de$FDR<.05),
      FDR05_effect025=sum(fit$de$FDR<.05&abs(fit$de$logFC)>=.25))
    if(species=="chick"&&ct=="Fibroblast") pb_models$unpaired<-fit$de
    if(species=="mouse") mouse_fits[[ct]]<-fit
  }
}
write_out(do.call(rbind,pb_summary),"pseudobulk_summary.csv")

# Test all six pool pairings; do not select a pairing by model fit.
perms<-rbind(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
m<-chick_m[chick_m$celltype=="Fibroblast",]
m<-m[match(c("C1","C2","C3","R1","R2","R3"),m$sample),]
m$group<-relevel(factor(m$group),ref="Control")
for(i in seq_len(nrow(perms))) {
  m$pair<-factor(c(1:3,perms[i,]))
  d<-model.matrix(~pair+group,m)
  fit<-fit_de(chick_x[,m$column,drop=FALSE],d,"groupRecovery")
  nm<-paste0("paired_",paste(perms[i,],collapse=""))
  pb_models[[nm]]<-add_symbols(fit$de,"chick")
  write_out(pb_models[[nm]],paste0("sensitivity_",nm,"_DE.csv"))
}
for(nm in c("fibro_strict","fibro_confident","fibro_downsample_1","fibro_downsample_42","fibro_downsample_20261014")) {
  x<-read_counts(paste0("chick_",nm,"_counts.csv.gz"))
  m<-read.csv(file.path(inputs,paste0("chick_",nm,"_metadata.csv")))
  m$group<-relevel(factor(m$group),ref="Control")
  d<-model.matrix(~group,m)
  fit<-fit_de(x[,m$column,drop=FALSE],d,"groupRecovery")
  pb_models[[nm]]<-add_symbols(fit$de,"chick")
  write_out(pb_models[[nm]],paste0("sensitivity_",nm,"_DE.csv"))
}
sen<-do.call(rbind,lapply(names(pb_models),function(nm) {
  z<-merge(pb_models$unpaired,pb_models[[nm]],by="gene_id",suffixes=c("_main","_sensitivity"))
  a<-pb_models$unpaired$gene_id[pb_models$unpaired$FDR<.05&abs(pb_models$unpaired$logFC)>=.25]
  b<-pb_models[[nm]]$gene_id[pb_models[[nm]]$FDR<.05&abs(pb_models[[nm]]$logFC)>=.25]
  data.frame(model=nm,genes_tested=nrow(pb_models[[nm]]),FDR05=sum(pb_models[[nm]]$FDR<.05),
    FDR05_effect025=length(b),logFC_correlation=cor(z$logFC_main,z$logFC_sensitivity),
    sign_agreement=mean(sign(z$logFC_main)==sign(z$logFC_sensitivity)),
    significant_overlap=length(intersect(a,b)),jaccard=length(intersect(a,b))/length(union(a,b)))
}))
write_out(sen,"discovery_sensitivity_summary.csv")
message("Pseudobulk and pairing sensitivity complete")

map<-read.csv(file.path(inputs,"chick_mouse_one2one.csv"))
map_high<-map[!is.na(map$is_high_confidence)&map$is_high_confidence==1,]
gpmap<-read.csv(file.path(inputs,"chick_gp_one2one.csv"))
make_sets<-function(chick_de,mp,mouse_de) {
  c<-chick_de[,c("gene_id","symbol","logFC","PValue","FDR")]
  colnames(c)<-c("chick_gene_id","symbol","chick_logFC","chick_PValue","chick_FDR")
  mm<-mouse_de[,c("gene_id","symbol","logFC","PValue","FDR")]
  colnames(mm)<-c("mouse_gene_id","mouse_symbol","mouse_logFC","mouse_PValue","mouse_FDR")
  j<-merge(merge(mp,c,by="chick_gene_id"),mm,by="mouse_gene_id")
  stopifnot(!anyDuplicated(j$chick_gene_id),!anyDuplicated(j$mouse_gene_id))
  j$module<-ifelse(j$chick_FDR<.05&j$chick_logFC>=.25,"Recovery_up",
                   ifelse(j$chick_FDR<.05&j$chick_logFC<=-.25,"Recovery_down","Background"))
  j$candidate<-j$module!="Background"&j$mouse_PValue<.05&abs(j$mouse_logFC)>=.1&
    sign(j$mouse_logFC)==-sign(j$chick_logFC)
  j
}
joint<-make_sets(pb_models$unpaired,map_high,mouse_fits$Fibroblast$de)
write_out(joint,"cross_species_common_genes.csv")
write_out(joint[joint$module!="Background",],"recovery_modules.csv")
write_out(joint[joint$candidate,],"opposite_candidates.csv")
scenarios<-lapply(pb_models,function(z) make_sets(z,map_high,mouse_fits$Fibroblast$de))
scenarios$all_confidence_one2one<-make_sets(pb_models$unpaired,map,mouse_fits$Fibroblast$de)
msx<-read_counts("mouse_fibro_strict_counts.csv.gz")
msm<-read.csv(file.path(inputs,"mouse_fibro_strict_metadata.csv"))
msm$group<-relevel(factor(msm$group),ref="Tropicamide")
msd<-model.matrix(~group,msm)
msfit<-fit_de(msx[,msm$column,drop=FALSE],msd,"groupSaline")
msfit$de<-add_symbols(msfit$de,"mouse")
write_out(msfit$de,"mouse_fibro_strict_DE.csv")
scenarios$mouse_strict<-make_sets(pb_models$unpaired,map_high,msfit$de)
scenarios$estimated_gene_correlation<-joint
scenario_members<-do.call(rbind,lapply(names(scenarios),function(nm) {
  x<-scenarios[[nm]];x$scenario<-nm;x
}))
write_out(scenario_members,"all_scenario_membership.csv")
write_out(do.call(rbind,lapply(names(scenarios),function(nm) {
  x<-scenarios[[nm]];a<-joint$chick_gene_id[joint$candidate];b<-x$chick_gene_id[x$candidate]
  data.frame(scenario=nm,common=nrow(x),up=sum(x$module=="Recovery_up"),
    down=sum(x$module=="Recovery_down"),candidates=sum(x$candidate),
    candidate_overlap=length(intersect(a,b)),candidate_jaccard=length(intersect(a,b))/length(union(a,b)))
})),"module_sensitivity_summary.csv")

camera_sets<-function(fit,j,species="chick",gp_high=TRUE,correlation=.01) {
  ids<-j$chick_gene_id
  if(species=="gp") {
    mp<-gpmap
    if(gp_high) mp<-mp[!is.na(mp$is_high_confidence)&mp$is_high_confidence==1,]
    ids<-mp$gp_gene_id[match(j$chick_gene_id,mp$chick_gene_id)]
  }
  sets<-list(Recovery_down=ids[j$module=="Recovery_down"],
    Recovery_up=ids[j$module=="Recovery_up"],Opposite_candidates=ids[j$candidate])
  idx<-lapply(sets,function(g)which(rownames(fit$y)%in%g))
  eligible<-lengths(idx)>=5
  if(!any(eligible))return(NULL)
  v<-voom(fit$y,fit$design,plot=FALSE)
  z<-camera(v,index=idx[eligible],design=fit$design,contrast=fit$coef,inter.gene.cor=correlation)
  if(!"Correlation"%in%colnames(z))z$Correlation<-correlation
  z$module<-rownames(z);z$FDR_within_contrast<-z$FDR;z$FDR<-NULL
  z
}
allcam<-list();allmeta<-list();bulkde<-list();bulk_summary<-list()
for(gse in c("GSE227724","GSE261232","GSE256200")) {
  filename<-switch(gse,GSE227724="GSE227724_rawCounts_dnrs_142samp.txt.gz",
    GSE261232="GSE261232_raw_counts_DNRS_144samp.txt.gz",GSE256200="GSE256200_count_table.txt.gz")
  raw<-read.delim(gzfile(file.path(inputs,filename)),check.names=FALSE)
  xall<-as.matrix(raw[,-1,drop=FALSE]);storage.mode(xall)<-"double"
  rownames(xall)<-sub("\\..*$","",raw[[1]])
  if(anyDuplicated(rownames(xall)))xall<-rowsum(xall,rownames(xall),reorder=FALSE)
  meta<-read.csv(file.path(inputs,paste0(gse,"_metadata.csv")))
  stopifnot(!anyDuplicated(meta$count_column))
  if(gse!="GSE256200") {
    meta<-meta[tolower(meta$tissue)=="choroid",]
    timefield<-grep("time",names(meta),value=TRUE)[1]
    birdfield<-grep("chick|bird|animal",names(meta),value=TRUE)[1]
    eyefield<-grep("eye.*treatment|treatment",names(meta),value=TRUE)[1]
    contrasts<-sort(unique(meta[[timefield]]))
  }else contrasts<-c("Choroid","Sclera")
  for(tp in contrasts) {
    if(gse!="GSE256200") {
      m<-meta[meta[[timefield]]==tp,]
      pairs<-table(m[[birdfield]],m[[eyefield]])
      good<-rownames(pairs)[apply(pairs,1,function(v)length(v)==2&&all(v==1))]
      m<-m[m[[birdfield]]%in%good,]
      m$bird<-factor(m[[birdfield]])
      m$eye<-relevel(factor(m[[eyefield]]),ref="open")
      d<-model.matrix(~bird+eye,m);coef<-grep("^eye",colnames(d),value=TRUE)[1]
      species<-"chick"
    }else {
      m<-meta[meta$tissue==tp&meta$strain%in%c("Albino Hyperopic","Albino Myopic"),]
      m$strain<-relevel(factor(m$strain),ref="Albino Hyperopic")
      d<-model.matrix(~strain,m);coef<-"strainAlbino Myopic";species<-"gp"
    }
    stopifnot(all(m$count_column%in%colnames(xall)))
    fit<-fit_de(xall[,m$count_column,drop=FALSE],d,coef)
    de<-fit$de;de$dataset<-gse;de$contrast<-tp
    if(species=="gp")de$chick_gene_id<-gpmap$chick_gene_id[match(de$gene_id,gpmap$gp_gene_id)] else de$chick_gene_id<-de$gene_id
    de$symbol<-keys$chick$symbol[match(de$chick_gene_id,keys$chick$gene_id)]
    write_out(de,paste0(gse,"_",tp,"_DE.csv"))
    bulkde[[paste(gse,tp)]]<-de
    m$dataset<-gse;m$analysis_contrast<-tp
    allmeta[[paste(gse,tp)]]<-m[,c("gsm","count_column","dataset","analysis_contrast")]
    bulk_summary[[paste(gse,tp)]]<-data.frame(dataset=gse,contrast=tp,samples=nrow(m),
      genes=nrow(de),FDR05=sum(de$FDR<.05))
    for(scenario in names(scenarios)) {
      cam<-camera_sets(fit,scenarios[[scenario]],species,gp_high=scenario!="all_confidence_one2one",
        correlation=if(scenario=="estimated_gene_correlation")NA_real_ else .01)
      if(is.null(cam))next
      cam$dataset<-gse;cam$contrast<-tp;cam$scenario<-scenario
      allcam[[paste(gse,tp,scenario)]]<-cam
    }
  }
  message("Bulk complete: ",gse)
}
cam<-do.call(rbind,allcam);rownames(cam)<-NULL
cam$FDR_across_time<-ave(cam$PValue,interaction(cam$dataset,cam$scenario),FUN=function(x)p.adjust(x,"BH"))
cam$FDR_across_chick_stages<-NA_real_
for(s in unique(cam$scenario)) {
  ix<-cam$scenario==s&cam$dataset%in%c("GSE227724","GSE261232")
  cam$FDR_across_chick_stages[ix]<-p.adjust(cam$PValue[ix],"BH")
}
write_out(cam,"CAMERA_all_scenarios.csv")
write_out(cam[cam$scenario=="unpaired",],"CAMERA_primary.csv")
write_out(do.call(rbind,allmeta),"bulk_samples_used.csv")
write_out(do.call(rbind,bulk_summary),"bulk_summary.csv")
deall<-do.call(rbind,bulkde)
candde<-deall[deall$chick_gene_id%in%joint$chick_gene_id[joint$candidate],]
write_out(candde,"candidate_all_bulk_effects.csv")

# Gene-distribution tests are descriptive context, not independent validation.
celltests<-list();cellvalues<-list()
for(ct in names(mouse_fits)) {
  z<-merge(joint[,c("mouse_gene_id","chick_logFC","module","candidate","symbol")],
    mouse_fits[[ct]]$de,by.x="mouse_gene_id",by.y="gene_id",suffixes=c("_chick","_mouse"))
  z$expected_effect<- -sign(z$chick_logFC)*z$logFC
  z$celltype<-ct;cellvalues[[ct]]<-z
  for(mod in c("Recovery_down","Recovery_up","Opposite_candidates")) {
    ix<-if(mod=="Opposite_candidates")z$candidate else z$module==mod
    effect<-if(mod=="Opposite_candidates")z$expected_effect else if(mod=="Recovery_down")z$logFC else -z$logFC
    if(sum(ix)<5||sum(!ix)<5)next
    w<-wilcox.test(effect[ix],effect[!ix],alternative="greater",exact=FALSE,correct=FALSE)
    celltests[[paste(ct,mod)]]<-data.frame(celltype=ct,module=mod,n_genes=sum(ix),
      median_expected_effect=median(effect[ix]),PValue=w$p.value)
  }
}
ct<-do.call(rbind,celltests);ct$FDR<-p.adjust(ct$PValue,"BH")
write_out(ct,"celltype_directional_tests.csv")
write_out(do.call(rbind,cellvalues),"celltype_gene_effects.csv")

clx<-read_counts("mouse_fibro_clusters_counts.csv.gz")
clm<-read.csv(file.path(inputs,"mouse_fibro_clusters_metadata.csv"))
clusters<-list();composition<-list()
for(cl in unique(clm$celltype)) {
  m<-clm[clm$celltype==cl,]
  if(nrow(m)!=6||any(m$n_cells<10))next
  m$group<-relevel(factor(m$group),ref="Tropicamide")
  d<-model.matrix(~group,m)
  f<-fit_de(clx[,m$column,drop=FALSE],d,"groupSaline")
  write_out(add_symbols(f$de,"mouse"),paste0("mouse_cluster_",cl,"_DE.csv"))
  clusters[[cl]]<-data.frame(cluster=cl,min_cells=min(m$n_cells),genes=nrow(f$de),
    FDR05=sum(f$de$FDR<.05),analysis_role=if(min(m$n_cells)>=30)"primary" else "exploratory")
  totals<-tapply(clm$n_cells,clm$sample,sum)
  m$logit<-qlogis((m$n_cells+.5)/(totals[m$sample]+1))
  co<-summary(lm(logit~group,m))$coefficients["groupSaline",]
  composition[[cl]]<-data.frame(cluster=cl,logit_effect=co[1],SE=co[2],PValue=co[4])
}
write_out(do.call(rbind,clusters),"mouse_subcluster_summary.csv")
cc<-do.call(rbind,composition);cc$FDR<-p.adjust(cc$PValue,"BH")
write_out(cc,"mouse_subcluster_composition.csv")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
message("REANALYSIS COMPLETE")
