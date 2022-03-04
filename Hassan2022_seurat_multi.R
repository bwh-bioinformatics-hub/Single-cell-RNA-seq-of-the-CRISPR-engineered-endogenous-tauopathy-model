# title: "scRNA-seq data analysis of dissected fly brain, Hassan project, 2022"
# author: "Tingting"
# date: "02/05/2022"
# usage: Rscript 2_Hassan2022_seurat_multi.R /data/bioinformatics/projects/hassan2022/scr/ /data/bioinformatics/projects/hassan2022/output/02_cellranger_count_expectedCells/ /data/bioinformatics/projects/hassan2022/output/04_seurat_multi/hassan/ /data/bioinformatics/projects/hassan2022/output/03_scrublet/ CRN00224913,CRN00224914,CRN00224915,CRN00224916,CRN00224917,CRN00224918,CRN00224919,CRN00224920 hassan2022 P251L,control,TauKO,P251L,P251L,TauKO,control,TauKO VAChT,VGlut,Gad1,Vmat,moody,ey,prt,sNPF,trio,hth,bsh,Eaat1,Lim3,svp,Vsx1,eya,Lim1,ato,acj6,Crz,SerT,Tdc2,ple,alrm,wrapper,Hml,otp,Oaz,C15,tnc,kn,CG14687,Poxn,Dh31,grn,HLH3B,Dr,Sox21b,CNMa,cry,Clk,lncRNA:CR45566,Pdfr,gl,Pdf,CG17777,CG18599,tim,vri,Dh44,Nplp1,CG2016,CG33777,CG5910,otk2,lncRNA:CR43856,Sulf1,ort,Pka-C3,Tk,Hmx,AstA,CG10257,Vsx2,CG14989,Ets65A,hbn,CG13698,Dll,sosie,Mmp2,tup,Eip63F-1,msi,CG10804,bru3,dati,mbl,lncRNA:CR44024,CG9650,CG4577,ap,scro,CG14757,DIP-theta,beat-Ia,side-IV,CG42750,Drgx,CNMaR,Sox102F,SoxN,CG14340,cv-c

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(stringr)
library(SoupX)
library(reticulate)
library(glmpca)
library(SeuratWrappers)
library(scry)
library(reticulate)
library(monocle3)
library(FlexDotPlot)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)

message("Step1: settings")
pwd = args[1]
indir = args[2]
outdir = args[3]
scrubletdir = args[4]
samples = unlist(strsplit(args[5], ','))
projectName = args[6]
treatment = unlist(strsplit(args[7], ','))
markers = unlist(strsplit(args[8], ','))

message("Set working dir")
setwd(pwd)

message("Options for cell clustering algorithm, 1=louvain, 2=GLMPCA, 3= leiden, 0=louvain and GLMPCA and leiden")
flag=1

message("If using Leiden algorithm in FindMarkers")
use_condaenv("r_leiden", required=TRUE)
py_config()


message("Step 2: Pre-processing")
message("Remove ambient RNA by SoupX")
data.10x = list()
for (sample in samples) {
  data.10x[[sample]] = load10X(paste0(indir, sample, "/outs/"))
  data.10x[[sample]] <- autoEstCont(data.10x[[sample]], forceAccept=TRUE) #
  print(sample)
  data.10x[[sample]] <- adjustCounts(data.10x[[sample]])
}
message("Create Seurat object after SoupX")
scrna.list = list()
for (sample in samples) {
    scrna.list[[sample]] = CreateSeuratObject(counts = data.10x[[sample]], min.cells=3, project=sample)
}
message("Remove raw data to save memory")
rm(data.10x)
message("Add percent.mt and percent.rb to cell level metadata")
for (sample in samples) {
  scrna.list[[sample]][["percent.mt"]] <- PercentageFeatureSet(scrna.list[[sample]], pattern = "^mt:") 
  scrna.list[[sample]][["percent.rb"]] <- PercentageFeatureSet(scrna.list[[sample]], pattern = "^Rp[LS]")
}
message("Run doublet detection scripts")
system2(command = "bash", args = c("run_scrublet_multi.sh"))
message("Read in doublet scores")
for (sample in samples){
  doublet_scores <- scan(paste0(scrubletdir, sample, "_srublet.score"))
  predicted_doublets <- scan(paste0(scrubletdir, sample, "_srublet.logic"))   
  ds <- as.data.frame(cbind(doublet_scores, predicted_doublets))
  ds$predicted_doublets <- as.logical(ds$predicted_doublets)
  rownames(ds) <- rownames(scrna.list[[sample]]@meta.data) 
  scrna.list[[sample]] <- AddMetaData(scrna.list[[sample]], ds)
  scrna.list[[sample]] <- subset(scrna.list[[sample]], subset=predicted_doublets == FALSE)
}

message("Add sample name")
for (sample in samples){
  treatment_labels <- sample(x = treatment, size = ncol(x = scrna.list[[sample]]), replace = TRUE)
  scrna.list[[sample]]$treatment <- treatment_labels
}

message("Step 3: QC")
qc_cutoff = 3
mito_cutoff = 10
for (sample in samples){
  mean.nCount <- mean(scrna.list[[sample]]@meta.data$nCount_RNA)
  sd.nCount <- sd(scrna.list[[sample]]@meta.data$nCount_RNA)
  mean.nFeature <- mean(scrna.list[[sample]]@meta.data$nFeature_RNA)
  sd.nFeature <- sd(scrna.list[[sample]]@meta.data$nFeature_RNA)
  scrna.list[[sample]] <- subset(scrna.list[[sample]], subset = nCount_RNA > mean.nCount - qc_cutoff*sd.nCount & nCount_RNA < mean.nCount + qc_cutoff*sd.nCount & nFeature_RNA > mean.nFeature - qc_cutoff*sd.nFeature & nFeature_RNA < mean.nFeature + qc_cutoff*sd.nFeature & percent.mt < mito_cutoff)
}


message("Step 4: Integration")
scrna.list <- lapply(X = scrna.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = scrna.list)
scrna.anchors <- FindIntegrationAnchors(object.list = scrna.list, anchor.features = features)
scrna.combined <- IntegrateData(anchorset = scrna.anchors)
DefaultAssay(scrna.combined) <- "integrated"
scrna.combined <- ScaleData(scrna.combined, verbose = FALSE)
scrna.combined <- RunPCA(scrna.combined, npcs = 15, verbose = FALSE)
scrna.combined <- RunUMAP(scrna.combined, reduction = "pca", dims = 1:15)
scrna.combined <- FindNeighbors(scrna.combined, reduction = "pca", dims = 1:15)
scrna.combined <- FindClusters(scrna.combined, resolution = 0.5)

p1 <- DimPlot(scrna.combined, reduction = "umap", group.by = "orig.ident")
p1
pdf(file=paste0(outdir, "combined.umap.colorBySample.pdf"))
print(p1 + coord_fixed())
dev.off()

p2 <- DimPlot(scrna.combined, reduction = "umap", label = TRUE, repel = TRUE)
p2
pdf(paste0(outdir, "combined.umap.colorByCluster.pdf"))
print(p2 + coord_fixed())
dev.off()

p3 <- DimPlot(scrna.combined, reduction = "umap", split.by = "orig.ident", ncol = 4)
p3
pdf(paste0(outdir, "combined.umap.samples.pdf"), width = 12, height = 8)
print(p3 + coord_fixed())
dev.off()

p4 <- DimPlot(scrna.combined, reduction = "umap", group.by = "treatment")
p4
pdf(file=paste0(outdir, "combined.umap.colorByTreatment.pdf"))
print(p4 + coord_fixed())
dev.off()

message("Step 5: Identify conserved cell type markers")
DefaultAssay(scrna.combined) <- "RNA"
nk.markers <- FindConservedMarkers(scrna.combined, ident.1 = 0, grouping.var = "orig.ident", verbose = FALSE)
head(nk.markers)

message("Step 6: Feature plot")
plots = FeaturePlot(scrna.combined, features = unlist(markers, use.names = F), ncol =4, by.col =F, combine = F)
pdf(paste0(outdir, "combined.markergenes.seperate.pdf"))
sapply(plots, print)
dev.off()

message("Step 7: Dotplot before annotation")
p1 <- DotPlot(scrna.combined, features = markers, cols = c("blue", "red", "green"), split.by = "treatment") + RotatedAxis()
p1
pdf(file=paste0(outdir, "combined.dotplot.beforeAnnotationA.pdf"), height = 20, width = 25)
print(p1 + coord_fixed())
dev.off()

p1 <- DotPlot(scrna.combined, features = markers) + RotatedAxis()
p1
pdf(file=paste0(outdir, "combined.dotplot.beforeAnnotationAa.pdf"), height = 20, width = 25)
print(p1 + coord_fixed())
dev.off()

message("Step 8: Dotplot before annotation with dendrogram")
# dengrogram function
new_dotplot <- function(object = NULL, features = NULL, group.by = NULL, genes.on.x = TRUE,
                        size.breaks.values = NULL, color.breaks.values = c(-3, -2, -1, 0, 1, 2, 3), shape.scale = 12,
                        dend_x_var = "Average expression", dend_y_var = "Average expression",
                        cols.use = c("lightgrey", "blue"), scale.by = "radius", col.min = -2.5, col.max = 2.5,
                        dot.min = 0) {
  scale.func <- switch(EXPR = scale.by, size = scale_size,
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  
  data.features <- FetchData(object = object, vars = features)
  object[[group.by, drop = TRUE]]
  data.features$id <- object[[group.by, drop = TRUE]]
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident,
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat:::PercentAbove,
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot ==
                                                     x, "avg.exp"]
                             data.use <- scale(x = data.use)
                             data.use <- MinMax(data = data.use, min = col.min,
                                                max = col.max)
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot,
                                    levels = rev(x = features))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  
  if(genes.on.x){
    data.final <- data.frame(features.plot = data.plot$features.plot, id = data.plot$id)
    data.final <- cbind(data.final, data.plot[,c(2,5)])
  }
  else {
    data.final <- data.frame(id = data.plot$id, features.plot = data.plot$features.plot)
    data.final <- cbind(data.final, data.plot[,c(2,5)])
  }
  colnames(data.final)[3:4] <- c("Percent expressed", "Average expression")
  dot_plot(data.final, size_var = "Percent expressed", "Average expression",
           dend_y_var = dend_y_var, dend_x_var = dend_x_var,
           dist = "euclidean", hclust_method = "ward.D2", x.lab.pos = "bottom",
           display_max_sizes = FALSE, size.breaks.values = size.breaks.values,
           shape.scale = shape.scale, color.breaks.values = color.breaks.values, cols.use = cols.use, y.lab.size.factor=0.05)
}  # modify y lab text size here


pdf(file=paste0(outdir, "combined.dotplot.beforeAnnotationB.pdf"), height = 20, width = 30)
dp = DotPlot(scrna.combined, features = markers) + RotatedAxis()
dotplot = dot_plot(dp$data[,c(3,4,1,2,5)], shape_var = "pct.exp", col_var = "avg.exp.scaled", shape_legend = "Percent Expressed", col_legend = "Average Expression", x.lab.pos = "bottom", dend_x_var = c("pct.exp","avg.exp.scaled"), dend_y_var = c("pct.exp","avg.exp.scaled"), hclust_method = "ward.D2", do.return=T, y.lab.size.factor=0.1)
plot1 <- new_dotplot(scrna.combined,
                       features = dotplot$plot$grobs[[10]]$label,
                       group.by = "seurat_clusters",
                       shape.scale = 6,
                       color.breaks.values = c(-2, -1, 0, 1, 2),
                       size.breaks.values = c(0, 25, 50, 75, 100)) +
    theme(axis.text = element_text(size = 12))
print(plot1)
dev.off()

message("Step 9: Cell type annotation. CHANGE HERE!!!")
scrna.combined <- RenameIdents(scrna.combined, `0` = "Dopaminergic", `1` = "Dopaminergic", `2` = "Unknown", `3` = "Dopaminergic", `4` = "Glutamatergic", `5` = "C2/C3", `6` = "Kenyon", `7` = "T4/T5", `8` = "L1/L2/L3/L4/L5", `9` = "Glutamatergic", `10` = "Unknown", `11` = "Mi1", `12` = "C2/C3", `13` = "Dopaminergic", `14` = "Lawf2", `15` = "Unknown", `16` = "Lawf2", `17` = "Lawf2", `18` = "Glial", `19` = "Unknown", `20` = "Perineurial", `21` = "Glutamatergic", `22` = "Glutamatergic", `23` = "TMY14", `24` = "Dm8/Tm5c", `25` = "Cholinergic", `26` = "Poxn")

message("Step 10: UMAP, cell annotated")
p1 <- DimPlot(scrna.combined, label = TRUE, repel = TRUE, pt.size = 0.5)
p1
pdf(file=paste0(outdir, "combined.umap.annotated.pdf"))
print(p1 + coord_fixed())
dev.off()

