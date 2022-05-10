# title: "scRNA-seq data analysis with Seurat"
# author: "Tingting Zhao"
# email: tzhao7@bwh.harvard.edu
# date: "02/05/2022"
  
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
library(googlesheets4)

args <- commandArgs(trailingOnly = TRUE)

pwd = args[1]
indir = args[2]
outdir = args[3]
scrubletdir = args[4]
samples = unlist(strsplit(args[5], ','))
projectName = args[6]
marker_link = args[7]
marker_sheet = args[8]
flag = args[9]
mtPattern = args[10]
rbPattern = args[11]
mitoCutoff = args[12]


message("Read in marker genes")
gsurl=marker_link
gs4_deauth()
markers = read_sheet(gsurl, sheet = marker_sheet) %>%
  select(Markers) %>% as.list()

message("Set working dir")
setwd(pwd)

message("If using Leiden algorithm in FindMarkers")
#use_condaenv("r_leiden", required=TRUE)
#py_config()


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
  scrna.list[[sample]][["percent.mt"]] <- PercentageFeatureSet(scrna.list[[sample]], pattern = mtPattern)
  scrna.list[[sample]][["percent.rb"]] <- PercentageFeatureSet(scrna.list[[sample]], pattern = rbPattern)
}

message("Rename nCount_RNA and nFeature_RNA")
for (sample in samples) {
  scrna.list[[sample]]$nUMI <- scrna.list[[sample]]$nCount_RNA
  scrna.list[[sample]]$nCount_RNA <- NULL
  scrna.list[[sample]]$nGene <- scrna.list[[sample]]$nFeature_RNA
  scrna.list[[sample]]$nFeature_RNA <- NULL
}

message("Run doublet detection scripts")
#system2(command = "bash", args = c("run_scrublet_multi.sh"))

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


message("Feature plot before QC")
pdf(file = paste0(outdir, "seurat_individual/", "QC.plot.before.pdf"), width = 8, height = 8)
for (sample in samples){
  plot1 <- VlnPlot(scrna.list[[sample]], features = c("nGene", "nUMI", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10)) + theme(aspect.ratio=40/10)
  print(plot1 + coord_fixed())
  plot2 <- FeatureScatter(scrna.list[[sample]], feature1 = "nUMI", feature2 = "percent.mt") + theme(aspect.ratio=10/10)
  print(plot2 + coord_fixed())
  plot3 <- FeatureScatter(scrna.list[[sample]], feature1 = "nUMI", feature2 = "nGene") + theme(aspect.ratio=10/10)
  print(plot3 + coord_fixed())
}
dev.off()

message("Filtered cells with 3SD of mean nCount and nFeature, percent of mito")
qc_cutoff = 3
mito_cutoff = mitoCutoff
for (sample in samples){
  mean.nCount <- mean(scrna.list[[sample]]@meta.data$nUMI)
  sd.nCount <- sd(scrna.list[[sample]]@meta.data$nUMI)
  mean.nFeature <- mean(scrna.list[[sample]]@meta.data$nGene)
  sd.nFeature <- sd(scrna.list[[sample]]@meta.data$nGene)
  scrna.list[[sample]] <- subset(scrna.list[[sample]], subset = nUMI > mean.nCount - qc_cutoff*sd.nCount & nUMI < mean.nCount + qc_cutoff*sd.nCount & nGene > mean.nFeature - qc_cutoff*sd.nFeature & nGene < mean.nFeature + qc_cutoff*sd.nFeature & percent.mt < mito_cutoff)
}

message("Feature plot after QC")
pdf(file = paste0(outdir, "seurat_individual/", "QC.plot.after.pdf"))
for (sample in samples){
  plot1 <- VlnPlot(scrna.list[[sample]], features = c("nGene", "nUMI", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10)) + theme(aspect.ratio=40/10)
  print(plot1 + coord_fixed())
  plot2 <- FeatureScatter(scrna.list[[sample]], feature1 = "nUMI", feature2 = "percent.mt") + theme(aspect.ratio=10/10)
  print(plot2 + coord_fixed())
  plot3 <- FeatureScatter(scrna.list[[sample]], feature1 = "nUMI", feature2 = "nGene") + theme(aspect.ratio=10/10)
  print(plot3 + coord_fixed())
}
dev.off()


message("Step 4: Sample processing, Normalization, Find variable features, Data scaling")
for (sample in samples){
  scrna.list[[sample]] <- NormalizeData(scrna.list[[sample]], normalization.method = "LogNormalize", scale.factor = 10000)
  scrna.list[[sample]] <- FindVariableFeatures(scrna.list[[sample]], selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(scrna.list[[sample]])
  scrna.list[[sample]] <- ScaleData(scrna.list[[sample]], features = all.genes, verbose = FALSE)
}


message("Step 5: Determine the ‘dimensionality’ of the dataset")
pdf(file = paste0(outdir, "seurat_individual/", "elbow.plot.pdf"), width = 8, height = 8)
for (sample in samples){
  scrna.list[[sample]] <- RunPCA(scrna.list[[sample]], verbos = FALSE)
  plot1 <- ElbowPlot(scrna.list[[sample]]) + ggtitle(sample) + theme(aspect.ratio=5/10) + theme(plot.margin = unit(c(3, 3, 3, 3), "cm"))
  print(plot1 + coord_fixed())
}
dev.off()


message("Step 6: Cell clustering")
scrna.default.list = list()
for (sample in samples){
  scrna.default.list[[sample]] <- scrna.list[[sample]]
}

message("Cell clustering using default settings: PCA, Louvain. CHANGE dims according to elbow plot !!!")
if (flag==1 | flag==0){
  pdf(file = paste0(outdir, "cluster.umap.louvain.pdf"), width = 6, height = 6)
for (sample in samples){
  scrna.list[[sample]] <- RunPCA(scrna.list[[sample]], npcs = 15, verbose = FALSE) 
  scrna.list[[sample]] <- FindNeighbors(scrna.list[[sample]], reduction = "pca", dims = 1:15) # picked 15
  scrna.list[[sample]] <- FindClusters(scrna.list[[sample]], resolution = 0.5)
  scrna.list[[sample]] <- RunUMAP(scrna.list[[sample]], reduction = "pca", dims = 1:15) 
  plot1 <- DimPlot(scrna.list[[sample]], reduction = "umap", label = TRUE, pt.size = 1) + ggtitle(sample)
  print(plot1 + coord_fixed())
}
dev.off()
}

message("Cell clustering using GLMPCA")
if (flag==2 | flag==0){
  pdf(file = paste0(outdir, "seurat_individual/", "cluster.umap.glmpca.pdf"), width = 6, height = 6)
scrna.glmpca.list = list()
for (sample in samples){
  scrna.glmpca.list[[sample]] <- RunGLMPCA(scrna.default.list[[sample]], L = 15)
  scrna.glmpca.list[[sample]] <- FindNeighbors(scrna.glmpca.list[[sample]], reduction = "glmpca", dims = 1:15)
  scrna.glmpca.list[[sample]] <- FindClusters(scrna.glmpca.list[[sample]], resolution = 0.5)
  scrna.glmpca.list[[sample]] <- RunUMAP(scrna.glmpca.list[[sample]], reduction = "glmpca", dims = 1:15)
  plot1 <- DimPlot(scrna.glmpca.list[[sample]], reduction = "umap", label = TRUE, pt.size = 1) + ggtitle(sample)
  print(plot1 + coord_fixed())
}
dev.off()
}

message("Cell clustering using Leiden")
if (flag==3 | flag==0){
pdf(file = paste0(outdir, "seurat_individual/", "cluster.umap.leiden.pdf"), width = 6, height = 6)
scrna.leiden.list = list()
for (sample in samples){
  scrna.leiden.list[[sample]] <- RunPCA(scrna.default.list[[sample]], npcs = 15, verbose = FALSE) 
  scrna.leiden.list[[sample]] <- FindNeighbors(scrna.leiden.list[[sample]], reduction = "pca", dims = 1:15) # picked 10
  scrna.leiden.list[[sample]] <- FindClusters(scrna.leiden.list[[sample]], resolution = 0.5, algorithm = 4)
  scrna.leiden.list[[sample]] <- RunUMAP(scrna.leiden.list[[sample]], reduction = "pca", dims = 1:15)
  plot1 <- DimPlot(scrna.leiden.list[[sample]], reduction = "umap", label = TRUE, pt.size = 1) + ggtitle(sample)
  print(plot1 + coord_fixed())
}
dev.off()  
}
