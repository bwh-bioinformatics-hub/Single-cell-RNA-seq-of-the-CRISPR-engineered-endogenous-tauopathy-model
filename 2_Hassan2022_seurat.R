# title: "scRNA-seq data analysis of dissected fly brain, Hassan project, 2022"
# author: "Tingting"
# date: "02/05/2022"
# usage: Rscript 0_Hassan2022_seurat.R /Volumes/BIOINFORMATICS/projects/hassan2022/scr/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/02_cellranger_count_expectedCells/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/03_seurat/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/03_scrublet/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/04_seurat/ CRN00224913,CRN00224914,CRN00224915,CRN00224916,CRN00224917,CRN00224918,CRN00224919,CRN00224920 hassan2022
  
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

args <- commandArgs(trailingOnly = TRUE)

# settings
pwd = args[1]
indir = args[2]
outdir = args[3]
scrubletdir = args[4]
seuratdir = args[5]
samples = args[6]
projectName = args[7]

# Set working dir
setwd(pwd)

# Options for cell clustering algorithm
flag=1 # 1=louvain, 2=GLMPCA, 3= leiden, 0=louvain and GLMPCA and leiden 

# If using Leiden algorithm in FindMarkers
use_condaenv("r_leiden", required=TRUE)
py_config()


# Step 2: Pre-processing
# Remove ambient RNA by SoupX
data.10x = list()
for (sample in samples) {
  data.10x[[sample]] = load10X(paste0(indir, sample, "/outs/"))
  data.10x[[sample]] <- autoEstCont(data.10x[[sample]], forceAccept=TRUE) #
  print(sample)
  data.10x[[sample]] <- adjustCounts(data.10x[[sample]])
}
# Create Seurat object after SoupX
scrna.list = list()
for (sample in samples) {
    scrna.list[[sample]] = CreateSeuratObject(counts = data.10x[[sample]], min.cells=3, project=sample)
}
# Remove raw data to save memory
rm(data.10x)
# Add percent.mt and percent.rb to cell level metadata
for (sample in samples) {
  scrna.list[[sample]][["percent.mt"]] <- PercentageFeatureSet(scrna.list[[sample]], pattern = "^mt:") 
  scrna.list[[sample]][["percent.rb"]] <- PercentageFeatureSet(scrna.list[[sample]], pattern = "^Rp[LS]")
}
# Run doublet detection scripts
system2(command = "bash", args = c("run_scrublet_multi.sh"))
# Read in doublet scores
for (sample in samples){
  doublet_scores <- scan(paste0(scrubletdir, sample, "_srublet.score"))
  predicted_doublets <- scan(paste0(scrubletdir, sample, "_srublet.logic"))   
  ds <- as.data.frame(cbind(doublet_scores, predicted_doublets))
  ds$predicted_doublets <- as.logical(ds$predicted_doublets)
  rownames(ds) <- rownames(scrna.list[[sample]]@meta.data) 
  scrna.list[[sample]] <- AddMetaData(scrna.list[[sample]], ds)
  scrna.list[[sample]] <- subset(scrna.list[[sample]], subset=predicted_doublets == FALSE)
}


# Step 3: QC
# Feature plot before QC
pdf(file = paste0(seuratdir, "QC.plot.before.pdf"), width = 8, height = 8)
for (sample in samples){
  plot1 <- VlnPlot(scrna.list[[sample]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10)) + theme(aspect.ratio=40/10)
  print(plot1 + coord_fixed())
  plot2 <- FeatureScatter(scrna.list[[sample]], feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(aspect.ratio=10/10)
  print(plot2 + coord_fixed())
  plot3 <- FeatureScatter(scrna.list[[sample]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(aspect.ratio=10/10)
  print(plot3 + coord_fixed())
}
dev.off()

# Filtered cells with 3SD of mean nCount and nFeature, percent of mito
qc_cutoff = 3
mito_cutoff = 10
for (sample in samples){
  mean.nCount <- mean(scrna.list[[sample]]@meta.data$nCount_RNA)
  sd.nCount <- sd(scrna.list[[sample]]@meta.data$nCount_RNA)
  mean.nFeature <- mean(scrna.list[[sample]]@meta.data$nFeature_RNA)
  sd.nFeature <- sd(scrna.list[[sample]]@meta.data$nFeature_RNA)
  scrna.list[[sample]] <- subset(scrna.list[[sample]], subset = nCount_RNA > mean.nCount - qc_cutoff*sd.nCount & nCount_RNA < mean.nCount + qc_cutoff*sd.nCount & nFeature_RNA > mean.nFeature - qc_cutoff*sd.nFeature & nFeature_RNA < mean.nFeature + qc_cutoff*sd.nFeature & percent.mt < mito_cutoff)
}

# Feature plot after QC
pdf(file = paste0(seuratdir, "QC.plot.after.pdf"))
for (sample in samples){
  plot1 <- VlnPlot(scrna.list[[sample]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10)) + theme(aspect.ratio=40/10)
  print(plot1 + coord_fixed())
  plot2 <- FeatureScatter(scrna.list[[sample]], feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(aspect.ratio=10/10)
  print(plot2 + coord_fixed())
  plot3 <- FeatureScatter(scrna.list[[sample]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(aspect.ratio=10/10)
  print(plot3 + coord_fixed())
}
dev.off()


# Step 4: Sample processing (Normalization, Find variable features, Data scaling)
for (sample in samples){
  scrna.list[[sample]] <- NormalizeData(scrna.list[[sample]], normalization.method = "LogNormalize", scale.factor = 10000)
  scrna.list[[sample]] <- FindVariableFeatures(scrna.list[[sample]], selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(scrna.list[[sample]])
  scrna.list[[sample]] <- ScaleData(scrna.list[[sample]], features = all.genes, verbose = FALSE)
}


# Step 5: Determine the ‘dimensionality’ of the dataset
pdf(file = paste0(seuratdir, "elbow.plot.pdf"), width = 8, height = 8)
for (sample in samples){
  scrna.list[[sample]] <- RunPCA(scrna.list[[sample]], verbos = FALSE)
  plot1 <- ElbowPlot(scrna.list[[sample]]) + ggtitle(sample) + theme(aspect.ratio=5/10) + theme(plot.margin = unit(c(3, 3, 3, 3), "cm"))
  print(plot1 + coord_fixed())
}
dev.off()


# Step 6: Cell clustering
# Save a copy of seurat object for GLMPCA and Leiden exploration
scrna.default.list = list()
for (sample in samples){
  scrna.default.list[[sample]] <- scrna.list[[sample]]
}
# Cell clustering using default settings: PCA, Louvain, CHANGE dims according to elbow plot !!!
# Run PCA with elbow plot determined PCs
# Find neighbors
# Find clusters
# Run UMAP
if flag==1 | flag==0{
  pdf(file = paste0(seuratdir, "cluster.umap.louvain.pdf"), width = 6, height = 6)
for (sample in samples){
  scrna.list[[sample]] <- RunPCA(scrna.list[[sample]], npcs = 10, verbose = FALSE) 
  scrna.list[[sample]] <- FindNeighbors(scrna.list[[sample]], reduction = "pca", dims = 1:15) # picked 15
  scrna.list[[sample]] <- FindClusters(scrna.list[[sample]], resolution = 0.5)
  scrna.list[[sample]] <- RunUMAP(scrna.list[[sample]], reduction = "pca", dims = 1:10) 
  plot1 <- DimPlot(scrna.list[[sample]], reduction = "umap", label = TRUE, pt.size = 1) + ggtitle(sample)
  print(plot1 + coord_fixed())
}
dev.off()
}

# Cell clustering using GLMPCA
if flag==2 | flag==0{
  pdf(file = paste0(seuratdir, "cluster.umap.glmpca.pdf"), width = 6, height = 6)
scrna.glmpca.list = list()
for (sample in samples){
  scrna.glmpca.list[[sample]] <- RunGLMPCA(scrna.default.list[[sample]], L = 10)
  scrna.glmpca.list[[sample]] <- FindNeighbors(scrna.glmpca.list[[sample]], reduction = "glmpca", dims = 1:10)
  scrna.glmpca.list[[sample]] <- FindClusters(scrna.glmpca.list[[sample]], resolution = 0.5)
  scrna.glmpca.list[[sample]] <- RunUMAP(scrna.glmpca.list[[sample]], reduction = "glmpca", dims = 1:10)
  plot1 <- DimPlot(scrna.glmpca.list[[sample]], reduction = "umap", label = TRUE, pt.size = 1) + ggtitle(sample)
  print(plot1 + coord_fixed())
}
dev.off()
}

# Cell clustering using Leiden
if flag==3 | flag==0{
pdf(file = paste0(seuratdir, "cluster.umap.leiden.pdf"), width = 6, height = 6)
scrna.leiden.list = list()
for (sample in samples){
  scrna.leiden.list[[sample]] <- RunPCA(scrna.default.list[[sample]], npcs = 10, verbose = FALSE) 
  scrna.leiden.list[[sample]] <- FindNeighbors(scrna.leiden.list[[sample]], reduction = "pca", dims = 1:10) # picked 10
  scrna.leiden.list[[sample]] <- FindClusters(scrna.leiden.list[[sample]], resolution = 0.5, algorithm = 4)
  scrna.leiden.list[[sample]] <- RunUMAP(scrna.leiden.list[[sample]], reduction = "pca", dims = 1:10)
  plot1 <- DimPlot(scrna.leiden.list[[sample]], reduction = "umap", label = TRUE, pt.size = 1) + ggtitle(sample)
  print(plot1 + coord_fixed())
}
dev.off()  
}


# Step 7: Save the data
for (sample in samples){
  saveRDS(scrna.list[[sample]], paste0(seuratdir, sample, ".seurat.rds"))
  scrna.list[[sample]] <- readRDS(paste0(seuratdir, sample, ".seurat.rds"))
}
```