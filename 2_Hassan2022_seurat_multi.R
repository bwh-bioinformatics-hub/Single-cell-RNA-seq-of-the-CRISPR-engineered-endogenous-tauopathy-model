# title: "scRNA-seq data analysis of dissected fly brain, Hassan project, 2022"
# author: "Tingting"
# date: "02/05/2022"
# usage: Rscript 2_Hassan2022_seurat_multi.R /Volumes/BIOINFORMATICS/projects/hassan2022/scr/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/02_cellranger_count_expectedCells/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/03_seurat/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/03_scrublet/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/04_seurat_multi/ CRN00224913,CRN00224914,CRN00224915,CRN00224916,CRN00224917,CRN00224918,CRN00224919,CRN00224920,CRN00224921,CRN00224922,CRN00224923,CRN00224924,CRN00224925,CRN00224926 hassan2022
  
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


# Step4: Merge Seurat objects into an object
scrna <- scrna.list[[1]]
y.list=c()
for(k in 2:length(samples)){
  y.list=c(y.list, scrna.list[[k]])
}
scrna <- merge(scrna, y=y.list, add.cell.ids = samples, project=projectName)
#rm(scrna.list)
rm(y.list)
str(scrna@meta.data)


# Step5: Integration
# normalize and identify variable features for each dataset independently
scrna.list <- lapply(X = scrna.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = scrna.list)
# Perform integration
scrna.anchors <- FindIntegrationAnchors(object.list = scrna.list, anchor.features = features)
scrna.combined <- IntegrateData(anchorset = scrna.anchors)
# Perform an integrated analysis
DefaultAssay(scrna.combined) <- "integrated"
scrna.combined <- ScaleData(scrna.combined, verbose = FALSE)
scrna.combined <- RunPCA(scrna.combined, npcs = 15, verbose = FALSE)
scrna.combined <- RunUMAP(scrna.combined, reduction = "pca", dims = 1:15)
scrna.combined <- FindNeighbors(scrna.combined, reduction = "pca", dims = 1:15)
scrna.combined <- FindClusters(scrna.combined, resolution = 0.5)

p1 <- DimPlot(scrna.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(scrna.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(scrna.combined, reduction = "umap", split.by = "orig.ident")

# Identify conserved cell type markers
DefaultAssay(scrna.combined) <- "RNA"
nk.markers <- FindConservedMarkers(scrna.combined, ident.1 = 6, grouping.var = "orig.ident", verbose = FALSE)
head(nk.markers)

