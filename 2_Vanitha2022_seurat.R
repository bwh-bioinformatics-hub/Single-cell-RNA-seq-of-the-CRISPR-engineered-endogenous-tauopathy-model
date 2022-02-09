# title: "scRNA-seq data analysis of dissected fly brain, Hassan project, 2022"
# author: "Tingting"
# date: "02/05/2022"
# usage: Rscript 0_Hassan2022_seurat.R /Volumes/BIOINFORMATICS/projects/hassan2022/scr/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/02_cellranger_count/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/03_seurat/ CRN00224921,CRN00224922,CRN00224923,CRN00224924,CRN00224925,CRN00224926
  
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
samples = args[4]

# Set working dir
setwd(pwd)

# Step1: Preparation
#srublet.score="/Volumes/bioinformatics/projects/hassan2021/output/05_Scrublet/srublet.score"
#srublet.logic="/Volumes/bioinformatics/projects/hassan2021/output/05_Scrublet/srublet.logic"

# If using Leiden algorithm in FindMarkers
use_condaenv("r_leiden", required=TRUE)
py_config()

# Step2: Read in count table
for (i in 1:length(samples)) {
  data.10x[[i]] = Read10X(data.dir = paste0(indir, samples[[i]], "/outs/filtered_feature_bc_matrix"));
}

# Step3: Create Seurat object
for (j in 1:length(data.10x)) {
  scrna.list[[j]] = CreateSeuratObject(counts = data.10x[[j]], min.cells=3, project=samples[j]);
  scrna.list[[j]][["DataSet"]] = samples[j];
}

# Remove raw data to save memory
rm(data.10x);

# Step4: Merge Seurat objects into an object
sc.hassan <- scrna.list[[1]]
for (k in 2:length(samples.hassan)) {
  sc.hassan <- merge(sc.hassan, y=c(scrna.list[[k]]), add.cell.ids = samples.hassan[1:k], project="hassan2022");
}
#sc.hassan <- merge(scrna.list[[1]], y=c(scrna.list[[2]], scrna.list[[3]], scrna.list[[4]], scrna.list[[5]], scrna.list[[6]], scrna.list[[7]], scrna.list[[8]]), add.cell.ids = samples.hassan, project="hassan2022");

sc.vanitha <- scrna.list[[1]]
for (l in 2:length(samples.vanitha)) {
  sc.vanitha <- merge(sc.vanitha, y=c(scrna.list[[l]]), add.cell.ids = samples.hassan[1:l], project="hassan2022");
}
#sc.vanitha <- merge(scrna.list[[9]], y=c(scrna.list[[10]], scrna.list[[11]], scrna.list[[12]], scrna.list[[13]], scrna.list[[14]]), add.cell.ids = samples.vanitha, project="vanitha2022");
rm(scrna.list); 
str(sc.hassan@meta.data);
str(sc.vanitha@meta.data);
saveRDS(sc.hassan, file = sprintf("%s/MergedSeuratObjectHassan.rds", outdir));
saveRDS(sc.vanitha, file = sprintf("%s/MergedSeuratObjectVanitha.rds", outdir))






