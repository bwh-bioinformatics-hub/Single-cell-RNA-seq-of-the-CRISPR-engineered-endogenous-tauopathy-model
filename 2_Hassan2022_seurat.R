# title: "scRNA-seq data analysis of dissected fly brain, Hassan project, 2022"
# author: "Tingting"
# date: "02/05/2022"
# usage: Rscript 0_Hassan2022_seurat.R /Volumes/BIOINFORMATICS/projects/hassan2022/scr/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/02_cellranger_count/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/03_seurat/ CRN00224913,CRN00224914,CRN00224915,CRN00224916,CRN00224917,CRN00224918,CRN00224919,CRN00224920 Hassan2022
  
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
projectName = args[5]

# Set working dir
setwd(pwd)

# Step 1: Preparation

# If using Leiden algorithm in FindMarkers
use_condaenv("r_leiden", required=TRUE)
py_config()

# Step 2: Ambient RNA removal
data.10x = list()
i = 1
  for (sample in samples) {
    data.10x[[i]] = load10X(paste0(indir, sample, "/outs/"))
    data.10x[[i]] <- autoEstCont(data.10x[[i]])
    print(i)
    print(sample)
    data.10x[[i]] <- adjustCounts(data.10x[[i]])
    i=i+1
  }

# Step 3: Create Seurat object after SoupX
scrna.list = list()
for (j in 1:length(data.10x)) {
    scrna.list[[j]] = CreateSeuratObject(counts = data.10x[[j]], min.cells=3, project=samples[j])
    scrna.list[[j]][["DataSet"]] = samples[j]
}

# Remove raw data to save memory
rm(data.10x)

# Step 4: QC, percent.mt and percent.rb
for (i in 1:length(samples)) {
  scrna.list[[i]][["percent.mt"]] <- PercentageFeatureSet(scrna.list[[i]], pattern = "^MT-") 
  scrna.list[[i]][["percent.rb"]] <- PercentageFeatureSet(scrna.list[[i]], pattern = "^Rp[LS]")
}

# Step 5: Doublet detection
ssystem2(command = "bash",
        args = c("run_scrublet_multi.sh"))


