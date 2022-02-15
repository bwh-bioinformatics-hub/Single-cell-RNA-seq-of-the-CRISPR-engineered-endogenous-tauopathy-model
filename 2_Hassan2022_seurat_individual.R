# title: "scRNA-seq data analysis of dissected fly brain, Hassan project, 2022"
# author: "Tingting"
# date: "02/05/2022"
# usage: Rscript 2_Hassan2022_seurat_individual.R /Volumes/BIOINFORMATICS/projects/hassan2022/scr/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/02_cellranger_count_expectedCells/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/04_seurat_individual/ /Volumes/BIOINFORMATICS/projects/hassan2022/output/03_scrublet/ CRN00224913,CRN00224914,CRN00224915,CRN00224916,CRN00224917,CRN00224918,CRN00224919,CRN00224920,CRN00224921,CRN00224922,CRN00224923,CRN00224924,CRN00224925,CRN00224926 hassan2022 VAChT,VGlut,Gad1,Vmat,moody,ey,prt,sNPF,trio,hth,bsh,Eaat1,Lim3,svp,Vsx1,eya,Lim1,ato,acj6,Crz,SerT,Tdc2,ple,alrm,wrapper,Hml,otp,Oaz,C15,tnc,kn,CG14687,Poxn,Dh31,grn,HLH3B,Dr,Sox21b,CNMa,cry,Clk,lncRNA:CR45566,Pdfr,gl,Pdf,CG17777,CG18599,tim,vri,Dh44,Nplp1,CG2016,CG33777,CG5910,otk2,lncRNA:CR43856,Sulf1,ort,Pka-C3,Tk,Hmx,AstA,CG10257,Vsx2,CG14989,Ets65A,hbn,CG13698,Dll,sosie,Mmp2,tup,Eip63F-1,msi,CG10804,CG343,bru3,dati,mbl,lncRNA:CR44024,CG9650,CG4577,ap,scro,CG14757,DIP-theta,beat-Ia,side-IV,CG42750,Drgx,CNMaR,Sox102F,SoxN,CG14340,cv-c
  
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

args <- commandArgs(trailingOnly = TRUE)

# settings
pwd = args[1]
indir = args[2]
outdir = args[3]
scrubletdir = args[4]
samples = args[5]
projectName = args[6]
markers = args[7]

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
pdf(file = paste0(outdir, "QC.plot.before.pdf"), width = 8, height = 8)
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
pdf(file = paste0(outdir, "QC.plot.after.pdf"))
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
pdf(file = paste0(outdir, "elbow.plot.pdf"), width = 8, height = 8)
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
if (flag==1 | flag==0){
  pdf(file = paste0(outdir, "cluster.umap.louvain.pdf"), width = 6, height = 6)
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
if (flag==2 | flag==0){
  pdf(file = paste0(outdir, "cluster.umap.glmpca.pdf"), width = 6, height = 6)
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
if (flag==3 | flag==0){
pdf(file = paste0(outdir, "cluster.umap.leiden.pdf"), width = 6, height = 6)
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
  saveRDS(scrna.list[[sample]], paste0(outdir, sample, ".seurat.rds"))
  scrna.list[[sample]] <- readRDS(paste0(outdir, sample, ".seurat.rds"))
}
```

# Step 8: Finding differentially expressed features
scrna.markders = list()
for (sample in samples){
  scrna.markders[[sample]] <- FindAllMarkers(scrna.list[[sample]], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.table(scrna.markders[[sample]], paste0(outdir, sample, ".FindAllMarkers.clusters.xls"), sep = "\t", col.names = NA)
}

# Heatmap of top10 (top 10 lines) marker genes
for (sample in samples){
  topN <- scrna.markders[[sample]] %>% group_by(cluster) %>% top_n(n = 10  , wt = avg_log2FC)
  DoHeatmap(scrna.list[[sample]], features = topN$gene, size = 2, draw.lines = T, angle = 45, hjust = 0.2) + theme(axis.text.y = element_text(size = 5)) + NoLegend()
  ggsave(paste0(outdir, sample, ".top10markergenes.heatmap.pdf"), width = 8, height = 6)
}

# Step 9: Feature plot
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

for (sample in samples){
  dp = DotPlot(scrna.list[[sample]], features = markers) + RotatedAxis()
  dotplot = dot_plot(dp$data[,c(3,4,1,2,5)], shape_var = "pct.exp", col_var = "avg.exp.scaled", shape_legend = "Percent Expressed", col_legend = "Average Expression", x.lab.pos = "bottom", dend_x_var = c("pct.exp","avg.exp.scaled"), dend_y_var = c("pct.exp","avg.exp.scaled"), hclust_method = "ward.D2", do.return=T, y.lab.size.factor=0.1)
#Idents(scrna.list[[sample]])=factor(Idents(scrna.list[[sample]]), levels=dotplot$plot$grobs[[4]]$label)
  p1 <- new_dotplot(scrna.list[[sample]], features = dotplot$plot$grobs[[10]]$label, group.by = "seurat_clusters", shape.scale = 6, color.breaks.values = c(-2, -1, 0, 1, 2), size.breaks.values = c(0, 25, 50, 75, 100)) + theme(axis.text = element_text(size = 12))
  p1
  pdf(paste0(outdir, sample, ".markergenes.dotPlot.clustered.pdf"), height = 6, width = 20)
  print(p1)
  dev.off()
}