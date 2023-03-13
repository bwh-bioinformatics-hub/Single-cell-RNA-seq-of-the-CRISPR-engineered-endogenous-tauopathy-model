#!/usr/bin/env Rscript

library(optparse)

options(stringsAsFactors=FALSE)

##################
# OPTION PARSING
##################


option_list <- list(
  make_option(opt_str = c("-l","--in_dir"),
              type = "character",
              default = NULL,
              help = "Input directory",
              metavar = "character"),
  make_option(opt_str = c("-s", "--in_scrublet_dir"),
              type = "character",
              default = NULL,
              help = "Input scrublet directory",
              metavar = "character"),
  make_option(opt_str = c("-k", "--out_dir"),
              type = "character",
              default = NULL,
              help = "output directory",
              metavar = "character"),
  make_option(opt_str = c("-j","--in_projectName"),
              type = "character",
              default = NULL,
              help = "Sample",
              metavar = "character"),
  make_option(opt_str = c("-r","--in_refdir"),
              type = "character",
              default = NULL,
              help = "ref",
              metavar = "character")
  )



opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)

# Load Libraries
quiet_library <- function(...) {
  suppressPackageStartupMessages(library(...))
}
quiet_library(dplyr)
quiet_library(Seurat)
quiet_library(patchwork)
quiet_library(ggplot2)
quiet_library(stringr)
quiet_library(SoupX)
quiet_library(reticulate)
quiet_library(glmpca)
quiet_library(SeuratWrappers)
quiet_library(scry)
quiet_library(reticulate)
quiet_library(monocle3)
quiet_library(FlexDotPlot)
quiet_library(cowplot)
quiet_library(tidyverse)
quiet_library(viridis)
quiet_library(scCustomize)
quiet_library(qs)
quiet_library(gridExtra)
quiet_library(plyr)
quiet_library(circlize)
quiet_library(ComplexHeatmap)
quiet_library(EnhancedVolcano)
quiet_library(H5MANIPULATOR)
quiet_library(rio)

# Arg Parse

if(is.null(args$in_dir)) {
  indir <- system.file("testdata/sample1_molecule_info.h5", package = "H5MANIPULATOR")
  scrubletdir <- system.file("testdata/sample1_metrics_summary.csv", package = "H5MANIPULATOR")
  refdir <- system.file("reference/SampleSheet_fallback.csv", package = "H5MANIPULATOR")
  projectName <- "B000-P0C0W0"
} else {
  indir <- args$in_dir
  scrubletdir <- args$in_scrublet_dir
  outdir <- args$out_dir
  refdir <- args$in_refdir
  projectName <- args$in_projectName
}

if(!file.exists(indir)) {
  stm(paste0("ERROR: Cannot find input directory:", in_dir))
  stop()
}
if(!file.exists(scrubletdir)) {
  stm(paste0("ERROR: Cannot find IN scrublet directory", in_scrublet_dir))
  stop()
}

flag=1 # 1=louvain, 2= leiden, 0=louvain and leiden

# If using Leiden algorithm in FindMarkers
# use_condaenv("r_leiden", required=TRUE)
# py_config()

# Create a vector of sample names, bad sequencing sample removed
samples = c("control1", "CRN00224919", "P251L1merged", "P251L2merged", "CTRL2merged", "LIB055588_CRN00233457")
treatments <- c("control", "control", "P251L", "P251L", "control", "P251L")

markers <-  c("VAChT", "VGlut", "Gad1", "Vmat", "ey", "prt", "sNPF", "trio", "hth", "bsh", "Eaat1", "Lim3", "svp", "eya", "Lim1", "acj6", "Crz", "ple", "alrm", "kn", "CG14687", "Poxn", "Dh31", "grn", "Dr", "Sox21b", "Pdfr", "tim", "vri", "Nplp1", "CG2016", "CG33777", "CG5910", "otk2", "lncRNA:CR43856", "Sulf1", "ort", "Pka-C3", "Tk", "CG14989", "Ets65A", "hbn", "CG13698", "Dll", "sosie", "Mmp2", "tup", "Eip63F-1", "msi", "CG10804", "bru3", "dati", "mbl", "lncRNA:CR44024", "CG9650", "CG4577", "ap", "scro", "DIP-theta", "beat-Ia", "side-IV", "CG42750", "Drgx", "CNMaR", "Sox102F", "SoxN", "CG14340", "cv-c")

# read in gene name table
geneTable <- read.csv(paste0(refdir, "geneAnnotationTable.csv"), header = T, row.names = 1)

# Clustered_DotPlot_relabel function from developer
# scCustermize code updated by developer, https://github.com/samuel-marsh/scCustomize/issues/27

Clustered_DotPlot_relabel <- function(
  seurat_object,
  features,
  new_row_labels,
  colors_use_exp = viridis_plasma_dark_high,
  exp_color_min = -2,
  exp_color_middle = NULL,
  exp_color_max = 2,
  print_exp_quantiles = FALSE,
  colors_use_idents = NULL,
  x_lab_rotate = TRUE,
  k = 1,
  row_km_repeats = 1000,
  column_km_repeats = 1000,
  row_label_size = 8,
  raster = FALSE,
  plot_km_elbow = TRUE,
  elbow_kmax = NULL,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  show_parent_dend_line = TRUE,
  ggplot_default_colors = FALSE,
  seed = 123
) {
  # Check for packages
  ComplexHeatmap_check <- PackageCheck("ComplexHeatmap", error = FALSE)
  if (!ComplexHeatmap_check[1]) {
    stop(
      "Please install the ComplexHeatmap package to use Clustered_DotPlot",
      "\nThis can be accomplished with the following commands: ",
      "\n----------------------------------------",
      "\ninstall.packages('BiocManager')",
      "\nBiocManager::install('ComplexHeatmap')",
      "\n----------------------------------------",
      call. = FALSE
    )
  }
  
  # Check Seurat
  scCustomize:::Is_Seurat(seurat_object = seurat_object)
  
  # Check unique features
  features_unique <- unique(x = features)
  
  if (length(x = features_unique) != length(x = features)) {
    warning("Feature list contains duplicates, making unique.")
  }
  
  # Check exp min/max set correctly
  if (!exp_color_min < exp_color_max) {
    stop("The value for 'exp_color_min': ", exp_color_min, ", must be less than the value for 'exp_color_max': ", exp_color_max, ".")
  }
  
  # Get DotPlot data
  seurat_plot <- DotPlot(object = seurat_object, features = features_unique, assay = assay, group.by = group.by, scale = TRUE, idents = idents, col.min = NULL, col.max = NULL)
  
  data <- seurat_plot$data
  
  # Get expression data
  exp_mat <- data %>%
    select(-pct.exp, -avg.exp) %>%
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
    as.data.frame()
  
  row.names(x = exp_mat) <- exp_mat$features.plot
  
  # Check NAs if idents
  if (!is.null(x = idents)) {
    # Find NA features and print warning
    excluded_features <- exp_mat[rowSums(is.na(x = exp_mat)) > 0,] %>%
      rownames()
    warning("The following features were removed as there is no scaled expression present in subset (`idents`) of object provided: ", glue_collapse_scCustom(input_string = excluded_features, and = TRUE), ".")
    
    # Extract good features
    good_features <- rownames(exp_mat)
    
    # Remove rows with NAs
    exp_mat <- exp_mat %>%
      filter(features.plot %in% good_features)
  }
  
  exp_mat <- exp_mat[,-1] %>%
    as.matrix()
  
  # Get percent expressed data
  percent_mat <- data %>%
    select(-avg.exp, -avg.exp.scaled) %>%
    pivot_wider(names_from = id, values_from = pct.exp) %>%
    as.data.frame()
  
  row.names(x = percent_mat) <- percent_mat$features.plot
  
  # Subset dataframe for NAs if idents so that exp_mat and percent_mat match
  if (!is.null(x = idents)) {
    percent_mat <- percent_mat %>%
      filter(features.plot %in% good_features)
  }
  
  percent_mat <- percent_mat[,-1] %>%
    as.matrix()
  
  # print quantiles
  if (print_exp_quantiles) {
    message("Quantiles of gene expression data are:")
    print(quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99)))
  }
  
  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = seurat_object)
  
  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }
  
  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use_idents) && ggplot_default_colors) {
    stop("Cannot provide both custom palette to `colors_use` and specify `ggplot_default_colors = TRUE`.")
  }
  if (is.null(x = colors_use_idents)) {
    # set default plot colors
    if (is.null(x = colors_use_idents)) {
      colors_use_idents <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
    }
  }
  
  # Pull Annotation and change colors to ComplexHeatmap compatible format
  Identity <- colnames(exp_mat)
  
  identity_colors <- DiscretePalette_scCustomize(num_colors = length(Identity), palette = "polychrome", shuffle_pal = F)
  names(identity_colors) <- Identity
  identity_colors_list <- list(Identity = identity_colors)
  
  # Create identity annotation
  column_ha <- ComplexHeatmap::HeatmapAnnotation(Identity = Identity,
                                                 col =  identity_colors_list,
                                                 na_col = "grey",
                                                 name = "Identity"
  )
  
  # Set middle of color scale if not specified
  if (is.null(x = exp_color_middle)) {
    exp_color_middle <- scCustomize:::Middle_Number(min = exp_color_min, max = exp_color_max)
  }
  
  palette_length <- length(colors_use_exp)
  palette_middle <- scCustomize:::Middle_Number(min = 0, max = palette_length)
  
  # Create palette
  col_fun = colorRamp2(c(exp_color_min, exp_color_middle, exp_color_max), colors_use_exp[c(1,palette_middle, palette_length)])
  
  # Calculate and plot Elbow
  if (plot_km_elbow) {
    # if elbow_kmax not NULL check it is usable
    if (!is.null(x = elbow_kmax) && elbow_kmax > (nrow(x = exp_mat) - 1)) {
      elbow_kmax <- nrow(x = exp_mat) - 1
      warning("The value provided for 'elbow_kmax' is too large.  Changing to (length(x = features)-1): ", elbow_kmax)
    }
    
    # if elbow_kmax is NULL set value based on input feature list
    if (is.null(x = elbow_kmax)) {
      # set to (length(x = features)-1) if less than 21 features OR to 20 if greater than 21 features
      if (nrow(x = exp_mat) > 21) {
        elbow_kmax <- 20
      } else {
        elbow_kmax <- nrow(x = exp_mat) - 1
      }
    }
    
    km_elbow_plot <- scCustomize:::kMeans_Elbow(data = exp_mat, k_max = elbow_kmax)
  }
  
  # prep heatmap
  if (raster) {
    layer_fun = function(j, i, x, y, w, h, fill) {
      grid.rect(x = x, y = y, width = w, height = h,
                gp = gpar(col = NA, fill = NA))
      grid.circle(x=x,y=y,r= sqrt(ComplexHeatmap::pindex(percent_mat, i, j)/100)  * unit(2, "mm"),
                  gp = gpar(fill = col_fun(ComplexHeatmap::pindex(exp_mat, i, j)), col = NA))
    }
  } else {
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.rect(x = x, y = y, width = w, height = h,
                gp = gpar(col = NA, fill = NA))
      grid.circle(x=x,y=y,r= sqrt(percent_mat[i, j]/100) * unit(2, "mm"),
                  gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))
    }
  }
  
  # Create legend for point size
  lgd_list = list(
    ComplexHeatmap::Legend(labels = c(0.25,0.5,0.75,1), title = "Percent Expressing",
                           graphics = list(
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                                              gp = gpar(fill = "black")))
    )
  )
  
  # Set x label roration
  if (is.numeric(x = x_lab_rotate)) {
    x_lab_rotate <- x_lab_rotate
  } else if (isTRUE(x = x_lab_rotate)) {
    x_lab_rotate <- 45
  } else {
    x_lab_rotate <- 0
  }
  
  # Create Plot
  set.seed(seed = seed)
  if (raster) {
    cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                heatmap_legend_param=list(title="Expression"),
                                                col=col_fun,
                                                rect_gp = gpar(type = "none"),
                                                layer_fun = layer_fun,
                                                row_names_gp = gpar(fontsize = row_label_size),
                                                row_km = k,
                                                row_km_repeats = row_km_repeats,
                                                border = "black",
                                                top_annotation = column_ha,
                                                column_km_repeats = column_km_repeats,
                                                show_parent_dend_line = show_parent_dend_line,
                                                column_names_rot = x_lab_rotate)
  } else {
    cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                heatmap_legend_param=list(title="Expression"),
                                                col=col_fun,
                                                rect_gp = gpar(type = "none"),
                                                cell_fun = cell_fun,
                                                row_names_gp = gpar(fontsize = row_label_size),
                                                row_km = k,
                                                row_km_repeats = row_km_repeats,
                                                border = "black",
                                                top_annotation = column_ha,
                                                column_km_repeats = column_km_repeats,
                                                show_parent_dend_line = show_parent_dend_line,
                                                column_names_rot = x_lab_rotate)
  }
  
  # Add pt.size legend & return plots
  if (plot_km_elbow) {
    return(list(km_elbow_plot, ComplexHeatmap::draw(cluster_dot_plot, annotation_legend_list = lgd_list)))
  }
  return(ComplexHeatmap::draw(cluster_dot_plot + rowAnnotation(rn= anno_text(new_row_labels)), annotation_legend_list = lgd_list))
}

# Step 2: Pre-processing
# Remove ambient RNA by SoupX
data.10x = list()
for (sample in samples){
  filt.matrix <- Read10X_h5(paste0(indir, sample, "/outs/filtered_feature_bc_matrix.h5"), use.names = F)
  raw.matrix <- Read10X_h5(paste0(indir, sample, "/outs/raw_feature_bc_matrix.h5"), use.names = F)
  srat <- CreateSeuratObject(counts = filt.matrix)
  soup.channel <- SoupChannel(raw.matrix, filt.matrix)
  srat <- SCTransform(srat, verbose = F)
  srat <- RunPCA(srat, verbose = F)
  srat <- RunUMAP(srat, dims = 1:30, verbose = F)
  srat <- FindNeighbors(srat, dims = 1:30, verbose = F)
  srat <- FindClusters(srat, verbose = T)
  meta <- srat@meta.data
  umap <- srat@reductions$umap@cell.embeddings
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- setDR(soup.channel, umap)
  soup.channel <- autoEstCont(soup.channel)
  data.10x[[sample]] <- adjustCounts(soup.channel, roundToInt = T)
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
#system2(command = "bash", args = c("run_scrublet_multi.sh"))

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

# add sample name
for(i in 1:length(samples)){
  sample=samples[i]; treatment=treatments[i];
  scrna.list[[sample]]$treatment <- treatment
}

rm(ds)
rm(filt.matrix)
rm(meta)
rm(raw.matrix)
rm(soup.channel)
rm(srat)
rm(doublet_scores)
rm(predicted_doublets)

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

# Step 4: Integration
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
rm(scrna.anchors)

# Perform an integrated analysis
DefaultAssay(scrna.combined) <- "integrated"
scrna.combined <- ScaleData(scrna.combined, verbose = FALSE)
scrna.combined <- RunPCA(scrna.combined, npcs = 15, verbose = FALSE)

# elbow plot
pdf(file = paste0(outdir, "elbow.plot.pdf"), width = 8, height = 8)
p1 <- ElbowPlot(scrna.combined) + ggtitle("Integrated") + theme(aspect.ratio=5/10) + theme(plot.margin = unit(c(3, 3, 3, 3), "cm"))
print(p1 + coord_fixed())
dev.off()

# Continue on integrated analysis
scrna.combined <- RunUMAP(scrna.combined, reduction = "pca", dims = 1:15)
scrna.combined <- FindNeighbors(scrna.combined, reduction = "pca", dims = 1:15)
scrna.combined <- FindClusters(scrna.combined, resolution = 0.5)

rm(scrna.list)

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

rm(p1)
rm(p2)
rm(p3)
rm(p4)

# Step 5: Identify conserved cell type markers
# For performing differential expression after integration, we switch back to the original data
DefaultAssay(scrna.combined) <- "RNA"
scrna.combined <- ScaleData(scrna.combined, verbose = FALSE)
nk.markers <- FindConservedMarkers(scrna.combined, ident.1 = 0, grouping.var = "orig.ident", verbose = FALSE)
head(nk.markers)

# FindAllMarkers
scrna.markers <- FindAllMarkers(scrna.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
names(scrna.markers)[names(scrna.markers) == "gene"] <- "geneID"

scrna.markers <- cbind(scrna.markers, geneSymbol=geneTable$geneSymbol[match(scrna.markers$geneID, geneTable$geneID)])
write.table(scrna.markers, paste0(outdir, "FindAllMarkers.clusters.xls"), sep = "\t", row.names = F)
scrna.markers.wide <- reshape(scrna.markers, idvar = c("geneID", "geneSymbol"), timevar = "cluster", direction = "wide")
write.table(scrna.markers.wide, paste0(outdir, "FindAllMarkers.clusters.wide.xls"), sep = "\t", row.names = F)

topN <- scrna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(topN, paste0(outdir, "FindAllMarkers.clusters.top10.xls"), sep = "\t", col.names = NA)

DefaultAssay(scrna.combined) <- "integrated"
DoHeatmap(scrna.combined, features = topN$geneID, size = 2, draw.lines = T, angle = 45, hjust = 0.2) + theme(axis.text.y = element_text(size = 5)) + NoLegend()
ggsave(paste0(outdir, "top10markers.heatmap.integrated.geneID.pdf"), width = 8, height = 12)

DoHeatmap(scrna.combined, features = topN$geneID, size = 2, draw.lines = T, angle = 45, hjust = 0.2) + theme(axis.text.y = element_text(size = 5)) + NoLegend() + scale_y_discrete(breaks=topN$geneID,
        labels=geneTable$geneSymbol[match(topN$geneID, geneTable$geneID)])
ggsave(paste0(outdir, "top10markers.heatmap.integrated.geneSymbol.pdf"), width = 8, height = 12)

DefaultAssay(scrna.combined) <- "RNA"
DoHeatmap(scrna.combined, features = topN$geneID, size = 2, draw.lines = T, angle = 45, hjust = 0.2) + theme(axis.text.y = element_text(size = 5)) + NoLegend()
ggsave(paste0(outdir, "top10markers.heatmap.RNA.geneID.pdf"), width = 8, height = 12)

DoHeatmap(scrna.combined, features = topN$geneID, size = 2, draw.lines = T, angle = 45, hjust = 0.2) + theme(axis.text.y = element_text(size = 5)) + NoLegend() + scale_y_discrete(breaks=topN$geneID,
        labels=geneTable$geneSymbol[match(topN$geneID, geneTable$geneID)])
ggsave(paste0(outdir, "top10markers.heatmap.RNA.geneSymbol.pdf"), width = 8, height = 12)

# Step 6: Top 3 identified genes, feature plot, dotplot
topN <- Extract_Top_Markers(scrna.markers, num_genes = 3, named_vector = FALSE, make_unique = TRUE, gene_column = "geneID")

# Feature plot
pdf(paste0(outdir, "combined.top3markers.geneID.pdf"))
ggp = list()
for (marker in topN){
    ggp[[marker]]=FeaturePlot(scrna.combined, features=marker)
    print(ggp[[marker]])
}
dev.off()

pdf(paste0(outdir, "combined.top3markers.geneSymbol.pdf"))
ggp = list()
for (marker in topN){
    ggp[[marker]]=FeaturePlot(scrna.combined, features=marker) + ggtitle(geneTable$geneSymbol[match(marker, geneTable$geneID)])
    print(ggp[[marker]])
}
dev.off()

# Dotplot
pdf(paste0(outdir, "combined.dotplot.DEtop3.geneID.pdf"), width = 20, height = 10)
p1 <- DotPlot_scCustom(scrna.combined, features = topN, x_lab_rotate = TRUE, colors_use = "blue")
print(p1)
dev.off()

pdf(paste0(outdir, "combined.dotplot.DEtop3.geneSymbol.pdf"), width = 20, height = 10)
p1 <- DotPlot_scCustom(scrna.combined, features = topN, x_lab_rotate = TRUE, colors_use = "blue") + scale_x_discrete(breaks=topN, labels=geneTable$geneSymbol[match(topN, geneTable$geneID)])
print(p1)
dev.off()

# developer, clustered
pdf(paste0(outdir, "combined.dotplot.DEtop3.clustered.geneSymbol.pdf"), width = 10, height = 15)
p1 <- Clustered_DotPlot_relabel(scrna.combined, features = topN, plot_km_elbow = F, new_row_labels = geneTable$geneSymbol[match(topN, geneTable$geneID)])
print(p1)
dev.off()

pdf(paste0(outdir, "combined.dotplot.DEtop3.clustered.geneID.pdf"), width = 10, height = 15)
p1 <- Clustered_DotPlot_relabel(scrna.combined, features = topN, x_lab_rotate = F, plot_km_elbow = FALSE, new_row_labels = topN)
print(p1)
dev.off()

# Step 7: Customer markers, feature plot, dotplot
# Feature plot
markers.ID = geneTable$geneID[match(markers, geneTable$geneSymbol)]

pdf(paste0(outdir, "combined.markers.geneID.pdf"))
ggp = list()
for (marker in markers.ID){
    ggp[[marker]]=FeaturePlot(scrna.combined, features=marker)
    print(ggp[[marker]])
}
dev.off()

pdf(paste0(outdir, "combined.markers.geneSymbol.pdf"))
ggp = list()
for (marker in markers.ID){
    ggp[[marker]]=FeaturePlot(scrna.combined, features=marker) + ggtitle(geneTable$geneSymbol[match(marker, geneTable$geneID)])
    print(ggp[[marker]])
}
dev.off()

# Dotplot
pdf(paste0(outdir, "combined.dotplot.geneID.pdf"), width = 30, height = 10)
p1 <- DotPlot_scCustom(scrna.combined, features = markers.ID, x_lab_rotate = TRUE, colors_use = "blue")
print(p1)
dev.off()

pdf(paste0(outdir, "combined.dotplot.geneSymbol.pdf"), width = 30, height = 10)
p1 <- DotPlot_scCustom(scrna.combined, features = markers.ID, x_lab_rotate = TRUE, colors_use = "blue") + scale_x_discrete(breaks=markers.ID, labels=geneTable$geneSymbol[match(markers.ID, geneTable$geneID)])
print(p1)
dev.off()

# developer, clustered
pdf(paste0(outdir, "combined.dotplot.clustered.geneSymbol.pdf"), width = 10, height = 15)
p1 <- Clustered_DotPlot_relabel(scrna.combined, features = markers.ID, plot_km_elbow = F, new_row_labels = geneTable$geneSymbol[match(markers.ID, geneTable$geneID)])
print(p1)
dev.off()

pdf(paste0(outdir, "combined.dotplot.clustered.geneID.pdf"), width = 10, height = 15)
p1 <- Clustered_DotPlot_relabel(scrna.combined, features = markers.ID, x_lab_rotate = F, plot_km_elbow = FALSE, new_row_labels = markers.ID)
print(p1)
dev.off()

# Step 7: Save the data
saveRDS(scrna.combined, paste0(outdir, "scrna.combined.seurat.", projectName, ".rds"))
#scrna.combined <- readRDS(paste0(outdir, "scrna.combined.seurat.", projectName, ".rds"))




# Cell type annotation. CHANGE HERE!!!
scrna.combined.beforeAnnotation <- scrna.combined
#scrna.combined <- RenameIdents(scrna.combined, `0` = "Mip_Serotonergic_Clock Neurons", `1` = "Tm20", `2` = "Photoreceptor like cells", `3` = "TmY5a/TmY14", `4` = "DM3", `5` = "PM1/Pm2", `6` = "T4/T5", `7` = "T2/T3", `8` = "Astrocyte-Like/Ensheating glia", `9` = "Gamma Kenyon cell", `10` = "Dopaminergic PAM neurons", `11` = "Mi1", `12` = "C2/C3 neurons", `13` = "Adult fat body/Glial cells", `14` = "Transmedullary neuron Tm9", `15` = "Alpha/beta Kenyon cell", `16` = "Mi15/Tmy8/T neurons", `17` = "Transmedullary Neurons", `18` = "MBON/Hemocytes", `19` = "LC12/LC17", `20` = "Mi15", `21` = "Columnar neuron T1", `22` = "Dm8", `23` = "Alpha/beta Kenyon cell", `24` = "Tm2/L neurons", `25` = "Lawf 1 and 2", `26` = "Dm8/Tm5c", `27` = "Tm1/Tm4/Tm20", `28` = "Pox neuron", `29` = "Lai/Tm1")


scrna.combined <- RenameIdents(scrna.combined, `0` = "Unannotated", `1` = "Tm20", `2` = "Unannotated", `3` = "TmY14", 
`4` = "Pm1-Pm3", `5` = "Dm3", `6` = "T4_T5", `7` = "Astrocytes", 
`8` = "T2_T3", `9` = "Tm3a", `10` = "G_KC", `11` = "Tmy8", `12` = "MBON", `13` = "C2_C3", `14` = "Mi1", 
`15` = "AB_kc", `16` = "Dm8-TM5c", `17` = "Mi4", `18` = "Tm9", 
`19` = "T1", `20` = "Mi15", `21` = "A_B_KC", `22` = "LC12_LC17", `23` = "Lawf1_2", `24` = "Perineurial glia",
`25` = "Dm9", `26` = "L1-L5", `27` = "Pox neurons", `28` = "Unannotated")


# UMAP, cell annotated
p1 <- DimPlot(scrna.combined, label = TRUE, repel = TRUE, pt.size = 0.5) + guides(color = guide_legend(override.aes = list(size=1), ncol=1) ) + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10))
p1
pdf(file=paste0(outdir, "combined.umap.annotated.pdf"))
print(p1 + coord_fixed())
dev.off()

p1 <- DimPlot(scrna.combined.beforeAnnotation, label = F, repel = TRUE, pt.size = 0.5) + guides(color = guide_legend(override.aes = list(size=1), ncol=1) ) + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10))
p1
pdf(file=paste0(outdir, "combined.umap.withoutAnnotation.pdf"))
print(p1 + coord_fixed())
dev.off()

# dotplot after annotation
#Idents(scrna.combined) <- factor(Idents(scrna.combined), levels = c("Mip_Serotonergic_Clock Neurons", "Tm20", "Photoreceptor like cells", "TmY5a/TmY14", "DM3", "PM1/Pm2", "T4/T5", "T2/T3", "Astrocyte-Like/Ensheating glia", "Gamma Kenyon cell", "Dopaminergic PAM neurons", "Mi1", "C2/C3 neurons", "Adult fat body/Glial cells", "Transmedullary neuron Tm9", "Alpha/beta Kenyon cell", "Mi15/Tmy8/T neurons", "Transmedullary Neurons", "MBON/Hemocytes", "LC12/LC17", "Mi15", "Columnar neuron T1", "Dm8", "Alpha/beta Kenyon cell", "Tm2/L neurons", "Lawf 1 and 2", "Dm8/Tm5c", "Tm1/Tm4/Tm20", "Pox neuron", "Lai/Tm1"))

#Idents(scrna.combined) <- factor(Idents(scrna.combined), levels = c("Mip_Serotonergic_Clock", "Tm20", "Photoreceptor like", "TmY5a_TmY14", "DM3", "PM1_Pm2", "T4_T5", "T2_T3", "Astrocyte-Like_Ensheating glia", "Gamma Kenyon cell", "Dopaminergic PAM neurons", "Mi1", "C2_C3 neurons", "Adult fat body_Glial cells", "Transmedullary neuron Tm9", "Alpha_beta Kenyon cell", "Mi15_Tmy8_T", "Transmedullary Neurons", "MBON_Hemocytes", "LC12_LC17", "Mi15", "Columnar neuron T1", "Dm8", "Alpha_beta Kenyon cell", "Tm2_L neurons", "Lawf 1 and 2", "Dm8_Tm5c", "Tm1_Tm4_Tm20", "Pox neuron", "Lai_Tm1"))

pdf(paste0(outdir, "combined.dotplot.afterAnnotation.clustered.geneSymbol.pdf"), width = 10, height = 15)
p1 <- Clustered_DotPlot_relabel(scrna.combined, features = markers.ID, plot_km_elbow = F, new_row_labels = geneTable$geneSymbol[match(markers.ID, geneTable$geneID)])
print(p1)
dev.off()

pdf(paste0(outdir, "combined.dotplot.afterAnnotation.clustered.geneID.pdf"), width = 10, height = 15)
p1 <- Clustered_DotPlot_relabel(scrna.combined, features = markers.ID, plot_km_elbow = F, new_row_labels = markers.ID)
print(p1)
dev.off()


# Step 11: Identify differentially expressed genes across treatments
theme_set(theme_cowplot())
scrna.combined$celltype <- Idents(scrna.combined)
scrna.combined$celltype.treatment <- paste(Idents(scrna.combined), scrna.combined$treatment, sep = "_")
list_cell <- levels(Idents(scrna.combined))

avg.list <- list()
for (cell in list_cell){
  g.cells <- subset(scrna.combined, idents = cell)
  Idents(g.cells) <- "treatment"
  avg.g.cells <- as.data.frame(log1p(AverageExpression(g.cells, verbose = FALSE)$RNA))
  avg.g.cells$geneID <- rownames(avg.g.cells)
  avg.g.cells <- avg.g.cells %>% select(geneID, everything())
  row.names(avg.g.cells) <- NULL
  colnames(avg.g.cells) <- c("geneID", "control_normalizedReadsCounts", "P251L_normalizedReadsCounts")
  avg.list[[cell]] <- avg.g.cells
}
openxlsx::write.xlsx(avg.list, paste0(outdir, "avgExpr.treatment.xls"))

Idents(scrna.combined) <- "celltype.treatment"
list_treament <- levels(Idents(scrna.combined))



de.list <- list()
pdf(paste0(outdir, "DE.P251L.vs.control.volcano.geneSymbol.pdf"))
for (cell in list_cell){
  mutant.responseA <- FindMarkers(scrna.combined, ident.1 = paste0(cell, "_P251L"), ident.2 = paste0(cell, "_control"), verbose = FALSE, test.use = "MAST")
  mutant.responseA$geneSymbol <- geneTable$geneSymbol[match(rownames(mutant.responseA), geneTable$geneID)]
  mutant.responseA$geneID <- rownames(mutant.responseA)
  mutant.responseA$p_val <- mutant.responseA$p_val+1e-200
  p1 <- EnhancedVolcano(mutant.responseA, lab = geneTable$geneSymbol[match(rownames(mutant.responseA), geneTable$geneID)], x = 'avg_log2FC', y = 'p_val', title = cell, FCcutoff = 0.25)
  print(p1)
  mutant.responseA <- mutant.responseA %>% select(geneID, geneSymbol, everything())
  row.names(mutant.responseA) <- NULL
  de.list[[cell]] <- mutant.responseA
} 
dev.off()
openxlsx::write.xlsx(de.list, paste0(outdir, "DE.P251L.vs.control.xls"))

pdf(paste0(outdir, "DE.P251L.vs.control.volcano.geneID.pdf"))
for (cell in list_cell){
  mutant.responseA <- FindMarkers(scrna.combined, ident.1 = paste0(cell, "_P251L"), ident.2 = paste0(cell, "_control"), verbose = FALSE, test.use = "MAST")
  mutant.responseA$p_val <- mutant.responseA$p_val+1e-200
  p1 <- EnhancedVolcano(mutant.responseA, lab = rownames(mutant.responseA), x = 'avg_log2FC', y = 'p_val', title = cell, FCcutoff = 0.25)
  print(p1)
} 
dev.off()

writeLines(capture.output(sessionInfo()), paste0(outdir, "session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
