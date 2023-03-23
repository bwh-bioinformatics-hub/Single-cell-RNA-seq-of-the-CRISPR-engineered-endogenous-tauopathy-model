### CODE FROM HASSAN 

#load libraries
library(mclust)
library(Seurat)
library(Matrix)
library(tidyverse)
library(viridis)
library(broom)
library(SingleCellExperiment)
library(slingshot)
library (ggplot2)
library(RColorBrewer)
library(scater)
library(TSCAN)
library(gam)
library(ggpubr)
library(clusterExperiment)
library(RColorBrewer)
library(tradeSeq)

#Data loading and wrangling
integrated = readRDS("/scrna.combined.seurat.hassan_2022.rds")
                     
#Simplify the cell types
integrated_celltype = integrated@meta.data$seurat_clusters

simple_celltype = case_when(
  integrated_celltype %in% c("23", "26") ~ "laminar",
  integrated_celltype %in% c("1","3","4","5","9","11","13","14","16","17","18","20","25") ~ "medullary",
  integrated_celltype %in% c("6", "8", "19", "22") ~ "lobular",
  integrated_celltype %in% c("10", "15", "21") ~ "Kenyon",
  integrated_celltype %in% c("12") ~ "MBON",
  integrated_celltype %in% c("7", "24") ~ "Glia",
  TRUE ~ "other")
integrated@meta.data$simple_celltype = simple_celltype
#convert into single cell experiemnt
sce = as.SingleCellExperiment(integrated)
# add dimension reduction embeddings to sce 
reducedDim(sce,"PCA",withDimnames=TRUE) <- integrated[['pca']]@cell.embeddings
reducedDim(sce,"UMAP",withDimnames=TRUE) <- integrated[['umap']]@cell.embeddings
# Run Slingshot on SCE
sl1 = slingshot(sce, clusterLabels="simple_celltype", reducedDim="UMAP", start.clus="MBON")
pal = brewer.pal(9, 'Set1')
pointcolors = pal[as.factor(sce$simple_celltype)]
# Plot the reduced dims of the slingshot dataset
pdf('/trajectory_output/trajectory/reduced_dims_plot.pdf')
plot(reducedDims(sce)$UMAP, col = pointcolors, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sl1), lwd=2, type = 'curves', col = 'black')
legendlabels = unique(as.factor(sce$simple_celltype))
legendcolors = unique(pointcolors)
legend("topleft", pch=16, legend=legendlabels, col=legendcolors)
dev.off()

#To get the pseudotime curves first get the embedded curves.
embedded <- embedCurves(sl1, "UMAP")
saveRDS(embedded, "/trajectory_output/trajectory/embedded.rds")


first calculate for the pseudotime 1
embedded1 <- slingCurves(embedded)[[1]] # only 1st path.
embedded1 <- data.frame(embedded1$s[embedded1$ord,])

pdf('//trajectory_output/trajectory/slingPseudotime_1.pdf')
plotUMAP(sl1, colour_by="slingPseudotime_1") +
  geom_path(data=embedded1, aes(x=UMAP_1, y=UMAP_2, cex =0.25), size=0.5)
#the output is with black curves. I draw red curves on illustrator and set the size to 0 to remove black curves for the naive picture.
dev.off()

#then calculate for the pseudotime 2
embedded2 <- slingCurves(embedded)[[2]] # for the 2nd path.
embedded2 <- data.frame(embedded2$s[embedded2$ord,])

pdf('/trajectory_output/trajectory/slingPseudotime_2.pdf')
plotUMAP(sl1, colour_by="slingPseudotime_2") +
  geom_path(data=embedded2, aes(x=UMAP_1, y=UMAP_2, cex =0.25), size=0.5)
#the output is with black curves. I draw red curves on illustrator and set the size to 0 to remove black curves for the naive picture.
dev.off()

#then calculate for the pseudotime 3
embedded3 <- slingCurves(embedded)[[3]] # for the 3rd path.
embedded3 <- data.frame(embedded3$s[embedded3$ord,])

pdf('/trajectory_output/trajectory/slingPseudotime_3.pdf')
plotUMAP(sl1, colour_by="slingPseudotime_3") +
  geom_path(data=embedded3, aes(x=UMAP_1, y=UMAP_2, cex =0.25), size=0.5)
#the output is with black curves. I draw red curves on illustrator and set the size to 0 to remove black curves for the naive picture.
dev.off()



#After calculating the pseudotimes, calculate the differntial genes by using the gam packages along each pseudotime lineage by using the following code.

# Only look at the 1,000 most variable genes when identifying temporally expressesd genes.
# Identify the variable genes by ranking all genes by their variance.
Y <- log2(logcounts(sl1) + 1)
var1K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:1000]
Y <- Y[var1K, ]  # only counts for variable genes

# Fit GAM for each gene using pseudotime as independent variable.
#simply replace sling Pseudtime for other pseudotimes and then save the gam.pval object for each pseudotimes
#gam.pval is saved for pseudotime 1, gam.pavl2 for 1, gam.pval3 for 3 and gam.pval4 for 4.
#it takes a lot of time to calculate these objects, ~ 10-15 hours.
t <- sl1$slingPseudotime_3
gam.pval3 <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})
saveRDS(gam.pval3,'/trajectory/hassan_traj/')

gam.pval3 <- readRDS("/home/acicalo/bioinf_hub/pipeline/hassan_2022/trajectory_output/trajectory/hassan_traj/gam.pval3.rds")
topgenes3 <- names(sort(gam.pval3, decreasing = FALSE))[1:100]  
write.csv(topgenes3,file='/home/acicalo/bioinf_hub/pipeline/hassan_2022/trajectory_output/trajectory/hassan_traj/topgenes3.csv', row.names=FALSE)

t <- sl1$slingPseudotime_3

# Prepare and plot a heatmap of the top genes that vary their expression over pseudotime.
require(clusterExperiment)
heatdata <- as.matrix(integrated@assays$RNA[rownames(integrated@assays$RNA) %in% topgenes3, order(t, na.last = NA)])
#since the original seurat object has flybase gene IDs. We will do a bit of data wrangling to replace the gene IDs with gene names.
#rownames1 <- rownames(heatdata)
#rownames1 <- write.csv(rownames1,file='/Volumes/Macintosh HD/Users/hassanbukhari/Desktop/Tau P251L Knock In Paper/Trajectory analysis/rownames1.csv')
#rownames2 <- read.csv("../Trajectory analysis/Heatmap1 annotated-1.txt")
#rownames3 <- rownames2 [, "x"]
#rownames(heatdata) <- rownames3
heatclus <- integrated@meta.data$simple_celltype[order(t, na.last = NA)]

# convert names of slingshot object from flybase id to gene symbol 
refdir <- '/mnt/data0/referenceGenome/Drosophila_Melanogaster/Ensembl/dm6/'
geneTable <- read.csv(paste0(refdir, "geneAnnotationTable.csv"), header = T, row.names = 1)
geneSymbols <-geneTable$geneSymbol[match(row.names(heatdata), geneTable$geneID)]
row.names(heatdata) <- geneSymbols

#png(paste0("/Volumes/Macintosh HD/Users/hassanbukhari/Desktop/Tau P251L Knock In Paper/Trajectory analysis/", "heatmap_time_genes-lineage1-Gene names.png"), width=10, height=10, units = "in", res=200)
png(paste0("trajectory_output/trajectory/hassan_traj/", "heatmap_time_genes-lineage3-Gene_names.png"), width=10, height=10, units = "in", res=200)
ce <- ClusterExperiment(heatdata, heatclus, transformation = log1p)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed', cexRow = 1.5, fontsize = 15)
dev.off()



# read in all of the top genes
topgenes1 <- read.csv('/trajectory_output/trajectory/hassan_traj/topgenes1.csv')
topgenes2 <- read.csv('/trajectory_output/trajectory/hassan_traj/topgenes2.csv')
topgenes2$geneSymbol <-geneTable$geneSymbol[match(topgenes2$x, geneTable$geneID)]

topgenes3 <- read.csv('/trajectory_output/trajectory/hassan_traj/topgenes3.csv')
topgenes3$geneSymbol <-geneTable$geneSymbol[match(topgenes3$x, geneTable$geneID)]

genes_in_common <- intersect(intersect(topgenes1$gene_symbol,topgenes2$geneSymbol),topgenes3$geneSymbol)
write.csv(genes_in_common,'/trajectory/trajectory_hassan_script/genesincommon_3lineages.csv')
geneID <-geneTable$geneID[match(genes_in_common, geneTable$geneSymbol)]


genes <- c(topgenes1$gene_symbol,topgenes2$geneSymbol,topgenes3$geneSymbol)

genes <- unique(genes)
write.csv(genes,'/trajectory/trajectory_output/topgenes_3lineages.csv',row.names = F)

# pull in biomaRt to take a look at these genes in common
library(biomaRt)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                           dataset = "dmelanogaster_gene_ensembl", 
                           host = "www.ensembl.org")

attr <- listAttributes(mart)


geneID <-geneTable$geneID[match(genes, geneTable$geneSymbol)]
go_list <- getBM(attributes=c("ensembl_gene_id","go_id", "name_1006", "namespace_1003"),
              filters = "ensembl_gene_id",
              values = geneID,
              mart=mart)
head(go_list)
go_list$geneSymbol <- geneTable$geneSymbol[match(go_list$ensembl_gene_id, geneTable$geneID)]
write.csv(go_list,'/home/acicalo/mydropbox/hassan2021/TOHassan_022023/trajectory/trajectory_hassan_script/genes_annotated_w_biomaRt.csv')
#After looking at top 100 differential genes in each lineage and the ontologies, select the genes to plot on the pseudotime lineages.
# the gene list in the top plot comes from the ontology term associative learning.
gene <- sce@rowRanges
cell <- sce@colData$simple_celltype
toplot = geneID
markerdat = logcounts(sce) %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  gather("cell", "logcount", -gene) %>%
  filter(gene %in% toplot)

sshot = SlingshotDataSet(sl1)

reducedDims <- reducedDims(sl1)

pcamarkers = data.frame(pseudotime1=sl1$slingPseudotime_1,
                        pseudotime2=sl1$slingPseudotime_2,
                        pseudotime3=sl1$slingPseudotime_3,
                        UMAP1=reducedDims$UMAP[, "UMAP_1"],
                        UMAP2=reducedDims$UMAP[, "UMAP_2"],
                        cell=colnames(sl1)) %>%
  left_join(markerdat, by="cell")



pcamarkers = pcamarkers %>%
  gather(c(pseudotime1, pseudotime2, pseudotime3), key="lineage", value="pseudotime")
pcamarkers = pcamarkers %>%
  mutate(lineagestring=ifelse(lineage == "pseudotime1", paste0(sshot@lineages$Lineage1, collapse="->"),
                              ifelse(lineage == "pseudotime2", paste0(sshot@lineages$Lineage2, collapse="->"),
                                            paste0(sshot@lineages$Lineage3, collapse="->"))))


# This makes the plot of gene expression as a function of pseudotime
png(paste0("/home/acicalo/bioinf_hub/pipeline/hassan_2022/trajectory_output/trajectory/hassan_traj/", ".png"), width=10, height=10, units = "in", res=200)
ggplot(pcamarkers, aes(pseudotime, logcount, group=gene, color=gene)) +
  geom_smooth(se=FALSE, lwd=0.5) +
  facet_wrap(~lineagestring, ncol=1)
dev.off()
png(paste0("/home/acicalo/bioinf_hub/pipeline/hassan_2022/trajectory_output/trajectory/hassan_traj/", ".png"), width=10, height=10, units = "in", res=200)
ggplot(pcamarkers, aes(pseudotime, logcount, color=lineagestring)) +
  geom_smooth(se=FALSE) +
  facet_wrap(~gene, ncol=1, strip.position="left", scale='free_y') +
  theme_pubr() +
  theme(panel.margin=unit(0.5, "lines"),
        strip.text.y = element_text(angle = 180),
        axis.text.y=element_text(size=8)) +
  scale_color_viridis(discrete=TRUE, begin=0) +
  ylab("log count") +
  labs(color="lineage") 
dev.off()
