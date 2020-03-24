

source('0_my_functions.R')
library(Seurat)
library(Antler)
library(cowplot)
library(clustree)
library(dplyr)
library(gridExtra)
library(grid)
library(pheatmap)
library(tidyverse)



data.dir = "~/dev/data/10x/cellranger_output/"
plot.path = "~/dev/output/10x_neural_tube/plots/1_Seurat_full/"
rds.path = '~/dev/output/10x_neural_tube/rds.files/1_Seurat_full/'

###### load data ##########
sample.paths<-data.frame(tissue = c("hh4", "hh6", "ss4", "ss8"),
                         path = c(paste0(data.dir, "cellranger_count_hh4/outs/filtered_feature_bc_matrix_chr_edit/"),
                                  paste0(data.dir, "cellranger_count_hh6/outs/filtered_feature_bc_matrix_chr_edit/"),
                                  paste0(data.dir, "cellranger_count_4ss/outs/filtered_feature_bc_matrix_chr_edit/"),
                                  paste0(data.dir, "cellranger_count_8ss/outs/filtered_feature_bc_matrix_chr_edit/")))

# Make Seurat objects for each of the different samples. The raw data for each sample is found in the relative directory assigned above in sample.paths.
for(i in 1:nrow(sample.paths["path"])){
  name<-paste(sample.paths[i,"tissue"])
  assign(name, CreateSeuratObject(counts= Read10X(data.dir = paste(sample.paths[i,"path"])), project = paste(sample.paths[i, "tissue"])))
}

# The four Seurat objects are then merged, before running CreateSeuratObject again on the output in order to apply the min.cells parameter on the final merged dataset.
temp <- merge(hh4, y = c(hh6, ss4, ss8), add.cell.ids = c("hh4", "hh6", "ss4", "ss8"), project = "chick.10x")
merged.data<-CreateSeuratObject(GetAssayData(temp), min.cells = 3, project = "chick.10x.mincells3")

# The original Seurat objects are then removed from the global environment
rm(hh4, hh6, ss4, ss8, sample.paths, temp)

# store mitochondrial percentage in object meta data
merged.data <- PercentageFeatureSet(merged.data, pattern = "^MT-", col.name = "percent.mt")


#####################################################################################################
#                           Filter data based on variable threshold                                 #
#####################################################################################################

# remove data which do not pass filter threshold
merged.data <- subset(merged.data, subset = c(nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 15))
# log normalize data and find variable features
norm.data <- NormalizeData(merged.data, normalization.method = "LogNormalize", scale.factor = 10000)
norm.data <- FindVariableFeatures(norm.data, selection.method = "vst", nfeatures = 2000)
# scale data and regress out MT content
norm.data <- ScaleData(norm.data, features = rownames(norm.data), vars.to.regress = "percent.mt")

# save RDS after scaling as this step takes time
saveRDS(norm.data, paste0(rds.path, "norm.data.RDS"))

#####################################################################################################
#                    Perform dimensionality reduction by PCA and UMAP embedding                    #
#####################################################################################################

# read in RDS data if needed
# norm.data <- readRDS(paste0(rds.path, "norm.data.RDS"))

curr.plot.path <- paste0(plot.path, "0_filt_data/")
dir.create(curr.plot.path)

# Run PCA analysis on the each set of data
norm.data <- RunPCA(object = norm.data, verbose = FALSE)
# Seurat's clustering algorithm is based on principle components, so we need to ensure that only the informative PCs are kept!
pdf(paste0(curr.plot.path, "dimHM.pdf"),width=15,height=25)
DimHeatmap(norm.data, dims = 1:30, balanced = TRUE, cells = 500)
dev.off()
# another heuristic method is ElbowPlot which ranks PCs based on the % variance explained by each PC
pdf(paste0(curr.plot.path, "elbowplot.pdf"),width=12,height=10)
print(ElbowPlot(norm.data, ndims = 40))
dev.off()

# Run clustering and UMAP at different PCA cutoffs - save this output to compare the optimal number of PCs to be used
PCA.level.comparison(norm.data, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
# Use PCA=15 as elbow plot is relatively stable across stages
# Use clustering resolution = 0.5 for filtering
norm.data <- FindNeighbors(norm.data, dims = 1:15, verbose = FALSE)
norm.data <- RunUMAP(norm.data, dims = 1:15, verbose = FALSE)
norm.data <- FindClusters(norm.data, resolution = 0.5, verbose = FALSE)

# plot UMAP for clusters and developmental stage
clust.stage.plot(norm.data)

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.data, only.pos = T, logfc.threshold = 0.25)
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
tenx_pheatmap(data = norm.data, metadata = c("orig.ident", "seurat_clusters"), used.genes = unique(top15$gene), width = 40, height = 20, main = "")


#####################################################################################################
#                                   Check for sex differences                                       #
#####################################################################################################
curr.plot.path <- paste0(plot.path, "1_sex_filt/")

dir.create(curr.plot.path)

# There is a strong sex effect - this plot shows DE genes between clusters 1 and 2 which are hh4 clusters. Clustering is driven by sex genes
tenx_pheatmap(data = norm.data[,rownames(norm.data@meta.data[norm.data$seurat_clusters == 1 | norm.data$seurat_clusters == 2,])],
              metadata = c("orig.ident", "seurat_clusters"), used.genes = rownames(FindMarkers(norm.data, ident.1 = 1, ident.2 = 2)),
              width = 40, height = 20, main = "", hclust_rows = T, basename = "sex.effect")

# Use W chromosome genes to K-means cluster the cells into male (zz) and female (zw)
W_genes <- as.matrix(norm.data@assays$RNA[grepl("W-", rownames(norm.data@assays$RNA)),])
k_clusters <- kmeans(t(W_genes), 2)
k_clusters <- data.frame(k_clusters$cluster)
norm.data@meta.data$k_clusters <- k_clusters[match(colnames(norm.data@assays$RNA), rownames(k_clusters)),]

# get rownames for kmeans clusters 1 and 2
k_clus_1 <- rownames(norm.data@meta.data[norm.data@meta.data$k_clusters == 1,])
k_clus_2 <- rownames(norm.data@meta.data[norm.data@meta.data$k_clusters == 2,])

# work out if cluster 1 is male or female
# sum of W genes is order of magnitude greater in cluster 2 - these are the female cells
sumclus1 <- sum(W_genes[,k_clus_1])
sumclus2 <- sum(W_genes[,k_clus_2])

if(sumclus1 < sumclus2){
  k.male <- k_clus_1
  k.female <- k_clus_2
} else {
  k.female <- k_clus_1
  k.male <- k_clus_2
}

cell.sex.ID <- list("male.cells" = k.male, "female.cells" = k.female)
saveRDS(object = cell.sex.ID, paste0(rds.path, 'sex_kmeans.RDS'))
#cell.sex.ID <- readRDS(file = paste0(rds.path, 'sex_kmeans.RDS'))

# add sex data to meta.data
norm.data@meta.data$sex <- unlist(lapply(rownames(norm.data@meta.data), function(x)
  if(x %in% k.male){"male"} else if(x %in% k.female){"female"} else{stop("cell sex is not assigned")}))


# Next - subset autosomal genes and Z genes - then calculate the sums/means/medians for each gene for both kmeans clustered cells
# and plot in order to compare whether the two groups show significant differences in their expression.

###### Z- chromosome genes ######

# this is a test on autosomal genes to try and calculate and compare FC betweeen clusters
# Calculating median is tricky as there are a lot of dropouts in 10x data so you end up with either 0s (when the median  = 0)
# or 1 (when the median expression in both clusters is the same - probably a result of normalisation resulting in a UMI of 0 or
# 1 being normalised to a nominal value)
mean.Z.male <- data.frame(Z.mean = apply(norm.data@assays$RNA[grepl("Z-", rownames(norm.data@assays$RNA)), k.male], 1, mean))
# add 1 before log2 as log2(1) = 0
mean.Z.male <- log2(mean.Z.male + 1)
mean.Z.female <- data.frame(Z.mean = apply(norm.data@assays$RNA[grepl("Z-", rownames(norm.data@assays$RNA)), k.female], 1, mean))
mean.Z.female <- log2(mean.Z.female + 1)

###### autosomal genes ######

# this is a test on autosomal genes to try and calculate and compare FC betweeen clusters
mean.auto.male <- data.frame(auto.mean = apply(norm.data@assays$RNA[!grepl("Z-", rownames(norm.data@assays$RNA)) & !grepl("W-", rownames(norm.data@assays$RNA)), k.male], 1, mean))
mean.auto.male <- log2(mean.auto.male + 1)
mean.auto.female <- data.frame(auto.mean = apply(norm.data@assays$RNA[!grepl("Z-", rownames(norm.data@assays$RNA)) & !grepl("W-", rownames(norm.data@assays$RNA)), k.female], 1, mean))
mean.auto.female <- log2(mean.auto.female + 1)

# calculate FC by subtracting log2 expression from each other
FC <- list()
FC$Z <- mean.Z.male - mean.Z.female
FC$auto <-  mean.auto.male - mean.auto.female

# plot hist/boxplot of data
hist(FC$Z$Z.mean)
hist(FC$auto$auto.mean)
pdf(paste0(curr.plot.path,"sex_kmeans_log2FC_boxplot.pdf"))
boxplot(c(FC$Z, FC$auto),  ylab = "male - female log2 FC (mean normalised UMI +1)", names = c("Z chromosome genes", "autosomal genes"))
dev.off()

# Z genes are upregulated within male genes relative to female genes whereas autosomal genes have a normal distribution of logFCs
# therefore Z genes should be filtered out - this could be done on a gene by gene basis, but can also be done through removing all Z genes

#####################################################################################################
#                                       Filter sex genes                                            #
#####################################################################################################

norm.data.sexfilt <- norm.data[rownames(norm.data)[!grepl("W-", rownames(norm.data)) & !grepl("Z-", rownames(norm.data))], ]

# re-run findvariablefeatures and scaling
norm.data.sexfilt <- FindVariableFeatures(norm.data.sexfilt, selection.method = "vst", nfeatures = 2000)
norm.data.sexfilt <- ScaleData(norm.data.sexfilt, features = rownames(norm.data.sexfilt), vars.to.regress = c("percent.mt", "sex"))

saveRDS(norm.data.sexfilt, paste0(rds.path, "norm.data.sexfilt.RDS"))


#####################################################################################################
#                    Perform dimensionality reduction by PCA and UMAP embedding                     #
#####################################################################################################

# read in RDS data if needed
# norm.data.sexfilt <- readRDS(paste0(rds.path, "norm.data.sexfilt.RDS"))

norm.data.sexfilt <- RunPCA(object = norm.data.sexfilt, verbose = FALSE)

pdf(paste0(curr.plot.path, "dimHM.pdf"),width=15,height=25)
DimHeatmap(norm.data.sexfilt, dims = 1:30, balanced = TRUE, cells = 500)
dev.off()

pdf(paste0(curr.plot.path, "elbowplot.pdf"),width=12,height=10)
print(ElbowPlot(norm.data.sexfilt, ndims = 40))
dev.off()

PCA.level.comparison(norm.data.sexfilt, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)

# Use PCA=15 as elbow plot is relatively stable across stages
norm.data.sexfilt <- FindNeighbors(norm.data.sexfilt, dims = 1:15, verbose = FALSE)
norm.data.sexfilt <- RunUMAP(norm.data.sexfilt, dims = 1:15, verbose = FALSE)
clust.res(seurat.obj = norm.data.sexfilt, multi.obj.list = F, by = 0.2)
# Use clustering resolution = 1.4 for subsequent filtering of poor quality clusters
# this increases the stringency of poor quality clusters, removing the least data possible
norm.data.sexfilt <- FindClusters(norm.data.sexfilt, resolution = 1.4, verbose = FALSE)

# plot UMAP for clusters and developmental stage
clust.stage.plot(norm.data.sexfilt)

QC_plot(norm.data.sexfilt)

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.data.sexfilt, only.pos = T, logfc.threshold = 0.25)
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
tenx_pheatmap(data = norm.data.sexfilt, metadata = c("orig.ident", "seurat_clusters"), used.genes = unique(top15$gene), width = 40, height = 30, main = "")


#####################################################################################################
#                 Identify and remove poor quality clusters / contamination (mesoderm)              #
#####################################################################################################

curr.plot.path <- past0(plot.path, "2_cluster_filt/")
dir.create(curr.plot.path)

# identify mesoderm
# UMAP plots GOI
genes <- c("EYA2", "SIX1", "TWIST1", "PITX2", "SOX17", "DAZL", "DND1", "CXCR4")
multi.feature.plot(seurat.obj = norm.data.sexfilt, multi.obj.list = F, gene.list = genes, cluster.res = 1.4, plot.clusters = T, 
                   plot.stage = T, n.col = 4, label = "", cluster.col = "RNA_snn_res.1.4")

# Dotplot for identifying PGCs, Early mesoderm and Late mesoderm
pdf(file = paste0(plot.path, "dotplot.GOI.pdf"), height = 8, width = 10)
DotPlot(norm.data.sexfilt, features = c( "SOX17", "CXCR4","EYA2", "TWIST1", "SIX1",  "PITX2", "DAZL"))
dev.off()

############################### Remove cells from clusters ########################################
# Clust 12, 14, 17 = poor quality cells
# Clust 18 = early mesoderm - expresses sox17, eya2, pitx2, cxcr4
# Clust 20 = Late mesoderm - expresses twist1, six1, eya2
# Clust 21 = PGC's - expresses dazl very highly
norm.data.clustfilt <- rownames(norm.data.sexfilt@meta.data)[norm.data.sexfilt@meta.data$seurat_clusters == 12 |
                                                               norm.data.sexfilt@meta.data$seurat_clusters == 14 | 
                                                               norm.data.sexfilt@meta.data$seurat_clusters == 17 |
                                                               norm.data.sexfilt@meta.data$seurat_clusters == 18 | 
                                                               norm.data.sexfilt@meta.data$seurat_clusters == 20 |
                                                               norm.data.sexfilt@meta.data$seurat_clusters == 21]

norm.data.clustfilt <- subset(norm.data.sexfilt, cells = norm.data.clustfilt, invert = T)

# re-run findvariablefeatures and scaling
norm.data.clustfilt <- FindVariableFeatures(norm.data.clustfilt, selection.method = "vst", nfeatures = 2000)
norm.data.clustfilt <- ScaleData(norm.data.clustfilt, features = rownames(norm.data.clustfilt), vars.to.regress = c("percent.mt", "sex"))

saveRDS(norm.data.clustfilt, paste0(rds.path, "norm.data.clustfilt.RDS"))

#####################################################################################################
#                    Perform dimensionality reduction by PCA and UMAP embedding                     #
#####################################################################################################

# read in RDS data if needed
# norm.data.clustfilt <- readRDS(paste0(rds.path, "norm.data.clustfilt.RDS"))

norm.data.clustfilt <- RunPCA(object = norm.data.clustfilt, verbose = FALSE)

pdf(paste0(curr.plot.path, "dimHM.pdf"),width=15,height=25)
DimHeatmap(norm.data.clustfilt, dims = 1:30, balanced = TRUE, cells = 500)
dev.off()
pdf(paste0(curr.plot.path, "elbowplot.pdf"),width=12,height=10)
print(ElbowPlot(norm.data.clustfilt, ndims = 40))
dev.off()

PCA.level.comparison(norm.data.clustfilt, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)

norm.data.clustfilt <- FindNeighbors(norm.data.clustfilt, dims = 1:15, verbose = FALSE)
norm.data.clustfilt <- RunUMAP(norm.data.clustfilt, dims = 1:15, verbose = FALSE)
clust.res(seurat.obj = norm.data.clustfilt, multi.obj.list = F, by = 0.2)
# Cluster data using resolution 0.8 (clusters are quite unstable)
norm.data.clustfilt <- FindClusters(norm.data.clustfilt, resolution = 0.8)

# plot UMAP for clusters and developmental stage
clust.stage.plot(norm.data.clustfilt)

####################################################################################
#                            Check for cell cycle effect                           #
####################################################################################
curr.plot.path <- paste0(plot.path, "3_cell_cycle/")
dir.create(curr.plot.path)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pre.cell.cycle.dat <- CellCycleScoring(norm.data.clustfilt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# re-run findvariablefeatures and scaling
norm.data.clustfilt.cc <- FindVariableFeatures(pre.cell.cycle.dat, selection.method = "vst", nfeatures = 2000)
norm.data.clustfilt.cc <- ScaleData(norm.data.clustfilt.cc, features = rownames(norm.data.clustfilt.cc), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

saveRDS(norm.data.clustfilt.cc, paste0(rds.path, "norm.data.clustfilt.cc.RDS"))

#####################################################################################################
#                    Perform dimensionality reduction by PCA and UMAP embedding                     #
#####################################################################################################

# read in RDS data if needed
# norm.data.clustfilt.cc <- readRDS(paste0(rds.path, "norm.data.clustfilt.cc.RDS"))

norm.data.clustfilt.cc <- RunPCA(object = norm.data.clustfilt.cc, verbose = FALSE)

pdf(paste0(curr.plot.path, "dimHM.pdf"),width=15,height=25)
DimHeatmap(norm.data.clustfilt.cc, dims = 1:30, balanced = TRUE, cells = 500)
dev.off()

pdf(paste0(curr.plot.path, "elbowplot.pdf"),width=12,height=10)
print(ElbowPlot(norm.data.clustfilt.cc, ndims = 40))
dev.off()

PCA.level.comparison(norm.data.clustfilt.cc, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)

norm.data.clustfilt.cc <- FindNeighbors(norm.data.clustfilt.cc, dims = 1:15, verbose = FALSE)
norm.data.clustfilt.cc <- RunUMAP(norm.data.clustfilt.cc, dims = 1:15, verbose = FALSE)
clust.res(seurat.obj = norm.data.clustfilt.cc, multi.obj.list = F, by = 0.2)
norm.data.clustfilt.cc <- FindClusters(norm.data.clustfilt.cc, resolution = 1.2)

# UMAP of cell cycle before and after regressing out
pdf(paste0(curr.plot.path, "cell.cycle.pdf"), width = 14, height = 7)
pre.plot <- DimPlot(pre.cell.cycle.dat, group.by = "Phase", reduction = "umap")
post.plot <- DimPlot(norm.data.clustfilt.cc, group.by = "Phase", reduction = "umap")
print(gridExtra::grid.arrange(pre.plot, post.plot, ncol = 2))
dev.off()

# plot UMAP for clusters and developmental stage
clust.stage.plot(norm.data.clustfilt.cc)

# plot cluster QC
QC_plot(norm.data.clustfilt.cc)

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.data.clustfilt.cc, only.pos = T, logfc.threshold = 0.25)
# re-order genes in top15 based on desired cluster order is subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = c(1,2,8,11,0,5,6,7,9,12,13,3,4,10,14)))
tenx_pheatmap(data = norm.data.clustfilt.cc, metadata = c("orig.ident", "seurat_clusters"),
              used.genes = unique(top15$gene), width = 40, height = 30, main = "", order.by = "seurat_clusters", custom_order = c(1,2,8,11,0,5,6,7,9,12,13,3,4,10,14))


#####################################################################################################
#                                        Cell type identification                                   #
#####################################################################################################

curr.plot.path <- paste0(plot.path, "4_cell_type_classification/")
dir.create(curr.plot.path)

###################### Nerual crest and placode cell type identification ###################### 
placNC_genes <- c(
  # delaminating NC
  "ETS1", "SOX10", "SOX8", "LMO4",  "TFAP2B", "SOX9",
  # NPB
  "TFAP2A", "DRAXIN", "MSX1", "CSRNP1", "PAX7", "BMP5", "MSX2",
  # NC
  "WNT6",
  # Placodes
  "PITX1", "PITX2", "ZNF385C",  "SIX1", "EYA2", "DLX6", "HOMER2"
)

# plot expression of NP genes
multi.feature.plot(norm.data.clustfilt.cc, stage.name = "all.stages", gene.list = placNC_genes, multi.obj.list = F,
                   plot.clusters = T, n.col = 5, basename = "plac.UMAP.GOI.pdf", label = "UMAP of genes used to identify non-NP clusters",
                   plot.stage = T)

# add neural crest and placodal cell types to metadata
norm.data.clustfilt.cc@meta.data$placNC_clust <- apply(norm.data.clustfilt.cc@meta.data, 1, function(x)
  if(x["seurat_clusters"] == 10){"Delaminating NC"} else if(x["seurat_clusters"] == 12){"NC Progenitors"}
  else if(x["seurat_clusters"] == 9){"Placodes"}else {NA})

# plot annotated neural crest and placodal clusters
png(paste0(curr.plot.path, "plac.clust.png"), width = 13, height = 10, res = 200, units = "cm")
DimPlot(norm.data.clustfilt.cc, group.by = "placNC_clust")  + ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

# plot dotplot for neural crest and placodal genes and clusters
png(paste0(curr.plot.path, "plac.dotplot.png"), width = 25, height = 10, res = 200, units = "cm")
DotPlot(norm.data.clustfilt.cc[, !is.na(norm.data.clustfilt.cc@meta.data$placNC_clust)], group.by = "placNC_clust", features = placNC_genes) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

###################### Neural tube cell identification ###################### 
NP_genes <- c(
  # HB
  "GBX2", "HOXA2", "HOXA3", "HOXB2", "KROX20",
  # MHB
  "WNT4", "PAX2", "FGF8", "WNT1",
  # Anterior
  "PAX6", "OTX2", "OLIG2" , "SIX3",
  # Ventral floor
  "SHH", "NKX2-2", "FOXA2"
)

# plot expression of NP genes
multi.feature.plot(norm.data.clustfilt.cc, stage.name = "all.stages", gene.list = NP_genes, multi.obj.list = F,
                   plot.clusters = T, n.col = 5, basename = "NP.UMAP.GOI.pdf", label = "UMAP of genes used to identify non-NP clusters",
                   plot.stage = T)

# add NP cell types to metadata
norm.data.clustfilt.cc@meta.data$NP_clust <- apply(norm.data.clustfilt.cc@meta.data, 1, function(x)
  if(x["seurat_clusters"] == 3){"Late Forebrain"} else if(x["seurat_clusters"] == 4){"Late Midbrain"}
  else if(x["seurat_clusters"] == 6){"Early Midbrain"} else if(x["seurat_clusters"] == 7){"Early Forebrain"}
  else if(x["seurat_clusters"] == 14){"Ventral Forebrain"} else if(x["seurat_clusters"] == 13){"Hindbrain"} 
  else {NA})


# plot annotated NP clusters
png(paste0(curr.plot.path, "NP.clust.png"), width = 13, height = 10, res = 200, units = "cm")
DimPlot(norm.data.clustfilt.cc, group.by = "NP_clust")  + ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

# plot dotplot for NP genes and clusters
norm.data.clustfilt.cc@meta.data$NP_clust <- factor(norm.data.clustfilt.cc@meta.data$NP_clust,
                                                    levels = c("Hindbrain", "Early Midbrain", "Late Midbrain", "Early Forebrain",
                                                               "Late Forebrain", "Ventral Forebrain"))

png(paste0(curr.plot.path, "NP.dotplot.png"), width = 25, height = 10, res = 200, units = "cm")
DotPlot(norm.data.clustfilt.cc[, !is.na(norm.data.clustfilt.cc@meta.data$NP_clust)], group.by = "NP_clust", features = NP_genes) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()






