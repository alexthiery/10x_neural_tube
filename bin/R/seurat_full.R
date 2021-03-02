#!/usr/bin/env Rscript

# Load packages
reticulate::use_python('/usr/bin/python3.7')
library(Seurat)

library(getopt)
library(future)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(Antler)
library(ggbeeswarm)
library(tidyverse)

# set arguments for Rscript
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer",
  'customFuncs', 'm', 2, "character",
  'networkGenes', 'd', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# set default location
if(length(commandArgs(trailingOnly = TRUE)) == 0){
  cat('No command line arguments provided, user defaults paths are set for running interactively in Rstudio on docker\n')
  opt$runtype = "user"
} else {
  if(is.null(opt$runtype)){
    stop("--runtype must be either 'user', 'nextflow' or 'docker'")
  }
  if(tolower(opt$runtype) != "docker" & tolower(opt$runtype) != "user" & tolower(opt$runtype) != "nextflow"){
    stop("--runtype must be either 'user', 'nextflow' or 'docker'")
  }
  if(tolower(opt$runtype) == "nextflow"){
    if(is.null(opt$customFuncs)){
      stop("--customFuncs path must be specified")
    }
    if(is.null(opt$networkGenes)){
      stop("--networkGenes path must be specified")
    }
  }
}

####################################################################
# user paths need to be defined here in order to run interactively #
####################################################################
if (opt$runtype == "user"){
  sapply(list.files('./bin/R/custom_functions/', full.names = T), source)
  
  plot.path = "./output/NF-downstream_analysis/plots/"
  rds.path = "./output/NF-downstream_analysis/RDS.files/"
  antler.dir = "./output/NF-downstream_analysis/antler.input/"
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  dir.create(antler.dir, recursive = T)
  
  ##################################
  # set path where data is located #
  ##################################
  data_path = "./output/NF-10x_alignment/cellrangerCounts"
  
  # read all files from dir
  files <- list.files(data_path, recursive = T, full.names = T)
  # remove file suffix
  file.path <- dirname(files)[!duplicated(dirname(files))]
  # make dataframe with tissue matching directory
  tissue = c("hh4", "hh6", "ss4", "ss8")
  matches <- sapply(tissue, function(x) file.path[grep(pattern = x, x = file.path)])
  sample.paths <- data.frame(tissue = names(matches), path = matches, row.names = NULL)
  
  # Read in network genes and remove 0 timepoint
  network_genes <- read.csv('./input_files/network_expression.csv') %>%
    filter(!timepoint == 0)

} else if (opt$runtype == "nextflow"){
  cat('pipeling running through nextflow\n')
  
  sapply(list.files(opt$customFuncs, full.names = T), source)
  
  plot.path = "plots/"
  dir.create(plot.path, recursive = T)
  rds.path = "RDS.files/"
  dir.create(rds.path, recursive = T)
  antler.dir = "antler.input/"
  dir.create(antler.dir, recursive = T)
  
  
  # read all files from folder and keep only those from chr_edit
  files <- list.files("./cellrangerCounts", recursive = T, full.names = T)
  # remove file suffix
  file.path <- dirname(files)[!duplicated(dirname(files))]
  # make dataframe with tissue matching directory
  tissue = c("hh4", "hh6", "ss4", "ss8")

  print(files)
  matches <- sapply(tissue, function(x) file.path[grep(pattern = x, x = file.path)])
  sample.paths <- data.frame(tissue = names(matches), path = matches, row.names = NULL)
  
  # Read in network genes and remove 0 timepoint
  network_genes <- read.csv(opt$networkGenes) %>%
    filter(!timepoint == 0)
}

# set number of cores to use for parallelisation
if(is.null(opt$cores)){ncores = 4}else{ncores= opt$cores}
cat(paste0("script ran with ", ncores, " cores\n"))

# Make Seurat objects for each of the different samples.
for(i in 1:nrow(sample.paths["path"])){
  name<-paste(sample.paths[i,"tissue"])
  assign(name, CreateSeuratObject(counts= Read10X(data.dir = paste(sample.paths[i,"path"])), project = paste(sample.paths[i, "tissue"])))
}

# The four Seurat objects are then merged, before running CreateSeuratObject again on the output in order to apply the min.cells parameter on the final merged dataset.
temp <- merge(hh4, y = c(hh6, ss4, ss8), add.cell.ids = c("hh4", "hh6", "ss4", "ss8"), project = "chick.10x")
merged.data<-CreateSeuratObject(GetAssayData(temp), min.cells = 3, project = "chick.10x.mincells3")


# make seurat object with ensembl names and save as separate dataframe for adding to misc slot
for(i in 1:nrow(sample.paths["path"])){
  name<-paste(sample.paths[i,"tissue"])
  assign(paste0(name, "_ensID"), CreateSeuratObject(counts= Read10X(data.dir = paste(sample.paths[i,"path"]), gene.column = 1), project = paste(sample.paths[i, "tissue"])))
}
temp <- merge(hh4_ensID, y = c(hh6_ensID, ss4_ensID, ss8_ensID), add.cell.ids = c("hh4", "hh6", "ss4", "ss8"), project = "chick.10x")
merged.data_ensID<-CreateSeuratObject(GetAssayData(temp), min.cells = 3, project = "chick.10x.mincells3")

# add gene IDs dataframe to merged data object
Misc(merged.data, slot = "geneIDs") <- cbind("gene_ID" = rownames(merged.data_ensID), "gene_name" =  rownames(merged.data))

# The original Seurat objects are then removed from the global environment
rm(hh4, hh6, ss4, ss8, sample.paths, temp, hh4_ensID, hh6_ensID, ss4_ensID, ss8_ensID, merged.data_ensID)

# Store mitochondrial percentage in object meta data
merged.data <- PercentageFeatureSet(merged.data, pattern = "^MT-", col.name = "percent.mt")


#####################################################################################################
#                           Filter data based on variable threshold                                 #
#####################################################################################################

# Remove data which do not pass filter threshold
merged.data <- subset(merged.data, subset = c(nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 15))

# Log normalize data and find variable features
norm.data <- NormalizeData(merged.data, normalization.method = "LogNormalize", scale.factor = 10000)
norm.data <- FindVariableFeatures(norm.data, selection.method = "vst", nfeatures = 2000)

# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

# Scale data and regress out MT content
norm.data <- ScaleData(norm.data, features = rownames(norm.data), vars.to.regress = "percent.mt")

# Save RDS after scaling as this step takes time
saveRDS(norm.data, paste0(rds.path, "norm.data.RDS"))

#####################################################################################################
#                    Perform dimensionality reduction by PCA and UMAP embedding                    #
#####################################################################################################

# Read in RDS data if needed
# norm.data <- readRDS(paste0(rds.path, "norm.data.RDS"))

# Change plot path
curr.plot.path <- paste0(plot.path, '0_filt_data/')
dir.create(curr.plot.path)

# Run PCA analysis on the each set of data
norm.data <- RunPCA(object = norm.data, verbose = FALSE)

#####################################################################################################
#   Seurat's clustering algorithm is based on principle components, so we need to ensure that only the informative PCs are kept!                   #
#####################################################################################################

# Plot heatmap of top variable genes across top principle components
png(paste0(curr.plot.path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(norm.data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

# another heuristic method is ElbowPlot which ranks PCs based on the % variance explained by each PC
png(paste0(curr.plot.path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(norm.data, ndims = 40))
graphics.off()

# Run clustering and UMAP at different PCA cutoffs - save this output to compare the optimal number of PCs to be used
png(paste0(curr.plot.path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(norm.data, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
# Use clustering resolution = 0.5 for filtering
norm.data <- FindNeighbors(norm.data, dims = 1:15, verbose = FALSE)
norm.data <- RunUMAP(norm.data, dims = 1:15, verbose = FALSE)
norm.data <- FindClusters(norm.data, resolution = 0.5, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(curr.plot.path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(norm.data)
graphics.off()

# Plot QC for each cluster
png(paste0(curr.plot.path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QC.plot(norm.data)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.data, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = norm.data, col.to.sort = seurat_clusters, sort.by = orig.ident)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(curr.plot.path, 'HM.top15.DE.png'), height = 50, width = 75, units = 'cm', res = 700)
tenx.pheatmap(data = norm.data, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster.order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters")
graphics.off()

#####################################################################################################
#     Heatmap clearly shows clusters segregate by sex - check this and regress out the sex effect   #
#####################################################################################################

# Change plot path
curr.plot.path <- paste0(plot.path, '1_sex_filt/')
dir.create(curr.plot.path)

# There is a strong sex effect - this plot shows DE genes between clusters 1 and 2 which are preodominantly hh4 clusters. Clustering is driven by sex genes
png(paste0(curr.plot.path, 'HM.top15.DE.pre-sexfilt.png'), height = 40, width = 70, units = 'cm', res = 500)
tenx.pheatmap(data = norm.data[,rownames(norm.data@meta.data[norm.data$seurat_clusters == 1 | norm.data$seurat_clusters == 2,])],
              metadata = c("seurat_clusters", "orig.ident"), selected_genes = rownames(FindMarkers(norm.data, ident.1 = 1, ident.2 = 2)),
              hclust_rows = T, gaps_col = "seurat_clusters")
graphics.off()

# Use W chromosome genes to K-means cluster the cells into male (zz) and female (zw)
W_genes <- as.matrix(norm.data@assays$RNA[grepl("W-", rownames(norm.data@assays$RNA)),])
k_clusters <- kmeans(t(W_genes), 2)
k_clusters <- data.frame(k_clusters$cluster)
norm.data@meta.data$k_clusters <- k_clusters[match(colnames(norm.data@assays$RNA), rownames(k_clusters)),]

# Get rownames for kmeans clusters 1 and 2
k_clus_1 <- rownames(norm.data@meta.data[norm.data@meta.data$k_clusters == 1,])
k_clus_2 <- rownames(norm.data@meta.data[norm.data@meta.data$k_clusters == 2,])

# K clustering identities are stochastic, so I need to identify which cluster is male and female
# Sum of W genes is order of magnitude greater in cluster 2 - these are the female cells
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
saveRDS(cell.sex.ID, paste0(rds.path, "sex_kmeans.RDS"))

# Add sex data to meta.data
norm.data@meta.data$sex <- unlist(lapply(rownames(norm.data@meta.data), function(x)
  if(x %in% k.male){"male"} else if(x %in% k.female){"female"} else{stop("cell sex is not assigned")}))

#####################################################################################################
#                             Remove W chromosome genes and regress sex effect                      #                      #
#####################################################################################################

# Init sexscale object
norm.data.sexscale <- norm.data[rownames(norm.data)[!grepl("W-", rownames(norm.data))],]

# Re-run findvariablefeatures and scaling
norm.data.sexscale <- FindVariableFeatures(norm.data.sexscale, selection.method = "vst", nfeatures = 2000)
# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

norm.data.sexscale <- ScaleData(norm.data.sexscale, features = rownames(norm.data.sexscale), vars.to.regress = c("percent.mt", "sex"))

# Save RDS
saveRDS(norm.data.sexscale, paste0(rds.path, "norm.data.sexscale.RDS"))

# Read in RDS data if needed
# norm.data.sexscale <- readRDS(paste0(rds.path, "norm.data.sexscale.RDS"))

# Set plot path
curr.plot.path <- paste0(plot.path, '1_sex_filt/')

# PCA
norm.data.sexscale <- RunPCA(object = norm.data.sexscale, verbose = FALSE)

png(paste0(curr.plot.path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(norm.data.sexscale, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(curr.plot.path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(norm.data.sexscale, ndims = 40))
graphics.off()

png(paste0(curr.plot.path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(norm.data.sexscale, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
norm.data.sexscale <- FindNeighbors(norm.data.sexscale, dims = 1:15, verbose = FALSE)
norm.data.sexscale <- RunUMAP(norm.data.sexscale, dims = 1:15, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(curr.plot.path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = norm.data.sexscale, by = 0.1)
graphics.off()

# Use clustering resolution = 0.5 to look for contamination clusters
norm.data.sexscale <- FindClusters(norm.data.sexscale, resolution = 0.5, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(curr.plot.path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(norm.data.sexscale)
graphics.off()

# Plot QC for each cluster
png(paste0(curr.plot.path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QC.plot(norm.data.sexscale)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.data.sexscale, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = norm.data.sexscale, col.to.sort = seurat_clusters, sort.by = orig.ident)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(curr.plot.path, 'HM.top15.DE.post-sexfilt.png'), height = 75, width = 100, units = 'cm', res = 500)
tenx.pheatmap(data = norm.data.sexscale, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster.order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters")
graphics.off()

#####################################################################################################
#                           Identify and remove contamination (mesoderm and PGCs)                   #
#####################################################################################################

# Change plot path
curr.plot.path <- paste0(plot.path, "2_contamination_filt/")
dir.create(curr.plot.path)

# Identify mesoderm and PGCs
# UMAP plots GOI
genes <- c("EYA2", "SIX1", "TWIST1", "PITX2", "SOX17", "DAZL", "CXCR4")

ncol = 3
png(paste0(curr.plot.path, "UMAP_GOI.png"), width = ncol*10, height = 10*ceiling((length(genes)+2)/ncol), units = "cm", res = 200)
multi.feature.plot(seurat.obj = norm.data.sexscale, gene.list = genes, plot.clusters = T,
                   plot.stage = T, label = "", cluster.col = "RNA_snn_res.0.5", n.col = ncol)
graphics.off()

# Dotplot for identifying PGCs, Early mesoderm and Late mesoderm
png(paste0(curr.plot.path, "dotplot.GOI.png"), width = 20, height = 12, units = "cm", res = 200)
DotPlot(norm.data.sexscale, features = c( "SOX17", "CXCR4","EYA2", "TWIST1", "SIX1",  "PITX2", "DAZL"))
graphics.off()

############################### Remove contaminating cells from clusters ########################################
# Clust 8 = early mesoderm - expresses sox17, eya2, pitx2, cxcr4
# Clust 10 = Late mesoderm - expresses twist1, six1, eya2
# Clust 11 = PGC's - expresses dazl very highly
norm.data.contamfilt <- rownames(norm.data.sexscale@meta.data)[norm.data.sexscale@meta.data$seurat_clusters ==  8 |
                                                                 norm.data.sexscale@meta.data$seurat_clusters == 10 |
                                                                 norm.data.sexscale@meta.data$seurat_clusters == 11]

norm.data.contamfilt <- subset(norm.data.sexscale, cells = norm.data.contamfilt, invert = T)

# Re-run findvariablefeatures and scaling
norm.data.contamfilt <- FindVariableFeatures(norm.data.contamfilt, selection.method = "vst", nfeatures = 2000)

# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

norm.data.contamfilt <- ScaleData(norm.data.contamfilt, features = rownames(norm.data.contamfilt), vars.to.regress = c("percent.mt", "sex"))

saveRDS(norm.data.contamfilt, paste0(rds.path, "norm.data.contamfilt.RDS"))

# Read in RDS data if needed
# norm.data.contamfilt <- readRDS(paste0(rds.path, "norm.data.contamfilt.RDS"))

# PCA
norm.data.contamfilt <- RunPCA(object = norm.data.contamfilt, verbose = FALSE)

png(paste0(curr.plot.path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(norm.data.contamfilt, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(curr.plot.path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(norm.data.contamfilt, ndims = 40))
graphics.off()

png(paste0(curr.plot.path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(norm.data.contamfilt, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
norm.data.contamfilt <- FindNeighbors(norm.data.contamfilt, dims = 1:15, verbose = FALSE)
norm.data.contamfilt <- RunUMAP(norm.data.contamfilt, dims = 1:15, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(curr.plot.path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = norm.data.contamfilt, by = 0.2)
graphics.off()

# Use clustering resolution = 1.4 in order to make lots of clusters and identify any remaining poor quality cells
norm.data.contamfilt <- FindClusters(norm.data.contamfilt, resolution = 1.4)

# Plot UMAP for clusters and developmental stage
png(paste0(curr.plot.path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(norm.data.contamfilt)
graphics.off()

# Plot QC for each cluster
png(paste0(curr.plot.path, "cluster.QC.png"), width=45, height=14, units = 'cm', res = 200)
QC.plot(norm.data.contamfilt)
graphics.off()



############################### Remove poor quality clusters ########################################

# Clust 11, 14, 16, 18 = poor quality cells
norm.data.clustfilt <- rownames(norm.data.contamfilt@meta.data)[norm.data.contamfilt@meta.data$seurat_clusters ==  11 |
                                                                  norm.data.contamfilt@meta.data$seurat_clusters == 14 |
                                                                  norm.data.contamfilt@meta.data$seurat_clusters == 16 |
                                                                  norm.data.contamfilt@meta.data$seurat_clusters == 18]

norm.data.clustfilt <- subset(norm.data.contamfilt, cells = norm.data.clustfilt, invert = T)

# Re-run findvariablefeatures and scaling
norm.data.clustfilt <- FindVariableFeatures(norm.data.clustfilt, selection.method = "vst", nfeatures = 2000)

# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

norm.data.clustfilt <- ScaleData(norm.data.clustfilt, features = rownames(norm.data.clustfilt), vars.to.regress = c("percent.mt", "sex"))

saveRDS(norm.data.clustfilt, paste0(rds.path, "norm.data.clustfilt.RDS"))

# Read in RDS data if needed
# norm.data.clustfilt <- readRDS(paste0(rds.path, "norm.data.clustfilt.RDS"))

# Change plot path
curr.plot.path <- paste0(plot.path, "3_cluster_filt/")
dir.create(curr.plot.path)

# PCA
norm.data.clustfilt <- RunPCA(object = norm.data.clustfilt, verbose = FALSE)

png(paste0(curr.plot.path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(norm.data.clustfilt, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(curr.plot.path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(norm.data.clustfilt, ndims = 40))
graphics.off()

png(paste0(curr.plot.path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(norm.data.clustfilt, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
norm.data.clustfilt <- FindNeighbors(norm.data.clustfilt, dims = 1:15, verbose = FALSE)
norm.data.clustfilt <- RunUMAP(norm.data.clustfilt, dims = 1:15, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(curr.plot.path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = norm.data.clustfilt, by = 0.2)
graphics.off()

# Use clustering resolution = 0.8
norm.data.clustfilt <- FindClusters(norm.data.clustfilt, resolution = 0.8)

# Plot UMAP for clusters and developmental stage
png(paste0(curr.plot.path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(norm.data.clustfilt)
graphics.off()

# Plot cluster QC
png(paste0(curr.plot.path, "cluster.QC.png"), width=45, height=14, units = 'cm', res = 200)
QC.plot(norm.data.clustfilt)
graphics.off()

####################################################################################
#                            Check for cell cycle effect                           #
####################################################################################

# Set plot path
curr.plot.path <- paste0(plot.path, "4_cell_cycle/")
dir.create(curr.plot.path)

# Calculate cell cycle for each cell
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pre.cell.cycle.dat <- CellCycleScoring(norm.data.clustfilt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Re-run findvariablefeatures and scaling
norm.data.clustfilt.cc <- FindVariableFeatures(pre.cell.cycle.dat, selection.method = "vst", nfeatures = 2000)
# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

norm.data.clustfilt.cc <- ScaleData(norm.data.clustfilt.cc, features = rownames(norm.data.clustfilt.cc), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

saveRDS(norm.data.clustfilt.cc, paste0(rds.path, "norm.data.clustfilt.cc.RDS"))

# Read in RDS data if needed
# norm.data.clustfilt.cc <- readRDS(paste0(rds.path, "norm.data.clustfilt.cc.RDS"))

# PCA
norm.data.clustfilt.cc <- RunPCA(object = norm.data.clustfilt.cc, verbose = FALSE)

png(paste0(curr.plot.path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(norm.data.clustfilt.cc, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(curr.plot.path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(norm.data.clustfilt.cc, ndims = 40))
graphics.off()

png(paste0(curr.plot.path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(norm.data.clustfilt.cc, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
norm.data.clustfilt.cc <- FindNeighbors(norm.data.clustfilt.cc, dims = 1:15, verbose = FALSE)
norm.data.clustfilt.cc <- RunUMAP(norm.data.clustfilt.cc, dims = 1:15, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(curr.plot.path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = norm.data.clustfilt.cc, by = 0.2)
graphics.off()

# Use clustering resolution = 1.2
norm.data.clustfilt.cc <- FindClusters(norm.data.clustfilt.cc, resolution = 1.2)

# UMAP of cell cycle before and after regressing out
png(paste0(curr.plot.path, "cell.cycle.png"), width=40, height=20, units = 'cm', res = 200)
pre.plot <- DimPlot(pre.cell.cycle.dat, group.by = "Phase", reduction = "umap")
post.plot <- DimPlot(norm.data.clustfilt.cc, group.by = "Phase", reduction = "umap")
print(gridExtra::grid.arrange(pre.plot, post.plot, ncol = 2))
graphics.off()

# Plot UMAP for clusters and developmental stage
png(paste0(curr.plot.path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(norm.data.clustfilt.cc)
graphics.off()

# Plot QC for each cluster
png(paste0(curr.plot.path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QC.plot(norm.data.clustfilt.cc)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.data.clustfilt.cc, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = norm.data.clustfilt.cc, col.to.sort = seurat_clusters, sort.by = orig.ident)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(curr.plot.path, 'HM.top15.DE.png'), height = 75, width = 100, units = 'cm', res = 500)
tenx.pheatmap(data = norm.data.clustfilt.cc, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster.order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters")
graphics.off()


########################################################################################################
#                                      Cell state classification                                    #
########################################################################################################

# Set plot path
curr.plot.path <- paste0(plot.path, "5_cell_state_classification/")
dir.create(curr.plot.path)


# Genes of interest used for cell state classification
GOI = rev(c( "EOMES", "ADMP","YEATS4", "MAFA", "ING5", "LIN28B", "AATF","SETD2",
             "GATA2", "DLX5", "TFAP2A", "BMP4", "SIX1", "EYA2",
             "MSX1", "PAX7", "CSRNP1", "SOX10",
             "SOX2", "SOX21", "SIX3", "OLIG2", "PAX6", "FOXA2", "SHH", "PAX2", "WNT4", "HOXB2", "HOXA2", "GBX2"))

# Change order or clusters for plotting dotplots
cluster_order <- factor(norm.data.clustfilt.cc$seurat_clusters, levels = c(12,5,2,1,
                                                                           3,8,0,
                                                                           11, 10,
                                                                           7,4,14,
                                                                           9,6,13))

# Set factor levels for plotting
norm.data.clustfilt.cc$seurat_clusters <- cluster_order

# Generate pie charts for cluster dev stage composition
venn_data <- norm.data.clustfilt.cc@meta.data %>%
  rownames_to_column('cell_name') %>%
  dplyr::select(c(cell_name, orig.ident, seurat_clusters)) %>%
  group_by(seurat_clusters) %>%
  count(orig.ident, .drop = FALSE)

venn_data <- venn_data %>%
  mutate(total_cells = sum(n)) %>%
  mutate(n = n/total_cells)


# Reverse levels to deal with seurat dotplot reversing y axis
norm.data.clustfilt.cc$seurat_clusters <- fct_rev(cluster_order)

# Generate DotPlot
dotplot <- DotPlot(norm.data.clustfilt.cc, group.by = "seurat_clusters", features = GOI) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position="bottom", legend.box = "horizontal",
        axis.title.x=element_blank(), legend.text=element_text(size=10),
        legend.title=element_text(size=12))

# Generate Pie charts
pies <- ggplot(venn_data, aes(x=as.numeric(total_cells)/2, y=as.numeric(n), fill=orig.ident, width = total_cells)) +
  geom_bar(position = 'fill', stat = "identity") +
  facet_wrap( ~ seurat_clusters, nrow = nlevels(norm.data.clustfilt.cc@meta.data$seurat_clusters)) +
  coord_polar("y") +
  theme_void() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        legend.position=c(0,0), legend.direction = "horizontal",
        plot.margin = margin(t = 14, r = 0, b = 119, l = -110, unit = "pt"),
        legend.margin=margin(t = 70, r = 0, b = -100, l = -200, unit = "pt"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12)) +
  labs(fill = "Stage")

# Plot dotplot and pies
png(paste0(curr.plot.path, "dotplot.png"), width = 32, height = 18, units = "cm", res = 200)
print(plot_grid(dotplot, pies, rel_widths = c(5,1)))
graphics.off()


# Plot dotplot with cell classification labels
dotplot <- DotPlot(norm.data.clustfilt.cc, group.by = "seurat_clusters", features = GOI) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position="bottom", legend.box = "horizontal",
        axis.title.x=element_blank(), legend.text=element_text(size=10),
        axis.title.y = element_blank(), 
        legend.title=element_text(size=12)) +
  scale_y_discrete(labels = rev(c('Node', 'HH4-1', 'HH4-2', 'HH4-3', 'Non-Neural Plate Progenitors',
                                  'Placodes', 'Neural Plate Progenitors', 'NC Progenitors', 'Delaminating NC',
                                  'Early Forebrain', 'Late Forebrain', 'Ventral Forebrain', 'Early Midbrain', 'Late Midbrain', 'Hindbrain')))

png(paste0(curr.plot.path, "dotplot_classified.png"), width = 36, height = 18, units = "cm", res = 200)
print(plot_grid(dotplot, pies, rel_widths = c(5,1)))
graphics.off()


########################################################################################################
#                            Load Antler data and generate gene modules                                #
########################################################################################################


# Set plot path
curr.plot.path <- paste0(plot.path, "6_gene_modules/")
dir.create(curr.plot.path, recursive = TRUE)

# Extract expression data and make dataset compatible with antler

# strip end of cell names as this is incorrectly reformated in Antler
norm.data.clustfilt.cc <- RenameCells(norm.data.clustfilt.cc, new.names = sub('-.*','',colnames(norm.data.clustfilt.cc)))

antler_data <- data.frame(row.names = colnames(norm.data.clustfilt.cc),
                          "timepoint" = as.numeric(substring(colnames(norm.data.clustfilt.cc), 3, 3)),
                          "treatment" = rep("null", ncol(norm.data.clustfilt.cc)),
                          "replicate_id" = rep(1, ncol(norm.data.clustfilt.cc))
)
# save pheno data
write.table(antler_data, file = paste0(antler.dir, "phenoData.csv"), row.names = T, sep = "\t", col.names = T)
# save count data
write.table(GetAssayData(norm.data.clustfilt.cc, assay = "RNA", slot = "counts"), file = paste0(antler.dir, "assayData.csv"), row.names = T, sep = "\t", col.names = T, quote = F)

# Load data into antler
antler <- Antler$new(output_folder = curr.plot.path, num_cores = 4)
antler$load_dataset(folder_path = antler.dir)

# Remove genes which do not have >= 1 UMI count in >= 10 cells
antler$exclude_unexpressed_genes(min_cells=10, min_level=1, verbose=T, data_status='Raw')

# Normalise data
antler$normalize(method = 'MR')

# Calculate unbiased gene modules
antler$gene_modules$identify(
  name                  = "unbiasedGMs",
  corr_t                = 0.3,  # the Spearman correlation treshold
  corr_min              = 3,    # min. number of genes a gene must correlate with
  mod_min_cell          = 10,   # min. number of cells expressing the module
  mod_consistency_thres = 0.4,  # ratio of expressed genes among "positive" cells
  process_plots         = TRUE)

# Get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = norm.data.clustfilt.cc, col.to.sort = seurat_clusters, sort.by = orig.ident)

# plot all gene modules
png(paste0(curr.plot.path, 'allmodules.png'), height = 150, width = 120, units = 'cm', res = 500)
GM.plot(data = norm.data.clustfilt.cc, metadata = c("seurat_clusters", "orig.ident"), gene_modules = antler$gene_modules$lists$unbiasedGMs$content,
        show_rownames = F, custom_order = cluster.order, custom_order_column = "seurat_clusters")
graphics.off()

# Plot gene modules with at least 50% of genes DE > 0.25 logFC & FDR < 0.001
# Find DEGs
DEgenes <- FindAllMarkers(norm.data.clustfilt.cc, only.pos = T, logfc.threshold = 0.25) %>% filter(p_val_adj < 0.001)

# Filter GMs with 50% genes DE logFC > 0.25 & FDR < 0.001
gms <- subset.gm(antler$gene_modules$lists$unbiasedGMs$content, selected_genes = DEgenes$gene, keep_mod_ID = T, selected_gene_ratio = 0.5)

png(paste0(curr.plot.path, 'DE.GM.png'), height = 160, width = 80, units = 'cm', res = 500)
GM.plot(data = norm.data.clustfilt.cc, metadata = c("seurat_clusters", "orig.ident"), gene_modules = gms, gaps_col = "seurat_clusters",
        show_rownames = T, custom_order = cluster.order, custom_order_column = "seurat_clusters")
graphics.off()


# Screen DE GMs for neural induction GRN genes
# Intersect DE GMs with network genes
filtered_gms <- lapply(gms, function(x) x[x %in% network_genes$gene])

# Remove empty list elements
filtered_gms <- filtered_gms[lapply(filtered_gms,length)>0]

png(paste0(curr.plot.path, 'network.GM.png'), height = 40, width = 80, units = 'cm', res = 500)
GM.plot(data = norm.data.clustfilt.cc, metadata = c("seurat_clusters", "orig.ident"), gene_modules = filtered_gms, gaps_col = "seurat_clusters",
        show_rownames = T, custom_order = cluster.order, custom_order_column = "seurat_clusters")
graphics.off()


#######################################################################################################
#                                     Subset neural clusters                                          #
#######################################################################################################

# Set plot path
curr.plot.path <- paste0(plot.path, "7_neural_subset/")
dir.create(curr.plot.path)

# Subset neural clusters
# Clust 12, 3, 11, 8, 10 = not neural clusters of interest
neural_subset <- rownames(norm.data.clustfilt.cc@meta.data)[norm.data.clustfilt.cc@meta.data$seurat_clusters ==  12 |
                                                                    norm.data.clustfilt.cc@meta.data$seurat_clusters == 3 |
                                                                    norm.data.clustfilt.cc@meta.data$seurat_clusters == 11 |
                                                                    norm.data.clustfilt.cc@meta.data$seurat_clusters == 8|
                                                                    norm.data.clustfilt.cc@meta.data$seurat_clusters == 10]

neural_subset <- subset(norm.data.clustfilt.cc, cells = neural_subset, invert = T)

# Re-run findvariablefeatures and scaling
neural_subset <- FindVariableFeatures(neural_subset, selection.method = "vst", nfeatures = 2000)

# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

neural_subset <- ScaleData(neural_subset, features = rownames(neural_subset), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

# PCA
neural_subset <- RunPCA(object = neural_subset, verbose = FALSE)

saveRDS(neural_subset, paste0(rds.path, "neural_subset.RDS"))

# Read in RDS data if needed
# neural_subset <- readRDS(paste0(rds.path, "neural_subset.RDS"))

png(paste0(curr.plot.path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(neural_subset, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(curr.plot.path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(neural_subset, ndims = 40))
graphics.off()

png(paste0(curr.plot.path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(neural_subset, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
neural_subset <- FindNeighbors(neural_subset, dims = 1:15, verbose = FALSE)
neural_subset <- RunUMAP(neural_subset, dims = 1:15, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(curr.plot.path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = neural_subset, by = 0.2)
graphics.off()

# Use clustering resolution = 1.2
neural_subset <- FindClusters(neural_subset, resolution = 1.2)

# Plot UMAP for clusters and developmental stage
png(paste0(curr.plot.path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(neural_subset)
graphics.off()

# Plot QC for each cluster
png(paste0(curr.plot.path, "cluster.QC.png"), width=40, height=14, units = 'cm', res = 200)
QC.plot(neural_subset)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(neural_subset, only.pos = T, logfc.threshold = 0.25)

# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = neural_subset, col.to.sort = seurat_clusters, sort.by = orig.ident)

# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(curr.plot.path, 'HM.top15.DE.png'), height = 75, width = 100, units = 'cm', res = 500)
tenx.pheatmap(data = neural_subset, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster.order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters")
graphics.off()



######################################################################################################
#                                           Pseudotime                                               #
######################################################################################################

# Set plot path
curr.plot.path <- paste0(plot.path, "8_pseudotime/")
dir.create(curr.plot.path)


# Extract PC1 values
pc1 <- neural_subset@meta.data[,'orig.ident', drop=F] %>%
  tibble::rownames_to_column(var = "cell_name") %>%
  mutate(pc1 = Embeddings(object = neural_subset[["pca"]])[, 1])


# Plot cell stage along PC1
png(paste0(curr.plot.path, 'pc1.png'), height = 18, width = 26, units = 'cm', res = 400)
ggplot(pc1, aes(x = pc1, y = orig.ident, colour = orig.ident)) +
  geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")
graphics.off()

# Developmental time negatively correlated with pc1 - inverse PC1 and calculate cell rank.
# Rank is then converted to 0-99 pseudotime scale
pc1_rank <- pc1 %>%
  mutate(rank = rank(-pc1)) %>%
  mutate(pseudotime = rank*(99/max(rank)))


# Filter scaled seurat counts by network_genes -> generate long df with corresponding pseudotime and normalised count
plot_data <- as.data.frame(t(as.matrix(GetAssayData(neural_subset, slot = 'scale.data')[rownames(neural_subset) %in% network_genes$gene_name,]))) %>%
  tibble::rownames_to_column(var = "cell_name") %>%
  dplyr::full_join(pc1_rank) %>%
  pivot_longer(!c(cell_name, orig.ident, pseudotime, pc1, rank), names_to = "gene_name", values_to = "scaled_expression") %>%
  dplyr::left_join(network_genes) %>%
  droplevels() %>%
  mutate(timepoint = factor(timepoint, levels = c(1,3,5,7,9,12)))


# mean and SE summary data
plot_data_summary <- plot_data %>%
  mutate(rank_bin = pseudotime - (pseudotime %% 2.5)) %>% 
  group_by(rank_bin, timepoint) %>% 
  summarise(mn = mean(scaled_expression),
            se = sd(scaled_expression)/sqrt(n()))

# plot gam for all stages without standard error
png(paste0(curr.plot.path, 'gam_pseudotime_allnetwork.png'), height = 18, width = 26, units = 'cm', res = 400)
ggplot(plot_data, aes(x = pseudotime, y = scaled_expression, colour = timepoint)) +
  scale_color_manual(values=c("#6600ff", "#0000ff", "#009900", "#ffcc00", "#ff8533", "#ee0000")) +
  geom_smooth(method="gam", se=FALSE) +
  theme_classic()
graphics.off()

# plot gam for all stages with standard error
png(paste0(curr.plot.path, 'gam_se_pseudotime_allnetwork.png'), height = 18, width = 26, units = 'cm', res = 400)
ggplot(plot_data, aes(x = pseudotime, y = scaled_expression, colour = timepoint)) +
  geom_errorbar(data = plot_data_summary, aes(x = rank_bin, y = mn, 
                                              ymax = mn + se, ymin = mn - se), width = 2) +
  geom_point(data = plot_data_summary, aes(x = rank_bin, y = mn)) +
  scale_color_manual(values=c("#6600ff", "#0000ff", "#009900", "#ffcc00", "#ff8533", "#ee0000")) +
  geom_smooth(method="gam", se=FALSE) +
  theme_classic()
graphics.off()


#