#!/usr/bin/env Rscript

# In order to be able to run the script from either Rstudio, local terminal, or cluster terminal, I add a switch which looks for command line arguments. This then sets the directory paths accordingly.
library('getopt')

# set arguments for Rscript
spec = matrix(c(
  'location', 'l', 2, "character",
  'cores'   , 'c', 2, "integer",
  'samples' , 's', 2, "character",
  'myfuncs', 'm', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# set default location
if(is.null(opt$location)){opt$location = "docker"}

if (opt$location == "docker"){
  cat('Script running in Rstudio on docker\n')
  
  sapply(list.files('./R/my_functions/', full.names = T), source)
  
  plot.path = "./results/plots/1_seurat_full/"
  rds.path = "./results/RDS.files/1_seurat_full/"
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  
  # read all files from folder and keep only those from chr_edit
  files <- Filter(function(x) grepl("chr_edit", x), list.files("./data/cellranger_output", recursive = T, full.names = T))
  
  # remove file suffix
  file.path <- dirname(files)[!duplicated(dirname(files))]
  
  # make dataframe with tissue matching directory
  tissue = c("hh4", "hh6", "4ss", "8ss")
  matches <- sapply(tissue, function(x) file.path[grep(pattern = x, x = file.path)])
  sample.paths <- data.frame(tissue = names(matches), path = matches, row.names = NULL)
  sample.paths$tissue = c("hh4", "hh6", "ss4", "ss8")
  
} else if (opt$location == "CAMP"){
  cat('data loaded from CAMP\n')
  
  sapply(list.files(opt$myfuncs, full.names = T), source)
  
  plot.path = "plots/"
  dir.create(plot.path, recursive = T)
  rds.path = "rds.files/"
  dir.create(rds.path, recursive = T)
  
  
  # read all files from folder and keep only those from chr_edit
  files <- Filter(function(x) grepl("chr_edit", x), list.files(opt$samples, recursive = T, full.names = T))
  
  # remove file suffix
  file.path <- dirname(files)[!duplicated(dirname(files))]
  
  # make dataframe with tissue matching directory
  tissue = c("hh4", "hh6", "4ss", "8ss")
  matches <- sapply(tissue, function(x) file.path[grep(pattern = x, x = file.path)])
  sample.paths <- data.frame(tissue = names(matches), path = matches, row.names = NULL)
  sample.paths$tissue = c("hh4", "hh6", "ss4", "ss8")
  
  # set number of cores to use for parallelisation
  if(is.null(opt$cores)){ncores = 4}else{ncores= opt$cores}
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  
} else {stop("Script can only be ran locally or on CAMP")}

# Load packages - packages are stored within renv in the repository
reticulate::use_python('/usr/bin/python3.7')
library(Seurat)

library(future)
library(dplyr)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)

# Make Seurat objects for each of the different samples.
for(i in 1:nrow(sample.paths["path"])){
  name<-paste(sample.paths[i,"tissue"])
  assign(name, CreateSeuratObject(counts= Read10X(data.dir = paste(sample.paths[i,"path"])), project = paste(sample.paths[i, "tissue"])))
}

# The four Seurat objects are then merged, before running CreateSeuratObject again on the output in order to apply the min.cells parameter on the final merged dataset.
temp <- merge(hh4, y = c(hh6, ss4, ss8), add.cell.ids = c("hh4", "hh6", "ss4", "ss8"), project = "chick.10x")
merged.data<-CreateSeuratObject(GetAssayData(temp), min.cells = 3, project = "chick.10x.mincells3")

# The original Seurat objects are then removed from the global environment
rm(hh4, hh6, ss4, ss8, sample.paths, temp)

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

# Scale data and regress out MT content
# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

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
png(paste0(curr.plot.path, "cluster.QC.png"), width=36, height=14, units = 'cm', res = 200)
QC.plot(norm.data)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.data, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = norm.data, col.to.sort = seurat_clusters, sort.by = orig.ident)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(curr.plot.path, 'HM.top15.DE.png'), height = 50, width = 75, units = 'cm', res = 200)
tenx.pheatmap(data = norm.data, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster.order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters")
graphics.off()

#####################################################################################################
#           Heatmap clearly shows clusters segregate by sex - check this and remove sex genes       #
#####################################################################################################

# Change plot path
curr.plot.path <- paste0(plot.path, '1_sex_filt/')
dir.create(curr.plot.path)

# There is a strong sex effect - this plot shows DE genes between clusters 1 and 2 which are preodominantly hh4 clusters. Clustering is driven by sex genes
png(paste0(curr.plot.path, 'HM.top15.DE.pre-sexfilt.png'), height = 40, width = 70, units = 'cm', res = 200)
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

# K clustering identities are stochastic, so I mneed to identify which cluster is male and female
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


###### Next - subset autosomal genes and Z genes - then calculate the average for each gene for both kmeans clustered cells and plot in order to compare whether the two groups show significant differences in their expression.

# This is a test on autosomal genes to try and calculate and compare FC betweeen clusters
# Calculating median is tricky as there are a lot of dropouts in 10x data so you end up with either 0s (when the median  = 0) or 1 (when the median expression in both clusters is the same - probably a result of normalisation resulting in a UMI of 0 or 1 being normalised to a nominal value)

# Make dataframe for mean Z expression in male cells
mean.Z.male <- data.frame(Z.mean = apply(norm.data@assays$RNA[grepl("Z-", rownames(norm.data@assays$RNA)), k.male], 1, mean))
# add 1 before log2 as log2(1) = 0
mean.Z.male <- log2(mean.Z.male + 1)

# Make dataframe for mean Z expression in female cells
mean.Z.female <- data.frame(Z.mean = apply(norm.data@assays$RNA[grepl("Z-", rownames(norm.data@assays$RNA)), k.female], 1, mean))
mean.Z.female <- log2(mean.Z.female + 1)

# Make dataframe for mean autosomal expression in male cells
mean.auto.male <- data.frame(auto.mean = apply(norm.data@assays$RNA[!grepl("Z-", rownames(norm.data@assays$RNA)) & !grepl("W-", rownames(norm.data@assays$RNA)), k.male], 1, mean))
mean.auto.male <- log2(mean.auto.male + 1)

# Make dataframe for mean autosomal expression in male cells
mean.auto.female <- data.frame(auto.mean = apply(norm.data@assays$RNA[!grepl("Z-", rownames(norm.data@assays$RNA)) & !grepl("W-", rownames(norm.data@assays$RNA)), k.female], 1, mean))
mean.auto.female <- log2(mean.auto.female + 1)

# Calculate FC by subtracting log2 expression from each other
FC <- list()
FC$Z <- mean.Z.male - mean.Z.female
FC$auto <-  mean.auto.male - mean.auto.female

# Plot boxplot of Z gene and autosomal expression in male vs female cells
png(paste0(curr.plot.path,"sex_kmeans_log2FC_boxplot.png"), height = 18, width = 18, units = "cm", res = 200)
boxplot(c(FC$Z, FC$auto),  ylab = "male - female log2 FC (mean normalised UMI +1)", names = c("Z chromosome genes", "autosomal genes"))
graphics.off()

# Z genes are upregulated within male genes relative to female genes whereas autosomal genes have a normal distribution of logFCs
# therefore Z genes should be filtered out - this could be done on a gene by gene basis, but can also be done through removing all Z genes

#####################################################################################################
#                                       Filter sex genes                                            #
#####################################################################################################

# Following subsetting of cells and/or genes, the same pipeline as above is repeated i.e.
# Find variable features -> Scale data/regress out confounding variables -> PCA -> Find neighbours -> Run UMAP -> Find Clusters -> Cluster QC -> Find top DE genes


# Remove sex genes
norm.data.nosexfilt <- norm.data

# Re-run findvariablefeatures and scaling
norm.data.nosexfilt <- FindVariableFeatures(norm.data.nosexfilt, selection.method = "vst", nfeatures = 2000)
# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

norm.data.nosexfilt <- ScaleData(norm.data.nosexfilt, features = rownames(norm.data.nosexfilt), vars.to.regress = c("percent.mt", "sex"))

# Save RDS
saveRDS(norm.data.nosexfilt, paste0(rds.path, "norm.data.nosexfilt.RDS"))




# Remove sex genes
norm.data.malefilt <- norm.data[rownames(norm.data)[!grepl("W-", rownames(norm.data))], ]

# Re-run findvariablefeatures and scaling
norm.data.malefilt <- FindVariableFeatures(norm.data.malefilt, selection.method = "vst", nfeatures = 2000)
# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

norm.data.malefilt <- ScaleData(norm.data.malefilt, features = rownames(norm.data.malefilt), vars.to.regress = c("percent.mt", "sex"))

# Save RDS
saveRDS(norm.data.malefilt, paste0(rds.path, "norm.data.malefilt.RDS"))






# Remove sex genes
norm.data.sexfilt <- norm.data[rownames(norm.data)[!grepl("W-", rownames(norm.data)) & !grepl("Z-", rownames(norm.data))], ]

# Re-run findvariablefeatures and scaling
norm.data.sexfilt <- FindVariableFeatures(norm.data.sexfilt, selection.method = "vst", nfeatures = 2000)
# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

norm.data.sexfilt <- ScaleData(norm.data.sexfilt, features = rownames(norm.data.sexfilt), vars.to.regress = c("percent.mt", "sex"))

# Save RDS
saveRDS(norm.data.sexfilt, paste0(rds.path, "norm.data.sexfilt.RDS"))

# Read in RDS data if needed
# norm.data.sexfilt <- readRDS(paste0(rds.path, "norm.data.sexfilt.RDS"))

# Set plot path
curr.plot.path <- paste0(plot.path, '1_sex_filt/')

# PCA
norm.data.sexfilt <- RunPCA(object = norm.data.sexfilt, verbose = FALSE)

png(paste0(curr.plot.path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(norm.data.sexfilt, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(curr.plot.path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(norm.data.sexfilt, ndims = 40))
graphics.off()

png(paste0(curr.plot.path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(norm.data.sexfilt, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
norm.data.sexfilt <- FindNeighbors(norm.data.sexfilt, dims = 1:15, verbose = FALSE)
norm.data.sexfilt <- RunUMAP(norm.data.sexfilt, dims = 1:15, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(curr.plot.path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = norm.data, by = 0.2)
graphics.off()

# Use clustering resolution = 1.4 for subsequent filtering of poor quality clusters this increases the stringency of poor quality clusters, removing the least data possible
norm.data.sexfilt <- FindClusters(norm.data.sexfilt, resolution = 1.4, verbose = FALSE)

# Plot UMAP for clusters and developmental stage
png(paste0(curr.plot.path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(norm.data.sexfilt)
graphics.off()

# Plot QC for each cluster
png(paste0(curr.plot.path, "cluster.QC.png"), width=36, height=14, units = 'cm', res = 200)
QC.plot(norm.data.sexfilt)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.data.sexfilt, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = norm.data.sexfilt, col.to.sort = seurat_clusters, sort.by = orig.ident)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(curr.plot.path, 'HM.top15.DE.post-sexfilt.png'), height = 75, width = 100, units = 'cm', res = 200)
tenx.pheatmap(data = norm.data.sexfilt, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster.order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters")
graphics.off()

#####################################################################################################
#                 Identify and remove poor quality clusters / contamination (mesoderm and PGCs)     #
#####################################################################################################

# Change plot path
curr.plot.path <- paste0(plot.path, "2_cluster_filt/")
dir.create(curr.plot.path)

# Identify mesoderm and PGCs
# UMAP plots GOI
genes <- c("EYA2", "SIX1", "TWIST1", "PITX2", "SOX17", "DAZL", "DND1", "CXCR4")

ncol = 4
png(paste0(curr.plot.path, "UMAP_GOI.png"), width = ncol*10, height = 10*ceiling((length(genes)+1)/ncol), units = "cm", res = 200)
multi.feature.plot(seurat.obj = norm.data.sexfilt, gene.list = genes, plot.clusters = T,
                   plot.stage = T, label = "", cluster.col = "RNA_snn_res.1.4", n.col = ncol)
graphics.off()

# Dotplot for identifying PGCs, Early mesoderm and Late mesoderm
png(paste0(curr.plot.path, "dotplot.GOI.png"), width = 20, height = 12, units = "cm", res = 200)
DotPlot(norm.data.sexfilt, features = c( "SOX17", "CXCR4","EYA2", "TWIST1", "SIX1",  "PITX2", "DAZL"))
graphics.off()

############################### Remove contaminating cells from clusters ########################################
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

# Re-run findvariablefeatures and scaling
norm.data.clustfilt <- FindVariableFeatures(norm.data.clustfilt, selection.method = "vst", nfeatures = 2000)

# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

norm.data.clustfilt <- ScaleData(norm.data.clustfilt, features = rownames(norm.data.clustfilt), vars.to.regress = c("percent.mt", "sex"))

saveRDS(norm.data.clustfilt, paste0(rds.path, "norm.data.clustfilt.RDS"))

# Read in RDS data if needed
# norm.data.clustfilt <- readRDS(paste0(rds.path, "norm.data.clustfilt.RDS"))

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

####################################################################################
#                            Check for cell cycle effect                           #
####################################################################################

# Set plot path
curr.plot.path <- paste0(plot.path, "3_cell_cycle/")
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
png(paste0(curr.plot.path, "cluster.QC.png"), width=36, height=14, units = 'cm', res = 200)
QC.plot(norm.data.clustfilt.cc)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.data.clustfilt.cc, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = norm.data.clustfilt.cc, col.to.sort = seurat_clusters, sort.by = orig.ident)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(curr.plot.path, 'HM.top15.DE.png'), height = 75, width = 100, units = 'cm', res = 200)
tenx.pheatmap(data = norm.data.clustfilt.cc, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster.order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters")
graphics.off()

#####################################################################################################
#                                        Cell type identification                                   #
#####################################################################################################

# Set plot path
curr.plot.path <- paste0(plot.path, "4_stage_split/")
dir.create(curr.plot.path)

# Split dataset into different stages
seurat_stage <- lapply(c('hh4', 'hh6', 'ss4', 'ss8'),
                       function(x) subset(norm.data.clustfilt.cc, cells = rownames(norm.data.clustfilt.cc@meta.data)[norm.data.clustfilt.cc$orig.ident == x]))
names(seurat_stage) = c('hh4', 'hh6', 'ss4', 'ss8')


# Re-run findvariablefeatures and scaling
seurat_stage <- lapply(seurat_stage, function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000))
# Enable parallelisation on camp
# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

seurat_stage <- lapply(seurat_stage, function(x) ScaleData(x, features = rownames(norm.data.clustfilt.cc),
                                                           vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score")))

saveRDS(seurat_stage, paste0(rds.path, "seurat_stage.RDS"))

# Read in RDS data if needed
# seurat_stage <- readRDS(paste0(rds.path, "seurat_stage.RDS"))

# PCA
seurat_stage <- lapply(seurat_stage, function(x) RunPCA(object = x, verbose = FALSE))

for(stage in names(seurat_stage)){
  png(paste0(curr.plot.path, "dimHM.",stage,".png"), width=30, height=50, units = 'cm', res = 200)
  DimHeatmap(seurat_stage[[stage]], dims = 1:30, balanced = TRUE, cells = 500)
  graphics.off()
  
  png(paste0(curr.plot.path, "elbowplot.", stage, ".png"), width=24, height=20, units = 'cm', res = 200)
  print(ElbowPlot(seurat_stage[[stage]], ndims = 40))
  graphics.off()
  
  png(paste0(curr.plot.path, "UMAP_PCA_comparison.", stage, ".png"), width=40, height=30, units = 'cm', res = 200)
  PCA.level.comparison(seurat_stage[[stage]], PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
  graphics.off()
}

# Use PCA=15 as elbow plot is relatively stable across stages
seurat_stage <- lapply(seurat_stage, function(x) FindNeighbors(x, dims = 1:15, verbose = FALSE))
seurat_stage <- lapply(seurat_stage, function(x) RunUMAP(x, dims = 1:15, verbose = FALSE))

# Find optimal cluster resolution
for(stage in names(seurat_stage)){
  png(paste0(curr.plot.path, "clustree.", stage, ".png"), width=70, height=35, units = 'cm', res = 200)
  clust.res(seurat.obj = seurat_stage[[stage]], by = 0.1)
  graphics.off()
}

# Use custom clustering resolution for each stage
res = c("hh4" = 0.5, "hh6" = 0.5, "ss4" = 0.3, "ss8" = 0.3)
seurat_stage <- lapply(names(res), function(x) FindClusters(seurat_stage[[x]], resolution = res[names(res) %in% x]))
names(seurat_stage) = names(res)

# Plot UMAP for clusters and developmental stage
for(stage in names(seurat_stage)){
  png(paste0(curr.plot.path, "UMAP.", stage, ".png"), width=20, height=10, units = 'cm', res = 200)
  clust.stage.plot(seurat_stage[[stage]])
  graphics.off()
}

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster at each stage
# lower FC threshold as some clusters have no DE genes
markers <- lapply(seurat_stage, function(x) FindAllMarkers(x, only.pos = T, logfc.threshold = 0.1))
top15 <- lapply(markers, function(x) x %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC))

for(stage in names(seurat_stage)){
  png(paste0(curr.plot.path, "HM.top15.DE.", stage, ".png"), height = 75, width = 100, units = 'cm', res = 200)
  tenx.pheatmap(data = seurat_stage[[stage]], metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
                selected_genes = unique(top15[[stage]]$gene), gaps_col = "seurat_clusters")
  graphics.off()
}

# Plot features listed below at each stage
GOI = list("hh4" = c("VGLL1", "EPAS1", "GRHL3", "MSX1", "DLX5", "GATA2",
                     "AATF", "MAFA", "ING5", "SETD2", "LIN28B", "YEATS4",
                     "TBXT", "EOMES", "ADMP"),
           "hh6" = c("DLX5", "SIX1", "GATA2", "MSX1", "BMP4", "GBX2", "SIX3", "SOX2", "SOX21"),
           "ss4" = c("SIX1", "EYA2", "CSRNP1", "PAX7", "WNT4", "SIX3", "OLIG2", "SOX2", "SOX21"),
           "ss8" = c("SIX1", "EYA2", "SOX10", "TFAP2A", "GBX2", "SIX3", "OLIG2", "SOX2", "SOX21"))

for(stage in names(GOI)){
  ncol = 3
  png(paste0(curr.plot.path, "UMAP_GOI.", stage, ".png"), width = ncol*10, height = 10*ceiling((length(unlist(GOI[names(GOI) %in% stage]))+1)/ncol), units = "cm", res = 200)
  print(multi.feature.plot(seurat_stage[[stage]], stage.name = stage, n.col = ncol, label = "", gene.list = unlist(GOI[names(GOI) %in% stage])))
  graphics.off()
}

# Change order or clusters for plotting dotplots
levels = list("hh4" = c(3,1,0,2), "hh6" = c(2,1,3,0), "ss4" = c(2,3,1,0), "ss8" = c(3,2,1,4,0))
for(stage in names(levels)){
  seurat_stage[[stage]]$seurat_clusters <- factor(seurat_stage[[stage]]$seurat_clusters, levels = unlist(levels[names(levels) %in% stage]))
}

# Plot dotplot to identify clusters
for(stage in names(GOI)){
  png(paste0(curr.plot.path, "dotplot.", stage, ".png"), width = 30, height = 12, units = "cm", res = 200)
  print(DotPlot(seurat_stage[[stage]], group.by = "seurat_clusters", features = unlist(GOI[names(GOI) %in% stage])) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  graphics.off()
}

# Save stage data after clustering
saveRDS(seurat_stage, paste0(rds.path, 'seurat_stage_out.RDS'))

# Read in RDS data if needed
# seurat_stage <- readRDS(paste0(rds.path, "seurat_stage_out.RDS"))

# Make list of clusters to subset
clust.sub = list("hh4" = c(0,1,2), "hh6" = c(0,3), "ss4" = c(0,1), "ss8" = c(0,1,4))

########## Subset neural cells from clear seurat data (norm.data.clustfilt.cc)

# Set plot path
curr.plot.path <- paste0(plot.path, "5_neural_subset/")
dir.create(curr.plot.path)

# Get cell IDs from each stage based on clusters to subset
cell.sub = unlist(lapply(names(clust.sub), function(x){
  rownames(seurat_stage[[x]]@meta.data)[seurat_stage[[x]]$seurat_clusters %in% unlist(clust.sub[names(clust.sub) %in% x])]
}))

# Subset neural cells
neural.seurat <- subset(norm.data.clustfilt.cc, cells = cell.sub)

# Re-run findvariablefeatures and scaling
neural.seurat <- FindVariableFeatures(neural.seurat, selection.method = "vst", nfeatures = 2000)

# Enable parallelisation
plan("multiprocess", workers = ncores)
options(future.globals.maxSize = 2000 * 1024^2)

neural.seurat <- ScaleData(neural.seurat, features = rownames(neural.seurat), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

saveRDS(neural.seurat, paste0(rds.path, "neural.seurat.RDS"))

# Read in RDS data if needed
# neural.seurat <- readRDS(paste0(rds.path, "neural.seurat.RDS"))

# PCA
neural.seurat <- RunPCA(object = neural.seurat, verbose = FALSE)

png(paste0(curr.plot.path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(neural.seurat, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(curr.plot.path, "elbowplot.png"), width=24, height=20, units = 'cm', res = 200)
print(ElbowPlot(neural.seurat, ndims = 40))
graphics.off()

png(paste0(curr.plot.path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
PCA.level.comparison(neural.seurat, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()

# Use PCA=15 as elbow plot is relatively stable across stages
neural.seurat <- FindNeighbors(neural.seurat, dims = 1:15, verbose = FALSE)
neural.seurat <- RunUMAP(neural.seurat, dims = 1:15, verbose = FALSE)

# Find optimal cluster resolution
png(paste0(curr.plot.path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
clust.res(seurat.obj = neural.seurat, by = 0.2)
graphics.off()

# Use clustering resolution = 1.2
neural.seurat <- FindClusters(neural.seurat, resolution = 1.2)

# Plot UMAP for clusters and developmental stage
png(paste0(curr.plot.path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
clust.stage.plot(neural.seurat)
graphics.off()

# Plot QC for each cluster
png(paste0(curr.plot.path, "cluster.QC.png"), width=36, height=14, units = 'cm', res = 200)
QC.plot(neural.seurat)
graphics.off()

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(neural.seurat, only.pos = T, logfc.threshold = 0.25)
# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = neural.seurat, col.to.sort = seurat_clusters, sort.by = orig.ident)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

png(paste0(curr.plot.path, 'HM.top15.DE.png'), height = 75, width = 100, units = 'cm', res = 200)
tenx.pheatmap(data = neural.seurat, metadata = c("seurat_clusters", "orig.ident"), custom_order_column = "seurat_clusters",
              custom_order = cluster.order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters")
graphics.off()


saveRDS(neural.seurat, paste0(rds.path, "neural.seurat.out.RDS"))


