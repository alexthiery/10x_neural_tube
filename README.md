---
title: '10x Chick Neural Tube'
author: 'Alex Thiery'
output:
  html_document:
    fig_width: 6
    fig_height: 4
    theme: simplex
    dev: jpeg
    keep_md: true
---
***
In order to be able to run the script from either Rstudio, local terminal, or cluster terminal, I add a switch which looks for command line arguments. This then sets the directory paths accordingly.


```{r, results='hide', message=FALSE, warning=FALSE}
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  cat('no arguments provided\n')
  
  sapply(list.files('/Users/alex/dev/repos/10x_neural_tube/R/my_functions/', full.names = T), source)
  
  plot.path = "/Users/alex/dev/output/10x_neural_tube/plots/1_seurat_full/"
  rds.path = "/Users/alex/dev/output/10x_neural_tube/RDS.files/1_seurat_full/"
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  
  ###### load data ##########
  sample.paths<-data.frame(tissue = c("hh4", "hh6", "ss4", "ss8"),
                           path = c("/Users/alex/dev/data/10x/cellranger_output/cellranger_count_hh4/outs/filtered_feature_bc_matrix_chr_edit/",
                                    "/Users/alex/dev/data/10x/cellranger_output/cellranger_count_hh6/outs/filtered_feature_bc_matrix_chr_edit/",
                                    "/Users/alex/dev/data/10x/cellranger_output/cellranger_count_4ss/outs/filtered_feature_bc_matrix_chr_edit/",
                                    "/Users/alex/dev/data/10x/cellranger_output/cellranger_count_8ss/outs/filtered_feature_bc_matrix_chr_edit/"))
} else if (length(args)==1) {
  if (args[1] == "CAMP") {
    cat('data loaded from CAMP\n')
    
    project.dir = "~/working/alexthiery/analysis/10x_scRNAseq_2019/"
    
    sapply(list.files(paste0(project.dir, 'repo/scripts/my_functions/'), full.names = T), source)

    plot.path = paste0(project.dir, "output/plots/1_seurat_full/")
    rds.path = paste0(project.dir, "output/RDS.files/1_seurat_full/")
    dir.create(plot.path, recursive = T)
    dir.create(rds.path, recursive = T)

    input.dat = paste0(project.dir, "cellranger_output/")

    ###### load data ##########
    sample.paths<-data.frame(tissue = c("hh4", "hh6", "ss4", "ss8"),
                             path = c('~/working/alexthiery/analysis/10x_scRNAseq_2019/cellranger_ouput/hh4_filtered_feature_bc_matrix_chr_edit/',
                                      '~/working/alexthiery/analysis/10x_scRNAseq_2019/cellranger_ouput/hh6_filtered_feature_bc_matrix_chr_edit/',
                                      '~/working/alexthiery/analysis/10x_scRNAseq_2019/cellranger_ouput/4ss_filtered_feature_bc_matrix_chr_edit/',
                                      '~/working/alexthiery/analysis/10x_scRNAseq_2019/cellranger_ouput/8ss_filtered_feature_bc_matrix_chr_edit/'))

  } else {stop("Only CAMP can be supplied as arguments")}
} else {stop("only one argument can be supplied")}
```

Load packages - packages are stored within renv in the repository

```{r, results='hide', message=FALSE, warning=FALSE}
library(dplyr)
library(Antler)
library(Seurat)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
```

Make Seurat objects for each of the different samples.
```{r, results='hide', message=FALSE, warning=FALSE}
for(i in 1:nrow(sample.paths["path"])){
  name<-paste(sample.paths[i,"tissue"])
  assign(name, CreateSeuratObject(counts= Read10X(data.dir = paste(sample.paths[i,"path"])), project = paste(sample.paths[i, "tissue"])))
}
```

The four Seurat objects are then merged, before running CreateSeuratObject again on the output in order to apply the min.cells parameter on the final merged dataset.
```{r, results='hide', message=FALSE, warning=FALSE}
temp <- merge(hh4, y = c(hh6, ss4, ss8), add.cell.ids = c("hh4", "hh6", "ss4", "ss8"), project = "chick.10x")
merged.data<-CreateSeuratObject(GetAssayData(temp), min.cells = 3, project = "chick.10x.mincells3")
```

The original Seurat objects are then removed from the global environment
```{r, results='hide', message=FALSE, warning=FALSE}
rm(hh4, hh6, ss4, ss8, sample.paths, temp)
```

Store mitochondrial percentage in object meta data
```{r, results='hide', message=FALSE, warning=FALSE}
merged.data <- PercentageFeatureSet(merged.data, pattern = "^MT-", col.name = "percent.mt")
```

### Filter data based on variable threshold
***

Remove data which do not pass filter threshold
```{r, results='hide', message=FALSE, warning=FALSE}
merged.data <- subset(merged.data, subset = c(nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 15))
```

Log normalize data and find variable features
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data <- NormalizeData(merged.data, normalization.method = "LogNormalize", scale.factor = 10000)
norm.data <- FindVariableFeatures(norm.data, selection.method = "vst", nfeatures = 2000)
```

Scale data and regress out MT content
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data <- ScaleData(norm.data, features = rownames(norm.data), vars.to.regress = "percent.mt")
```

Save RDS after scaling as this step takes time
```{r, results='hide', message=FALSE, warning=FALSE}
saveRDS(norm.data, paste0(rds.path, "norm.data.RDS"))
```

### Perform dimensionality reduction by PCA and UMAP embedding

Read in RDS data if needed
```{r, results='hide', message=FALSE, warning=FALSE}
# norm.data <- readRDS(paste0(rds.path, "norm.data.RDS"))
```

Change plot path
```{r, results='hide', message=FALSE, warning=FALSE}
curr.plot.path <- paste0(plot.path, '0_filt_data/')
dir.create(curr.plot.path)
```

Run PCA analysis on the each set of data
```{r}
norm.data <- RunPCA(object = norm.data, verbose = FALSE)
```


### Seurat's clustering algorithm is based on principle components, so we need to ensure that only the informative PCs are kept!
***
Plot heatmap of top variable genes across top principle components
```{r, results='hide', message=FALSE, warning=FALSE}
#pdf(paste0(curr.plot.path, "dimHM.pdf"),width=15,height=25)
DimHeatmap(norm.data, dims = 1:30, balanced = TRUE, cells = 500)
#graphics.off()
```

Another heuristic method is ElbowPlot which ranks PCs based on the % variance explained by each PC
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "elbowplot.pdf"),width=12,height=10)
print(ElbowPlot(norm.data, ndims = 40))
graphics.off()
```

Run clustering and UMAP at different PCA cutoffs - save this output to compare the optimal number of PCs to be used
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, 'UMAP_PCA_comparison.pdf'), width= 20, height= 15)
PCA.level.comparison(norm.data, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()
```

Use PCA=15 as elbow plot is relatively stable across stages
Use clustering resolution = 0.5 for filtering
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data <- FindNeighbors(norm.data, dims = 1:15, verbose = FALSE)
norm.data <- RunUMAP(norm.data, dims = 1:15, verbose = FALSE)
norm.data <- FindClusters(norm.data, resolution = 0.5, verbose = FALSE)
```

Plot UMAP for clusters and developmental stage
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "UMAP.pdf"), width=10, height=5)
clust.stage.plot(norm.data)
graphics.off()
```

Plot QC for each 
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "cluster.QC.pdf"), height = 7, width = 18)
QC.plot(norm.data)
graphics.off()
```

Find differentially expressed genes and plot heatmap of top DE genes for each cluster
```{r, results='hide', message=FALSE, warning=FALSE}
markers <- FindAllMarkers(norm.data, only.pos = T, logfc.threshold = 0.25)
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)

png(paste0(curr.plot.path, "HM.top15.DE.png"), width=75, height=50, units = "cm", res = 200)
tenx.pheatmap(data = norm.data, metadata = c("orig.ident", "seurat_clusters"), used.genes = unique(top15$gene), main = "")
graphics.off()
```

### Heatmap clearly shows clusters segregate by sex - check this and remove sex genes
***

Change plot path
```{r, results='hide', message=FALSE, warning=FALSE}
curr.plot.path <- paste0(plot.path, '1_sex_filt/')
dir.create(curr.plot.path)
```

There is a strong sex effect - this plot shows DE genes between clusters 1 and 2 which are hh4 clusters. Clustering is driven by sex genes
```{r, results='hide', message=FALSE, warning=FALSE}
png(paste0(curr.plot.path, "HM.top15.DE.png"), width=70, height=40, units = "cm", res = 200)
tenx.pheatmap(data = norm.data[,rownames(norm.data@meta.data[norm.data$seurat_clusters == 1 | norm.data$seurat_clusters == 2,])],
              metadata = c("orig.ident", "seurat_clusters"), used.genes = rownames(FindMarkers(norm.data, ident.1 = 1, ident.2 = 2)),
              main = "", hclust_rows = T)
graphics.off()
```

Use W chromosome genes to K-means cluster the cells into male (zz) and female (zw)
```{r, results='hide', message=FALSE, warning=FALSE}
W_genes <- as.matrix(norm.data@assays$RNA[grepl("W-", rownames(norm.data@assays$RNA)),])
k_clusters <- kmeans(t(W_genes), 2)
k_clusters <- data.frame(k_clusters$cluster)
norm.data@meta.data$k_clusters <- k_clusters[match(colnames(norm.data@assays$RNA), rownames(k_clusters)),]
```

Get rownames for kmeans clusters 1 and 2
```{r, results='hide', message=FALSE, warning=FALSE}
k_clus_1 <- rownames(norm.data@meta.data[norm.data@meta.data$k_clusters == 1,])
k_clus_2 <- rownames(norm.data@meta.data[norm.data@meta.data$k_clusters == 2,])
```

K clustering identities are stochastic, so I mneed to identify which cluster is male and female
Sum of W genes is order of magnitude greater in cluster 2 - these are the female cells
```{r, results='hide', message=FALSE, warning=FALSE}
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
```

Add sex data to meta.data
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data@meta.data$sex <- unlist(lapply(rownames(norm.data@meta.data), function(x)
  if(x %in% k.male){"male"} else if(x %in% k.female){"female"} else{stop("cell sex is not assigned")}))
```

#### Next - subset autosomal genes and Z genes - then calculate the average for each gene for both kmeans clustered cells and plot in order to compare whether the two groups show significant differences in their expression.

This is a test on autosomal genes to try and calculate and compare FC betweeen clusters
Calculating median is tricky as there are a lot of dropouts in 10x data so you end up with either 0s (when the median  = 0) or 1 (when the median expression in both clusters is the same - probably a result of normalisation resulting in a UMI of 0 or 1 being normalised to a nominal value)

Make dataframe for mean Z expression in male cells
```{r, results='hide', message=FALSE, warning=FALSE}
mean.Z.male <- data.frame(Z.mean = apply(norm.data@assays$RNA[grepl("Z-", rownames(norm.data@assays$RNA)), k.male], 1, mean))
```
Add 1 before log2 as log2(1) = 0
```{r, results='hide', message=FALSE, warning=FALSE}
mean.Z.male <- log2(mean.Z.male + 1)
```

Make dataframe for mean Z expression in female cells
```{r, results='hide', message=FALSE, warning=FALSE}
mean.Z.female <- data.frame(Z.mean = apply(norm.data@assays$RNA[grepl("Z-", rownames(norm.data@assays$RNA)), k.female], 1, mean))
mean.Z.female <- log2(mean.Z.female + 1)
```

Make dataframe for mean autosomal expression in male cells
```{r, results='hide', message=FALSE, warning=FALSE}
mean.auto.male <- data.frame(auto.mean = apply(norm.data@assays$RNA[!grepl("Z-", rownames(norm.data@assays$RNA)) & !grepl("W-", rownames(norm.data@assays$RNA)), k.male], 1, mean))
mean.auto.male <- log2(mean.auto.male + 1)
```

Make dataframe for mean autosomal expression in male cells
```{r, results='hide', message=FALSE, warning=FALSE}
mean.auto.female <- data.frame(auto.mean = apply(norm.data@assays$RNA[!grepl("Z-", rownames(norm.data@assays$RNA)) & !grepl("W-", rownames(norm.data@assays$RNA)), k.female], 1, mean))
mean.auto.female <- log2(mean.auto.female + 1)
```

Calculate FC by subtracting log2 expression from each other
```{r, results='hide', message=FALSE, warning=FALSE}
FC <- list()
FC$Z <- mean.Z.male - mean.Z.female
FC$auto <-  mean.auto.male - mean.auto.female
```

Plot boxplot of Z gene and autosomal expression in male vs female cells
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path,"sex_kmeans_log2FC_boxplot.pdf"))
boxplot(c(FC$Z, FC$auto),  ylab = "male - female log2 FC (mean normalised UMI +1)", names = c("Z chromosome genes", "autosomal genes"))
graphics.off()
```

Z genes are upregulated within male genes relative to female genes whereas autosomal genes have a normal distribution of logFC therefore Z genes should be filtered out - this could be done on a gene by gene basis, but can also be done through removing all Z genes

### Filter sex genes
***
Following subsetting of cells and/or genes, the same pipeline as above is repeated i.e.
Find variable features -> Scale data/regress out confounding variables -> PCA -> Find neighbours -> Run UMAP -> Find Clusters -> Cluster QC -> Find top DE genes

Remove sex genes
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.sexfilt <- norm.data[rownames(norm.data)[!grepl("W-", rownames(norm.data)) & !grepl("Z-", rownames(norm.data))], ]
```

Re-run findvariablefeatures and scaling
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.sexfilt <- FindVariableFeatures(norm.data.sexfilt, selection.method = "vst", nfeatures = 2000)
norm.data.sexfilt <- ScaleData(norm.data.sexfilt, features = rownames(norm.data.sexfilt), vars.to.regress = c("percent.mt", "sex"))
```

Save RDS
```{r, results='hide', message=FALSE, warning=FALSE}
saveRDS(norm.data.sexfilt, paste0(rds.path, "norm.data.sexfilt.RDS"))
```

Read in RDS data if needed
```{r, results='hide', message=FALSE, warning=FALSE}
# norm.data.sexfilt <- readRDS(paste0(rds.path, "norm.data.sexfilt.RDS"))
```

Set plot path
```{r, results='hide', message=FALSE, warning=FALSE}
curr.plot.path <- paste0(plot.path, '1_sex_filt/')
```

PCA
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.sexfilt <- RunPCA(object = norm.data.sexfilt, verbose = FALSE)

pdf(paste0(curr.plot.path, "dimHM.pdf"),width=15,height=25)
DimHeatmap(norm.data.sexfilt, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

pdf(paste0(curr.plot.path, "elbowplot.pdf"),width=12,height=10)
print(ElbowPlot(norm.data.sexfilt, ndims = 40))
graphics.off()

pdf(paste0(curr.plot.path, 'UMAP_PCA_comparison.pdf'), width= 20, height= 15)
PCA.level.comparison(norm.data.sexfilt, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()
```

Use PCA=15 as elbow plot is relatively stable across stages
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.sexfilt <- FindNeighbors(norm.data.sexfilt, dims = 1:15, verbose = FALSE)
norm.data.sexfilt <- RunUMAP(norm.data.sexfilt, dims = 1:15, verbose = FALSE)
```

Find optimal cluster resolution
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "clustree.pdf"), width= 25, height= 15, onefile = F)
clust.res(seurat.obj = norm.data.sexfilt, by = 0.2)
graphics.off()
```

Use clustering resolution = 1.4 for subsequent filtering of poor quality clusters this increases the stringency of poor quality clusters, removing the least data possible
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.sexfilt <- FindClusters(norm.data.sexfilt, resolution = 1.4, verbose = FALSE)
```

Plot UMAP for clusters and developmental stage
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "UMAP.pdf"), width=10, height=5)
clust.stage.plot(norm.data.sexfilt)
graphics.off()
```

Plot QC for each cluster
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "cluster.QC.pdf"), height = 7, width = 18)
QC.plot(norm.data.sexfilt)
graphics.off()
```

Find differentially expressed genes and plot heatmap of top DE genes for each cluster
```{r, results='hide', message=FALSE, warning=FALSE}
markers <- FindAllMarkers(norm.data.sexfilt, only.pos = T, logfc.threshold = 0.25)
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)

png(paste0(curr.plot.path, "HM.top15.DE.png"), width=100, height=75, units = "cm", res = 200)
tenx.pheatmap(data = norm.data.sexfilt, metadata = c("orig.ident", "seurat_clusters"), used.genes = unique(top15$gene), main = "")
graphics.off()
```

### Identify and remove poor quality clusters / contamination (mesoderm and PGCs)
***

Change plot path
```{r, results='hide', message=FALSE, warning=FALSE}
curr.plot.path <- paste0(plot.path, "2_cluster_filt/")
dir.create(curr.plot.path)
```

Identify mesoderm and PGCs
UMAP plots GOI
```{r, results='hide', message=FALSE, warning=FALSE}
genes <- c("EYA2", "SIX1", "TWIST1", "PITX2", "SOX17", "DAZL", "DND1", "CXCR4")

ncol = 3
pdf(paste0(curr.plot.path, "UMAP_GOI.pdf"), width = ncol*4, height = 5*ceiling(length(genes)/ncol))
multi.feature.plot(seurat.obj = norm.data.sexfilt, gene.list = genes, plot.clusters = T,
                   plot.stage = T, label = "", cluster.col = "RNA_snn_res.1.4", n.col = ncol)
graphics.off()
```

Dotplot for identifying PGCs, Early mesoderm and Late mesoderm
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "dotplot.GOI.pdf"), height = 8, width = 10)
DotPlot(norm.data.sexfilt, features = c( "SOX17", "CXCR4","EYA2", "TWIST1", "SIX1",  "PITX2", "DAZL"))
graphics.off()
```

### Remove contaminating cells from clusters
***

Clust 12, 14, 17 = poor quality cells
Clust 18 = early mesoderm - expresses sox17, eya2, pitx2, cxcr4
Clust 20 = Late mesoderm - expresses twist1, six1, eya2
Clust 21 = PGC's - expresses dazl very highly

```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.clustfilt <- rownames(norm.data.sexfilt@meta.data)[norm.data.sexfilt@meta.data$seurat_clusters == 12 |
                                                               norm.data.sexfilt@meta.data$seurat_clusters == 14 |
                                                               norm.data.sexfilt@meta.data$seurat_clusters == 17 |
                                                               norm.data.sexfilt@meta.data$seurat_clusters == 18 |
                                                               norm.data.sexfilt@meta.data$seurat_clusters == 20 |
                                                               norm.data.sexfilt@meta.data$seurat_clusters == 21]

norm.data.clustfilt <- subset(norm.data.sexfilt, cells = norm.data.clustfilt, invert = T)
```

Re-run findvariablefeatures and scaling
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.clustfilt <- FindVariableFeatures(norm.data.clustfilt, selection.method = "vst", nfeatures = 2000)
norm.data.clustfilt <- ScaleData(norm.data.clustfilt, features = rownames(norm.data.clustfilt), vars.to.regress = c("percent.mt", "sex"))

saveRDS(norm.data.clustfilt, paste0(rds.path, "norm.data.clustfilt.RDS"))
```

Read in RDS data if needed
```{r, results='hide', message=FALSE, warning=FALSE}
# norm.data.clustfilt <- readRDS(paste0(rds.path, "norm.data.clustfilt.RDS"))
```

PCA
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.clustfilt <- RunPCA(object = norm.data.clustfilt, verbose = FALSE)

pdf(paste0(curr.plot.path, "dimHM.pdf"),width=15,height=25)
DimHeatmap(norm.data.clustfilt, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

pdf(paste0(curr.plot.path, "elbowplot.pdf"),width=12,height=10)
print(ElbowPlot(norm.data.clustfilt, ndims = 40))
graphics.off()

pdf(paste0(curr.plot.path, 'UMAP_PCA_comparison.pdf'), width= 20, height= 15)
PCA.level.comparison(norm.data.clustfilt, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()
```

Use PCA=15 as elbow plot is relatively stable across stages
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.clustfilt <- FindNeighbors(norm.data.clustfilt, dims = 1:15, verbose = FALSE)
norm.data.clustfilt <- RunUMAP(norm.data.clustfilt, dims = 1:15, verbose = FALSE)
```

Find optimal cluster resolution
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "clustree.pdf"), width= 25, height= 15, onefile = F)
clust.res(seurat.obj = norm.data.clustfilt, by = 0.2)
graphics.off()
```

Use clustering resolution = 0.8
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.clustfilt <- FindClusters(norm.data.clustfilt, resolution = 0.8)
```

Plot UMAP for clusters and developmental stage
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "UMAP.pdf"), width=10, height=5)
clust.stage.plot(norm.data.clustfilt)
graphics.off()
```

### Check for cell cycle effect

Set plot path
```{r, results='hide', message=FALSE, warning=FALSE}
curr.plot.path <- paste0(plot.path, "3_cell_cycle/")
dir.create(curr.plot.path)
```

Calculate cell cycle for each cell
```{r, results='hide', message=FALSE, warning=FALSE}
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pre.cell.cycle.dat <- CellCycleScoring(norm.data.clustfilt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
```

Re-run findvariablefeatures and scaling
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.clustfilt.cc <- FindVariableFeatures(pre.cell.cycle.dat, selection.method = "vst", nfeatures = 2000)
norm.data.clustfilt.cc <- ScaleData(norm.data.clustfilt.cc, features = rownames(norm.data.clustfilt.cc), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

saveRDS(norm.data.clustfilt.cc, paste0(rds.path, "norm.data.clustfilt.cc.RDS"))
```

Read in RDS data if needed
```{r, results='hide', message=FALSE, warning=FALSE}
# norm.data.clustfilt.cc <- readRDS(paste0(rds.path, "norm.data.clustfilt.cc.RDS"))
```

PCA
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.clustfilt.cc <- RunPCA(object = norm.data.clustfilt.cc, verbose = FALSE)

pdf(paste0(curr.plot.path, "dimHM.pdf"),width=15,height=25)
DimHeatmap(norm.data.clustfilt.cc, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

pdf(paste0(curr.plot.path, "elbowplot.pdf"),width=12,height=10)
print(ElbowPlot(norm.data.clustfilt.cc, ndims = 40))
graphics.off()

pdf(paste0(curr.plot.path, 'UMAP_PCA_comparison.pdf'), width= 20, height= 15)
PCA.level.comparison(norm.data.clustfilt.cc, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()
```

Use PCA=15 as elbow plot is relatively stable across stages
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.clustfilt.cc <- FindNeighbors(norm.data.clustfilt.cc, dims = 1:15, verbose = FALSE)
norm.data.clustfilt.cc <- RunUMAP(norm.data.clustfilt.cc, dims = 1:15, verbose = FALSE)
```

Find optimal cluster resolution
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "clustree.pdf"), width= 25, height= 15, onefile = F)
clust.res(seurat.obj = norm.data.clustfilt.cc, by = 0.2)
graphics.off()
```

Use clustering resolution = 1.2
```{r, results='hide', message=FALSE, warning=FALSE}
norm.data.clustfilt.cc <- FindClusters(norm.data.clustfilt.cc, resolution = 1.2)
```

UMAP of cell cycle before and after regressing out
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "cell.cycle.pdf"), width = 14, height = 7)
pre.plot <- DimPlot(pre.cell.cycle.dat, group.by = "Phase", reduction = "umap")
post.plot <- DimPlot(norm.data.clustfilt.cc, group.by = "Phase", reduction = "umap")
print(gridExtra::grid.arrange(pre.plot, post.plot, ncol = 2))
graphics.off()
```

Plot UMAP for clusters and developmental stage
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "UMAP.pdf"), width=10, height=5)
clust.stage.plot(norm.data.clustfilt.cc)
graphics.off()
```

Plot QC for each cluster
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "cluster.QC.pdf"), height = 7, width = 18)
QC.plot(norm.data.clustfilt.cc)
graphics.off()
```

Find differentially expressed genes and plot heatmap of top DE genes for each cluster
```{r, results='hide', message=FALSE, warning=FALSE}
markers <- FindAllMarkers(norm.data.clustfilt.cc, only.pos = T, logfc.threshold = 0.25)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = c(1,2,11,7,0,6,4,8,9,10,3,5,12,13,14)))

png(paste0(curr.plot.path, "HM.top15.DE.png"), width=100, height=75, units = "cm", res = 200)
tenx.pheatmap(data = norm.data.clustfilt.cc, metadata = c("orig.ident", "seurat_clusters"), used.genes = unique(top15$gene),
              main = "", order.by = "seurat_clusters", custom_order = c(1,2,11,7,0,6,4,8,9,10,3,5,12,13,14))
graphics.off()
```

### Cell type identification
***

Set plot path
```{r, results='hide', message=FALSE, warning=FALSE}
curr.plot.path <- paste0(plot.path, "4_stage_split/")
dir.create(curr.plot.path)
```

Split dataset into different stages
```{r, results='hide', message=FALSE, warning=FALSE}
seurat_stage <- lapply(c('hh4', 'hh6', 'ss4', 'ss8'),
       function(x) subset(norm.data.clustfilt.cc, cells = rownames(norm.data.clustfilt.cc@meta.data)[norm.data.clustfilt.cc$orig.ident == x]))
names(seurat_stage) = c('hh4', 'hh6', 'ss4', 'ss8')
```

Re-run findvariablefeatures and scaling
```{r, results='hide', message=FALSE, warning=FALSE}
seurat_stage <- lapply(seurat_stage, function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000))
seurat_stage <- lapply(seurat_stage, function(x) ScaleData(x, features = rownames(norm.data.clustfilt.cc),
                                                           vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score")))

saveRDS(seurat_stage, paste0(rds.path, "seurat_stage.RDS"))
```

Read in RDS data if needed
```{r, results='hide', message=FALSE, warning=FALSE}
# seurat_stage <- readRDS(paste0(rds.path, "seurat_stage.RDS"))
```

PCA
```{r, results='hide', message=FALSE, warning=FALSE}
seurat_stage <- lapply(seurat_stage, function(x) RunPCA(object = x, verbose = FALSE))

for(stage in names(seurat_stage)){
  pdf(paste0(curr.plot.path, "dimHM.", stage, ".pdf"),width=15,height=25)
  DimHeatmap(seurat_stage[[stage]], dims = 1:30, balanced = TRUE, cells = 500)
  graphics.off()
  
  pdf(paste0(curr.plot.path, "elbowplot.", stage, ".pdf"),width=12,height=10)
  print(ElbowPlot(seurat_stage[[stage]], ndims = 40))
  graphics.off()
  
  pdf(paste0(curr.plot.path, "UMAP_PCA_comparison.", stage, ".pdf"), width= 20, height= 15)
  PCA.level.comparison(seurat_stage[[stage]], PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
  graphics.off()
}
```

Use PCA=15 as elbow plot is relatively stable across stages
```{r, results='hide', message=FALSE, warning=FALSE}
seurat_stage <- lapply(seurat_stage, function(x) FindNeighbors(x, dims = 1:20, verbose = FALSE))
seurat_stage <- lapply(seurat_stage, function(x) RunUMAP(x, dims = 1:20, verbose = FALSE))
```

Find optimal cluster resolution
```{r, results='hide', message=FALSE, warning=FALSE}
for(stage in names(seurat_stage)){
  pdf(paste0(curr.plot.path, "clustree.", stage, ".pdf"), width= 25, height= 15, onefile = F)
  clust.res(seurat.obj = seurat_stage[[stage]], by = 0.1)
  graphics.off()
}
```

Use custom clustering resolution for each stage
```{r, results='hide', message=FALSE, warning=FALSE}
res = c("hh4" = 0.5, "hh6" = 0.8, "ss4" = 0.3, "ss8" = 0.3)
seurat_stage <- lapply(names(res), function(x) FindClusters(seurat_stage[[x]], resolution = res[names(res) %in% x]))
names(seurat_stage) = names(res)
```

Plot UMAP for clusters and developmental stage
```{r, results='hide', message=FALSE, warning=FALSE}
for(stage in names(seurat_stage)){
  pdf(paste0(curr.plot.path, "UMAP.", stage, ".pdf"), width=10, height=5)
  clust.stage.plot(seurat_stage[[stage]])
  graphics.off()
}
```

Plot features listed below at each stage
```{r, results='hide', message=FALSE, warning=FALSE}
GOI = list("hh6" = c("DLX5", "SIX1", "GATA2", "MSX1", "BMP4", "SIX3", "GBX2", "SOX2", "SOX21"),
            "ss4" = c("SIX1", "EYA2", "CSRNP1", "PAX7", "WNT4", "SIX3", "OLIG2", "SOX2", "SOX21"),
            "ss8" = c("SIX1", "EYA2", "SOX10", "TFAP2A", "GBX2", "SIX3", "OLIG2", "SOX2", "SOX21"))

for(stage in names(GOI)){
  ncol = 3
  pdf(paste0(curr.plot.path, "UMAP_GOI.", stage, ".pdf"), width = ncol*4, height = 5*ceiling(length(genes)/ncol))
  multi.feature.plot(seurat_stage[[stage]], stage.name = stage, n.col = ncol, label = "", gene.list = unlist(GOI[names(GOI) %in% stage]))
  graphics.off()
}
```

Change order or clusters for plotting dotplots
```{r, results='hide', message=FALSE, warning=FALSE}
levels = list("hh6" = c(5,0,3,2,1,4), "ss4" = c(2,3,1,4,0), "ss8" = c(3,2,1,4,0))
for(stage in names(levels)){
  seurat_stage[[stage]]$seurat_clusters <- factor(seurat_stage[[stage]]$seurat_clusters, levels = unlist(levels[names(levels) %in% stage]))
}
```

Plot dotplot to identify clusters
```{r, results='hide', message=FALSE, warning=FALSE}
for(stage in names(GOI)){
  pdf(paste0(curr.plot.path, "dotplot.", stage, ".pdf"), width = 15, height = 6)
  print(DotPlot(seurat_stage[[stage]], group.by = "seurat_clusters", features = unlist(GOI[names(GOI) %in% stage])) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  graphics.off()
}
```

Make list of clusters to subset
```{r, results='hide', message=FALSE, warning=FALSE}
clust.sub = list("hh6" = c(2,1,4), "ss4" = c(0,1,4), "ss8" = c(0,1,4))
```

### Subset neural cells from clear seurat data (norm.data.clustfilt.cc)
***

Set plot path
```{r, results='hide', message=FALSE, warning=FALSE}
curr.plot.path <- paste0(plot.path, "5_neural_subset/")
dir.create(curr.plot.path)
```

Get cell IDs from each stage based on clusters to subset
```{r, results='hide', message=FALSE, warning=FALSE}
cell.sub = unlist(lapply(names(clust.sub), function(x){
  rownames(seurat_stage[[x]]@meta.data)[seurat_stage[[x]]$seurat_clusters %in% unlist(clust.sub[names(clust.sub) %in% x])]
}))
```

Subset neural cells
```{r, results='hide', message=FALSE, warning=FALSE}
neural.seurat <- subset(norm.data.clustfilt.cc, cells = cell.sub)
```

Re-run findvariablefeatures and scaling
```{r, results='hide', message=FALSE, warning=FALSE}
neural.seurat <- FindVariableFeatures(neural.seurat, selection.method = "vst", nfeatures = 2000)
neural.seurat <- ScaleData(neural.seurat, features = rownames(neural.seurat), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

saveRDS(neural.seurat, paste0(rds.path, "norm.data.clustfilt.cc.RDS"))
```

Read in RDS data if needed
```{r, results='hide', message=FALSE, warning=FALSE}
# neural.seurat <- readRDS(paste0(rds.path, "neural.seurat.RDS"))
```

PCA
```{r, results='hide', message=FALSE, warning=FALSE}
neural.seurat <- RunPCA(object = neural.seurat, verbose = FALSE)

pdf(paste0(curr.plot.path, "dimHM.pdf"),width=15,height=25)
DimHeatmap(neural.seurat, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

pdf(paste0(curr.plot.path, "elbowplot.pdf"),width=12,height=10)
print(ElbowPlot(neural.seurat, ndims = 40))
graphics.off()

pdf(paste0(curr.plot.path, 'UMAP_PCA_comparison.pdf'), width= 20, height= 15)
PCA.level.comparison(neural.seurat, PCA.levels = c(7, 10, 15, 20), cluster_res = 0.5)
graphics.off()
```

Use PCA=15 as elbow plot is relatively stable across stages
```{r, results='hide', message=FALSE, warning=FALSE}
neural.seurat <- FindNeighbors(neural.seurat, dims = 1:20, verbose = FALSE)
neural.seurat <- RunUMAP(neural.seurat, dims = 1:20, verbose = FALSE)
```

Find optimal cluster resolution
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "clustree.pdf"), width= 25, height= 15, onefile = F)
clust.res(seurat.obj = neural.seurat, by = 0.2)
graphics.off()
```

Use clustering resolution = 1.2
```{r, results='hide', message=FALSE, warning=FALSE}
neural.seurat <- FindClusters(neural.seurat, resolution = 1.2)
```

Plot UMAP for clusters and developmental stage
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "UMAP.pdf"), width=10, height=5)
clust.stage.plot(neural.seurat)
graphics.off()
```

Plot QC for each cluster
```{r, results='hide', message=FALSE, warning=FALSE}
pdf(paste0(curr.plot.path, "cluster.QC.pdf"), height = 7, width = 18)
QC.plot(neural.seurat)
graphics.off()
```

Find differentially expressed genes and plot heatmap of top DE genes for each cluster
```{r, results='hide', message=FALSE, warning=FALSE}
markers <- FindAllMarkers(neural.seurat, only.pos = T, logfc.threshold = 0.25)
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)

png(paste0(curr.plot.path, "HM.top15.DE.png"), width=100, height=75, units = "cm", res = 200)
tenx.pheatmap(data = neural.seurat, metadata = c("orig.ident", "seurat_clusters"), used.genes = unique(top15$gene),
              main = "", order.by = "seurat_clusters")
graphics.off()

saveRDS(neural.seurat, paste0(rds.path, "neural.seurat.out.RDS"))
```
