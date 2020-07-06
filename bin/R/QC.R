# In order to be able to run the script from either Rstudio, local terminal, or cluster terminal, I add a switch which looks for command line arguments. This then sets the directory paths accordingly.
library('getopt')

# set arguments for Rscript
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer",
  'customFuncs', 'm', 2, "character",
  'networkGenes', 'd', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

list.files("./", recursive = T, full.names = T)

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
  
  plot.path = "./results/plots/"
  rds.path = "./results/RDS.files/"
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  
  ##################################
  # set path where data is located #
  ##################################
  data_path = "./alignmentOut/cellrangerCounts_renamed"
  
  # read all files from dir
  files <- list.files(data_path, recursive = T, full.names = T)
  # remove file suffix
  file.path <- dirname(files)[!duplicated(dirname(files))]
  # make dataframe with tissue matching directory
  tissue = c("hh4", "hh6", "ss4", "ss8")
  matches <- sapply(tissue, function(x) file.path[grep(pattern = x, x = file.path)])
  sample.paths <- data.frame(tissue = names(matches), path = matches, row.names = NULL)
  
} else if (opt$runtype == "nextflow"){
  cat('pipeling running through nextflow\n')
  
  sapply(list.files(opt$customFuncs, full.names = T), source)
  
  plot.path = "plots/"
  dir.create(plot.path, recursive = T)
  rds.path = "RDS.files/"
  dir.create(rds.path, recursive = T)
  
  # read all files from folder and keep only those from chr_edit
  files <- list.files("./", recursive = T, full.names = T)
  # remove file suffix
  file.path <- dirname(files)[!duplicated(dirname(files))]
  # make dataframe with tissue matching directory
  tissue = c("hh4", "hh6", "ss4", "ss8")
  matches <- sapply(tissue, function(x) file.path[grep(pattern = x, x = file.path)])
  sample.paths <- data.frame(tissue = names(matches), path = matches, row.names = NULL)
  
} else if (opt$runtype == "docker"){
  cat('R script running through docker\n')
  
  sapply(list.files('/home/bin/R/custom_functions/', full.names = T), source)
  
  plot.path = "/home/results/plots/"
  rds.path = "/home/results/RDS.files/"
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  
  data_path = "/home/alignmentOut/cellrangerCounts_renamed"
  # read all files from dir
  files <- list.files(data_path, recursive = T, full.names = T)
  # remove file suffix
  file.path <- dirname(files)[!duplicated(dirname(files))]
  # make dataframe with tissue matching directory
  tissue = c("hh4", "hh6", "ss4", "ss8")
  matches <- sapply(tissue, function(x) file.path[grep(pattern = x, x = file.path)])
  sample.paths <- data.frame(tissue = names(matches), path = matches, row.names = NULL)
}

# set number of cores to use for parallelisation
if(is.null(opt$cores)){ncores = 4}else{ncores= opt$cores}
cat(paste0("script ran with ", ncores, " cores\n"))

# Load packages - packages are stored within renv in the repository
reticulate::use_python('/usr/bin/python3.7')
library(Seurat)

library(future)
library(dplyr)
library(Antler)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(reshape2)
library(viridis)

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

# make dataframe with different filtering parameters which can be put into a loop for carrying out downstream analysis
filt_crit_df<-data.frame(gene_min = c(0, 750, 1000, 1250), gene_max = c(Inf, 7000, 6000, 5000), MT_max = c(Inf, 17.5, 15, 12.5))
rownames(filt_crit_df)<-c("unfilt", "low", "med", "high")

# save image of filter conditions table
filt_tab <- t(filt_crit_df)[c(3,1,2), 2:4]
rownames(filt_tab) <- c("% max MT content", "min gene count", "max gene count")
pdf(paste0(plot.path, "filter_conditions.pdf"), height = 1.5, width = 3.5)
grid.arrange(top=textGrob("Filter Conditions",gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.3, vjust = 1), tableGrob(filt_tab, theme = ttheme_minimal()))
graphics.off()
rm(filt_tab)


# Add columns of metadata for each set of filtering parameters in order to label which cells should be kept/removed
for(condition in rownames(filt_crit_df)){
  merged.data <- add.filter.to.meta(data = merged.data, filt.param.df = filt_crit_df, filt.param.name = condition)
}

# Stats on no. genes expressed in each cell at each stage and no. of genes in >3 cells - pre and post filering
filt.stats <- alternate.cols(gene.stats(data=merged.data, min.cell = 3, filter_type = "unfilt", group_var = "orig.ident"),
                             gene.stats(data=merged.data, min.cell = 3, filter_type = "low", group_var = "orig.ident"),
                             gene.stats(data=merged.data, min.cell = 3, filter_type = "med", group_var = "orig.ident"),
                             gene.stats(data=merged.data, min.cell = 3, filter_type = "high", group_var = "orig.ident"))

############# plot filtering statistics ############

# save image of remaining cells stats after filtering
filt_tab <- sapply(c("unfilt", "low", "med", "high"),
                   function(x) filt.stats["no. remaining cells", grep(x, colnames(filt.stats))])
rownames(filt_tab) <- sapply(strsplit(rownames(filt_tab), " "), "[", 1)
# plot table
pdf(paste0(plot.path, "remaining_cells.pdf"), height = 2, width = 3.5)
grid.arrange(top=textGrob("Remaining Cell Count",gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.3, vjust = 0.5), tableGrob(filt_tab, theme = ttheme_minimal()))
graphics.off()
# plot bar plot
filt_tab <- melt(filt_tab[1:4,])
g <- ggplot(filt_tab, aes(x = Var2, y = value, fill = Var1)) +
  xlab("Filter Condition") +
  ylab("Cell Count") +
  labs(fill = "Stage") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = viridis(4)) +
  ggtitle("Cell count after filtering") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
pdf(paste0(plot.path, "remaining_cells_bar.pdf"), height = 5, width = 7)
plot(g)
graphics.off()
rm(filt_tab)

# save image of how many cells are remaining after filtering
filt_tab <- sapply(c("unfilt", "low", "med", "high"),
                   function(x) filt.stats["median genes per cell", grep(x, colnames(filt.stats))])
rownames(filt_tab) <- sapply(strsplit(rownames(filt_tab), " "), "[", 1)
filt_tab <- melt(filt_tab[1:4,])

g <- ggplot(filt_tab, aes(x = Var2, y = value, group = Var1)) +
  xlab("Filter Condition") +
  ylab("Median Gene Count") +
  geom_line(aes(colour = Var1)) +
  scale_colour_manual(values = viridis(4)) +
  ggtitle("Median gene count per cell after filtering") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
pdf(paste0(plot.path, "median_rem_genes.pdf"), height = 5, width = 7)
plot(g)
graphics.off()
rm(filt_tab)


############# QC plots ############

# plot number of remaining cells and median genes per cell at simulated minimum cutoff thresholds
plots <- gene.count.sim(seurat_data = merged.data, sep.group = T, sep.group.var = "orig.ident")
pdf(paste0(plot.path, "gene.count.sim.pdf"))
print(plot_grid(plotlist = plots))
graphics.off()


# plot histograms of gene number after applying different thresholds
plots <- list()
temp.dat<-data.frame("gene_number" = colSums(as.matrix(merged.data@assays$RNA@counts) > 0))
plots$unfilt <- ggplot(temp.dat, aes(x = gene_number)) +
  xlab("Number of genes / cell") +
  geom_histogram(aes(y=..density..), alpha=0.5, colour = "#A9A9A9", fill = "#A9A9A9", bins = 50) +
  geom_density(alpha=0.6, colour = "#A9A9A9", fill = "#A9A9A9") +
  geom_vline(data=temp.dat, aes(xintercept=filt_crit_df$gene_max[2]) , linetype="dashed", colour = "red") +
  geom_text(aes(x=filt_crit_df$gene_max[2] - 125, label="low", y=0.0005), angle = 90, colour="red", size=3) +
  geom_vline(data=temp.dat, aes(xintercept=filt_crit_df$gene_max[3]) , linetype="dashed", colour = "red") +
  geom_text(aes(x=filt_crit_df$gene_max[3] - 125, label="med", y=0.0005), angle = 90, colour="red", size=3) +
  geom_vline(data=temp.dat, aes(xintercept=filt_crit_df$gene_max[4]) , linetype="dashed", colour = "red") +
  geom_text(aes(x=filt_crit_df$gene_max[4] - 125, label="high", y=0.0005), angle = 90, colour="red", size=3) +
  geom_vline(data=temp.dat, aes(xintercept=filt_crit_df$gene_min[2]), linetype="dashed", colour = "blue") +
  geom_text(aes(x=filt_crit_df$gene_min[2] - 125, label="low", y=0.0005), angle = 90, colour="blue", size=3) +
  geom_vline(data=temp.dat, aes(xintercept=filt_crit_df$gene_min[3]), linetype="dashed", colour = "blue") +
  geom_text(aes(x=filt_crit_df$gene_min[3] - 125, label="med", y=0.0005), angle = 90, colour="blue", size=3) +
  geom_vline(data=temp.dat, aes(xintercept=filt_crit_df$gene_min[4]), linetype="dashed", colour = "blue") +
  geom_text(aes(x=filt_crit_df$gene_min[4] - 125, label="high", y=0.0005), angle = 90, colour="blue", size=3) +
  ggtitle("Unfiltered cells") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

for(i in rownames(filt_crit_df)){
  if (i == "unfilt") {
    next
  } else {}
  if (i == "low") {
    title <- "Low filter"
  } else if (i == "med") {
    title <- "Medium filter"
  } else if (i == "high") {
    title <- "High filter"
  }
  colsum <- colSums(as.matrix(merged.data@assays$RNA@counts) > 0)
  temp.dat <- data.frame("gene_number" = colsum[colsum > filt_crit_df[i, "gene_min"] & colsum < filt_crit_df[i, "gene_max"]])
  plots[[i]] <- ggplot(temp.dat, aes(x = gene_number)) +
    xlab("Number of genes / cell") +
    geom_histogram(aes(y=..density..), alpha=0.5, colour = "#A9A9A9", fill = "#A9A9A9", bins = 50) +
    geom_density(alpha=0.6, colour = "#A9A9A9", fill = "#A9A9A9") +
    ggtitle(title) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
}
pdf(paste0(plot.path, "gene.count.hist.pdf"))
grid.arrange(plots$unfilt, plots$low, plots$med, plots$high, nrow = 3, 
             layout_matrix = cbind(c(1,1,2), cbind(c(1,1,3)), c(1,1,4)))
graphics.off()


# For each of the different filtering columns in the metadata plot the QC plots before filtering
for(condition in rownames(filt_crit_df)){
  temp.filt_crit_df<-filt_crit_df[condition,]
  # violin plots of cell gene no, UMI count and MT content post filtering
  pdf(paste0(plot.path, "filt_param_", condition, ".filtered_vioplot.pdf"),width=30,height=15)
  p1 = vio.plot(seurat_data = merged.data, y.dat = "nCount_RNA", x_split = condition, y.lab = "UMI Count", group.var = "orig.ident")
  p2 = vio.plot(seurat_data = merged.data, y.dat = "nFeature_RNA", x_split = condition, y.lab = "Gene Count", group.var = "orig.ident")
  p3 = vio.plot(seurat_data = merged.data, y.dat = "percent.mt", x_split = condition, y.lab = "Percent Mitochondrial Genes", group.var = "orig.ident")
  p = plot_grid(p1, p2, p3, ncol = 3)
  title <- ggdraw() + draw_label(paste0("Filtered Cells (<", temp.filt_crit_df[,"gene_min"], " or ", temp.filt_crit_df[,"gene_max"], " genes and >", temp.filt_crit_df[,"MT_max"], "% MT genes)"),
                                 fontface='bold', size = 20)
  print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  graphics.off()
  
  # same plot again, however this time the filtering criteria are combined and coloured instead so that you can see where the cells are on the original data
  pdf(paste0(plot.path, "filt_param_", condition, ".unfiltered_vioplot.pdf"),width=30,height=15)
  p1 = vio.plot(seurat_data = merged.data, y.dat = "nCount_RNA", dot.col = condition, y.lab = "UMI Count", fill.vio = F, group.var = "orig.ident")
  p2 = vio.plot(seurat_data = merged.data, y.dat = "nFeature_RNA", dot.col = condition, y.lab = "Gene Count", fill.vio = F, group.var = "orig.ident")
  p3 = vio.plot(seurat_data = merged.data, y.dat = "percent.mt", dot.col = condition, y.lab = "Percent Mitochondrial Genes", fill.vio = F, group.var = "orig.ident")
  p = plot_grid(p1, p2, p3, ncol = 3)
  title <- ggdraw() + draw_label(paste0("Unfiltered Cells (<", temp.filt_crit_df[,"gene_min"], " or ", temp.filt_crit_df[,"gene_max"], " genes and >", temp.filt_crit_df[,"MT_max"], "% MT genes)"),
                                 fontface='bold', size = 20)
  print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  graphics.off()
  
  # plot feature scatter plots before filtering so that I can see which cells will be removed
  pdf(paste0(plot.path, "filt_param_", condition, ".featscat_prefilt_UMI.MT.pdf"),width=15,height=8)
  plot1 <- FeatureScatter(merged.data, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = condition, pt.size = 0.1)
  plot2 <- FeatureScatter(merged.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = condition, pt.size = 0.1)
  print(CombinePlots(plots = list(plot1, plot2)))
  graphics.off()
}

########### Filter subset ##########

# Subset data based on the different set of filtering parameters and plot QC after filtering
filt.object.list <- list()
for(condition in rownames(filt_crit_df)){
  temp.filt_crit_df<-filt_crit_df[condition,]
  # subset data
  filt.data<-subset(merged.data, subset = c(nFeature_RNA > temp.filt_crit_df[,"gene_min"] & nFeature_RNA < temp.filt_crit_df[,"gene_max"] & percent.mt < temp.filt_crit_df[,"MT_max"]))
  
  # plot only remaining cells
  pdf(paste0(plot.path, "postfiltered_vioplot.pdf"),width=30,height=15)
  p1 = vio.plot(seurat_data = filt.data, y.dat = "nCount_RNA", y.lab = "UMI Count", group.var = "orig.ident")
  p2 = vio.plot(seurat_data = filt.data, y.dat = "nFeature_RNA", y.lab = "Gene Count", group.var = "orig.ident")
  p3 = vio.plot(seurat_data = filt.data, y.dat = "percent.mt", y.lab = "Percent Mitochondrial Genes", group.var = "orig.ident")
  p = plot_grid(p1, p2, p3, ncol = 3)
  title <- ggdraw() + draw_label(paste0("Post Filtering (<", temp.filt_crit_df[,"gene_min"], " or ", temp.filt_crit_df[,"gene_max"], " genes and >", temp.filt_crit_df[,"MT_max"], "% MT genes)"), fontface='bold', size = 20)
  print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  graphics.off()
  
  # replot feature scatter plots after filtering
  pdf(paste0(plot.path, "featscat.postfilt_UMI.MT.pdf"),width=15,height=8)
  plot1 <- FeatureScatter(filt.data, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = condition, pt.size = 0.1)
  plot2 <- FeatureScatter(filt.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = condition, pt.size = 0.1)
  print(CombinePlots(plots = list(plot1, plot2)))
  graphics.off()
  
  # save object in loop into the environment
  name<-paste0(condition, ".filt.data")
  filt.object.list[[name]]<-filt.data
}

# plot violin plots of all data after applying different filtering parameters - allows cross comparison of filtering levels
pdf(paste0(plot.path, "filter.comparison.vioplot.pdf"),width=30,height=15)
p1 = vio.plot.filt(seurat_data = merged.data, filt.param.df = filt_crit_df, y.dat = "nCount_RNA", y.lab = "UMI", x.lab = "Filtering Threshold")
p2 = vio.plot.filt(seurat_data = merged.data, filt.param.df = filt_crit_df, y.dat = "nFeature_RNA", y.lab = "Gene Count", x.lab = "Filtering Threshold")
p3 = vio.plot.filt(seurat_data = merged.data, filt.param.df = filt_crit_df, y.dat = "percent.mt", y.lab = "Percent Mitochondrial Genes", x.lab = "Filtering Threshold")
p = plot_grid(p1, p2, p3, ncol = 3)
title <- ggdraw() + draw_label(paste0("Filtering Threshold QC Comparison (all stages)"), fontface='bold', size = 20)
print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
graphics.off()



########### Normalise and scale data ##########

# Log normalize data and find variable features
norm.unfilt <- NormalizeData(filt.object.list[["unfilt.filt.data"]], normalization.method = "LogNormalize", scale.factor = 10000)
norm.unfilt <- FindVariableFeatures(norm.unfilt, selection.method = "vst", nfeatures = 2000)

# Scale data and regress out MT content
# Enable parallelisation on camp
if (opt$location == "CAMP") {
  plan("multiprocess", workers = ncores)
  options(future.globals.maxSize = 2000 * 1024^2)
} else {}
norm.unfilt <- ScaleData(norm.unfilt, features = rownames(norm.unfilt), vars.to.regress = "percent.mt")

# Save RDS after scaling as this step takes time
saveRDS(norm.unfilt, paste0(rds.path, "norm.unfilt.RDS"))
#norm.unfilt <- readRDS(paste0(rds.path, "norm.unfilt.RDS"))

########### PCA, UMAP and clustering ##########

# Run PCA analysis on the each set of data
norm.unfilt <- RunPCA(object = norm.unfilt, verbose = FALSE)

# Seurat's clustering algorithm is based on principle components, so we need to ensure that only the informative PCs are kept!
pdf(paste0(plot.path, "unfilt_dimHM.pdf"),width=15,height=25)
DimHeatmap(norm.unfilt, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

# another heuristic method is ElbowPlot which ranks PCs based on the % variance explained by each PC
pdf(paste0(plot.path, "unfilt_elbowplot.pdf"),width=12,height=10)
print(ElbowPlot(norm.unfilt, ndims = 40))
dev.off()

# FindNeighbors and FindClusters carry out the clustering
# FindNeighbors builds the SNN (shared nearest neighbour) graph
norm.unfilt <- FindNeighbors(norm.unfilt, dims = 1:15, verbose = FALSE)

# FindClusters runs community detection on the SNN graph - uses PCA data by default in order to identify clusters resolution arg can be altered to increase or reduce no. of clusters
norm.unfilt <- FindClusters(norm.unfilt, resolution = 0.5, verbose = FALSE)

# RunUMAP carries out non-linear dimensionality reduction for data visualisation
norm.unfilt <- RunUMAP(object = norm.unfilt, dims = 1:15, verbose = FALSE)

########### UMAP feature plots ##########

# UMAP plots of QC features to compare different filtering thresholds
p1 <- DimPlot(norm.unfilt, label = TRUE)
p2 <- DimPlot(norm.unfilt, group.by = "orig.ident", label = TRUE)
p3 <- FeaturePlot(norm.unfilt, features = "nFeature_RNA") + theme(plot.title = element_blank())
p4 <- FeaturePlot(norm.unfilt, features = "percent.mt") + theme(plot.title = element_blank())
p5 <- FeaturePlot(norm.unfilt, features = "nCount_RNA") + theme(plot.title = element_blank())
p6 <- DimPlot(norm.unfilt, group.by = "low")
p7 <- DimPlot(norm.unfilt, group.by = "med")
p8 <- DimPlot(norm.unfilt, group.by = "high")
lab = c("Clusters", "Developmental Stage", "Gene Count", "Percentage Mitochondrial Content", "UMI Count", "Low Filter Threshold",
        "Medium Filter Threshold", "High Filter Threshold")
p <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4, labels = lab, label_x = 0.1, hjust = 0, label_size = 17)
title <- ggdraw() + draw_label(paste0("UMAP plots of different QC metrics in the unfiltered normalised dataset"), fontface='bold', size = 20)
pdf(paste0(plot.path, "UMAP_filter_QC.pdf"),width=36,height=16)
print(plot_grid(title, p, ncol=1, rel_heights=c(0.2, 1)))
dev.off()

#################### DGEA ########################
# Differential expression test
markers <- FindAllMarkers(norm.unfilt, only.pos = T, logfc.threshold = 0.25)
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
top15 <- unique(subset(top15$gene, top15$gene %in% rownames(x = GetAssayData(object = norm.unfilt, slot = "scale.data"))))
clusters <- norm.unfilt[["seurat_clusters"]][order(norm.unfilt[["seurat_clusters"]]),, drop=FALSE]

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(norm.unfilt, only.pos = T, logfc.threshold = 0.25)

# get automated cluster order based on percentage of cells in adjacent stages
cluster.order = order.cell.stage.clust(seurat_object = norm.unfilt, col.to.sort = seurat_clusters, sort.by = orig.ident)

# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) %>% arrange(factor(cluster, levels = cluster.order))

# Specify colors
man_colours = list(
  low = c("Remaining Cells" = "#2e9112", "Filtered Cells"  = "#c80e0e") ,
  med = c("Remaining Cells" = "#2e9112", "Filtered Cells"  = "#c80e0e") ,
  high = c("Remaining Cells" = "#2e9112", "Filtered Cells"  = "#c80e0e"))

png(paste0(plot.path, 'HM.top15.DE.png'), height = 50, width = 75, units = 'cm', res = 200)
tenx.pheatmap(data = norm.unfilt, metadata = c("seurat_clusters", "orig.ident", "low", "med", "high"), custom_order = cluster.order, custom_order_column = "seurat_clusters",
              selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", cell_man_cols = man_colours)
graphics.off()

