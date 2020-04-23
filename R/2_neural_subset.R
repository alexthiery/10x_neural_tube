# In order to be able to run the script from either Rstudio, local terminal, or cluster terminal, I add a switch which looks for command line arguments. This then sets the directory paths accordingly.

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  cat('no arguments provided\n')
  
  sapply(list.files('/Users/alex/dev/repos/10x_neural_tube/R/my_functions/', full.names = T), source)
  
  prev.rds.path = "/Users/alex/dev/output/10x_neural_tube/RDS.files/1_seurat_full/"
  
  plot.path = "/Users/alex/dev/output/10x_neural_tube/plots/2_neural_subset/"
  rds.path = "/Users/alex/dev/output/10x_neural_tube/RDS.files/2_neural_subset/"
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  
} else if (length(args)==1) {
  if (args[1] == "CAMP") {
    cat('data loaded from CAMP\n')
    
    project.dir = "~/working/alexthiery/analysis/10x_neural_tube/"
    
    sapply(list.files(paste0(project.dir, 'repo/scripts/my_functions/'), full.names = T), source)
    prev.rds.path = paste0(project.dir, "output/RDS.files/1_seurat_full/")
    
    plot.path = paste0(project.dir, "output/plots/2_neural_subset/")
    rds.path = paste0(project.dir, "output/plots/2_neural_subset/")
    dir.create(plot.path, recursive = T)
    dir.create(rds.path, recursive = T)
    

  } else {stop("Only CAMP can be supplied as arguments")}
} else {stop("only one argument can be supplied")}


# Load packages - packages are stored within renv in the repository

library(dplyr)
library(Antler)
library(Seurat)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)


all.seurat <- readRDS(paste0(prev.rds.path, "norm.data.clustfilt.cc.RDS"))

neural.seurat <- readRDS(paste0(prev.rds.path, "neural.seurat.out.RDS"))

hh4_genes <- read.table("~/dev/data/10x_neural/neural_induction_related_genes_HH4.txt", stringsAsFactors = F)

hh6_genes <- read.table("~/dev/data/10x_neural/neural_induction_related_genes_HH6.txt", stringsAsFactors = F)




# plot features for both neural subset and all data
pdf(paste0(plot.path, "neural_sub_hh4genes.featureplot.pdf"), width=10, height=100)
multi.feature.plot(neural.seurat, gene.list = hh4_genes$V1[hh4_genes$V1 %in% rownames(neural.seurat)])
graphics.off()

pdf(paste0(plot.path, "neural_sub_hh6genes.featureplot.pdf"), width=10, height=100)
multi.feature.plot(neural.seurat, gene.list = hh6_genes$V1[hh6_genes$V1 %in% rownames(neural.seurat)])
graphics.off()


pdf(paste0(plot.path, "allseurat_hh4genes.featureplot.pdf"), width=10, height=100)
multi.feature.plot(all.seurat, gene.list = hh4_genes$V1[hh4_genes$V1 %in% rownames(neural.seurat)])
graphics.off()

pdf(paste0(plot.path, "allseurat_hh6genes.featureplot.pdf"), width=10, height=100)
multi.feature.plot(all.seurat, gene.list = hh6_genes$V1[hh6_genes$V1 %in% rownames(neural.seurat)])
graphics.off()





# plot heatmap for both neural subset and all data
png(paste0(plot.path, "neural_sub_hh4genes.HM.png"), width=75, height=100, units = "cm", res = 200)
tenx.pheatmap(data = neural.seurat, metadata = c("orig.ident", "seurat_clusters"), order.by = "orig.ident",
              used.genes = hh4_genes$V1[hh4_genes$V1 %in% rownames(neural.seurat)], main = "", hclust_rows = T)
graphics.off()


png(paste0(plot.path, "neural_sub_hh6genes.HM.png"), width=75, height=100, units = "cm", res = 200)
tenx.pheatmap(data = neural.seurat, metadata = c("orig.ident", "seurat_clusters"), order.by = "orig.ident",
              used.genes = hh6_genes$V1[hh6_genes$V1 %in% rownames(neural.seurat)], main = "", hclust_rows = T)
graphics.off()



png(paste0(plot.path, "allseurat_hh4genes.HM.png"), width=75, height=100, units = "cm", res = 200)
tenx.pheatmap(data = all.seurat, metadata = c("orig.ident", "seurat_clusters"), order.by = "orig.ident",
              used.genes = hh4_genes$V1[hh4_genes$V1 %in% rownames(neural.seurat)], main = "", hclust_rows = T)
graphics.off()


png(paste0(plot.path, "allseurat_hh6genes.HM.png"), width=75, height=100, units = "cm", res = 200)
tenx.pheatmap(data = all.seurat, metadata = c("orig.ident", "seurat_clusters"), order.by = "orig.ident",
              used.genes = hh6_genes$V1[hh6_genes$V1 %in% rownames(neural.seurat)], main = "", hclust_rows = T)
graphics.off()

