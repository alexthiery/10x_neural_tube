#!/usr/bin/env Rscript

# In order to be able to run the script from either Rstudio, local terminal, or cluster terminal, I add a switch which looks for command line arguments. This then sets the directory paths accordingly.
library('getopt')

# set arguments for Rscript
spec = matrix(c(
  'location', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# set default location
if(is.null(opt$location)){opt$location = "local"}

if (opt$location == "local"){
  cat('Script running locally\n')
  
  sapply(list.files('/Users/alex/dev/repos/10x_neural_tube/R/my_functions/', full.names = T), source)
  
  prev.rds.path = "/Users/alex/dev/repos/10x_neural_tube/RDS.files/1_seurat_full/"
  
  plot.path = "/Users/alex/dev/repos/10x_neural_tube/output/plots/2_neural_subset/"
  rds.path = "/Users/alex/dev/repos/10x_neural_tube/output/RDS.files/2_neural_subset/"
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  
  # Read in favourite genes
  hh4_genes <- read.table("~/dev/data/10x_neural/neural_induction_related_genes_HH4.txt", stringsAsFactors = F)[,1]
  hh6_genes <- read.table("~/dev/data/10x_neural/neural_induction_related_genes_HH6.txt", stringsAsFactors = F)[,1]
  
} else if (opt$location == "CAMP"){
  cat('data loaded from CAMP\n')
    
  project.dir = "/camp/home/thierya/working/analysis/10x_neural_tube/"
  
  sapply(list.files(paste0(project.dir, 'repo/scripts/my_functions/'), full.names = T), source)
  prev.rds.path = paste0(project.dir, "output/RDS.files/1_seurat_full/")
  
  plot.path = paste0(project.dir, "output/plots/2_neural_subset/")
  rds.path = paste0(project.dir, "output/plots/2_neural_subset/")
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  
  # Read in favourite genes
  hh4_genes <- read.table("/camp/home/thierya/working/analysis/10x_neural_tube/repo/data/neural_induction_related_genes_HH4.txt", stringsAsFactors = F)[,1]
  hh6_genes <- read.table("/camp/home/thierya/working/analysis/10x_neural_tube/repo/data/neural_induction_related_genes_HH6.txt", stringsAsFactors = F)[,1]
  
  # set number of cores to use for parallelisation
  if(is.null(opt$cores)){ncores = 4}else{ncores= opt$cores}
  
  cat(paste0("script ran with ", ncores, " cores\n"))

} else {stop("Script can only be ran locally or on CAMP")}


# Load packages - packages are stored within renv in the repository
library(dplyr)
library(Antler)
library(Seurat)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)

# Plot UMAPs of favourite genes at each stage
seurat_stage <- readRDS(paste0(prev.rds.path, "seurat_stage_out.RDS"))
curr.plot.path = paste0(plot.path, "UMAPs_stage_data/")
dir.create(curr.plot.path)

# plot genes from hh4 gene list at each stage and zip
lapply(names(seurat_stage), function(x) plot.genes.zip(seurat_stage[[x]], hh4_genes, paste0(curr.plot.path, "hh4_genes_UMAPs/", x, "/")))
# plot genes from hh6 gene list at each stage and zip
lapply(names(seurat_stage), function(x) plot.genes.zip(seurat_stage[[x]], hh6_genes, paste0(curr.plot.path, "hh6_genes_UMAPs/", x, "/")))


# Plot UMAPs of favourite genes at neural subset
neural.seurat <- readRDS(paste0(prev.rds.path, "neural.seurat.out.RDS"))
curr.plot.path = paste0(plot.path, "UMAPs_neural_subset/")
dir.create(curr.plot.path)

# plot genes from hh4 gene list in neural subset and zip
plot.genes.zip(neural.seurat, hh4_genes, paste0(curr.plot.path, "hh4_genes_UMAPs/"))
# plot genes from hh6 gene list in neural subset and zip
plot.genes.zip(neural.seurat, hh4_genes, paste0(curr.plot.path, "hh6_genes_UMAPs/"))

# Plot heatmap for hh4 genes in neural subset
png(paste0(plot.path, "neural.seurat_hh4genes.HM.png"), width=75, height=100, units = "cm", res = 200)
tenx.pheatmap(data = neural.seurat, metadata = c("orig.ident", "seurat_clusters"), selected_genes = hh4_genes[hh4_genes %in% rownames(neural.seurat)],
              hclust_rows = T, gaps_col = "orig.ident", col_ann_order = c("orig.ident", "seurat_clusters"))
graphics.off()

# Plot heatmap for hh6 genes in neural subset
png(paste0(plot.path, "neural.seurat_hh6genes.HM.png"), width=75, height=100, units = "cm", res = 200)
tenx.pheatmap(data = neural.seurat, metadata = c("orig.ident", "seurat_clusters"), selected_genes = hh6_genes[hh6_genes %in% rownames(neural.seurat)],
              hclust_rows = T, gaps_col = "orig.ident", col_ann_order = c("orig.ident", "seurat_clusters"))
graphics.off()


