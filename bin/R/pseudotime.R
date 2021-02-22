#!/usr/bin/env Rscript

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
  
  plot.path = "./output/R_results/plots/8_pseudotime/"
  rds.path = "./output/R_results/RDS.files/8_pseudotime/"
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  
} else if (opt$runtype == "nextflow"){
  cat('pipeling running through nextflow\n')
  
  sapply(list.files(opt$customFuncs, full.names = T), source)
  
  plot.path = "plots/"
  dir.create(plot.path, recursive = T)
  rds.path = "RDS.files/"
  dir.create(rds.path, recursive = T)
}
  
 
# set number of cores to use for parallelisation
if(is.null(opt$cores)){ncores = 4}else{ncores= opt$cores}
cat(paste0("script ran with ", ncores, " cores\n"))

# Load packages - packages are stored within renv in the repository
reticulate::use_python('/usr/bin/python3.7')
library(Seurat)

library(ggplot2)
library(ggbeeswarm)
library(tidyverse)


# Plot profile for selected early, middle and late genes from network

# Read in seurat object
neural_seurat_out <- readRDS('./output/R_results/RDS.files/neural.seurat.out.RDS')

# Rank cells according to position on PC1
pc1_rank <- neural_seurat_out@meta.data[,'orig.ident', drop=F] %>%
  tibble::rownames_to_column(var = "cell_name") %>%
  mutate(pc1 = Embeddings(object = neural_seurat_out[["pca"]])[, 1]) %>%
  mutate(rank = rank(pc1)) %>%
  mutate(pseudotime = rank*(100/max(rank)))

# Plot cell stage along PC1
png(paste0(plot.path, 'pc1.png'), height = 18, width = 26, units = 'cm', res = 400)
ggplot(pc1_rank, aes(x = pc1, y = orig.ident, colour = orig.ident)) +
  geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")
graphics.off()

# Get gene names and timepoint first expressed - filter 0 and down-regulated timepoints
goi <- read.csv('./input_files/network_expression.csv') %>%
  filter(!timepoint == 0)

# Filter scaled seurat counts by goi -> generate long df with corresponding pseudotime and normalised count
plot_data <- as.data.frame(t(as.matrix(GetAssayData(neural_seurat_out, slot = 'scale.data')[rownames(neural_seurat_out) %in% goi$gene_name,]))) %>%
  tibble::rownames_to_column(var = "cell_name") %>%
  dplyr::full_join(pc1_rank) %>%
  pivot_longer(!c(cell_name, orig.ident, pseudotime, pc1, rank), names_to = "gene_name", values_to = "scaled_expression") %>%
  dplyr::left_join(goi) %>%
  droplevels() %>%
  mutate(timepoint = factor(timepoint, levels = c(1,3,5,7,9,12)))


# mean and SE summary data
plot_data_summary <- plot_data %>%
  mutate(rank_bin = pseudotime - (pseudotime %% 2.5)) %>% 
  group_by(rank_bin, timepoint) %>% 
  summarise(mn = mean(scaled_expression),
            se = sd(scaled_expression)/sqrt(n()))


# plot gam for all stages without standard error
png(paste0(plot.path, 'gam_pseudotime_allnetwork.png'), height = 18, width = 26, units = 'cm', res = 400)
ggplot(plot_data, aes(x = pseudotime, y = scaled_expression, colour = timepoint)) +
  geom_smooth(method="gam", se=FALSE) +
  theme_classic()
graphics.off()

# plot gam for all stages with standard error
png(paste0(plot.path, 'gam_se_pseudotime_allnetwork.png'), height = 18, width = 26, units = 'cm', res = 400)
ggplot(plot_data, aes(x = pseudotime, y = scaled_expression, colour = timepoint)) +
  geom_errorbar(data = plot_data_summary, aes(x = rank_bin, y = mn, 
                                 ymax = mn + se, ymin = mn - se), width = 2) +
  geom_point(data = plot_data_summary, aes(x = rank_bin, y = mn)) +
  geom_smooth(method="gam", se=FALSE) +
  theme_classic()
graphics.off()