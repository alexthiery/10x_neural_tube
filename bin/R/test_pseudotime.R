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
  
  plot.path = "./results/R_results/plots/"
  rds.path = "./results/R_results/RDS.files/"
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


neural.seurat.out <- readRDS('./output/R_results/RDS.files/neural.seurat.out.RDS')

# plot_dat <- neural.seurat.out@meta.data[,'orig.ident', drop=F] %>%
#   tibble::rownames_to_column(var = "cell_name") %>%
#   mutate(rank = Embeddings(object = neural.seurat.out[["pca"]])[, 1]) %>%
#   mutate(rank = rank(rank))

# Plot cell stage along PC1
plot_dat <- neural.seurat.out@meta.data[,'orig.ident', drop=F] %>%
  tibble::rownames_to_column(var = "cell_name") %>%
  mutate(pc1 = Embeddings(object = neural.seurat.out[["pca"]])[, 1])


ggplot(plot_dat, aes(x = pc1, y = orig.ident, colour = orig.ident)) +
  geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")




norm.expression <- GetAssayData(neural.seurat.out, slot = 'scale.data')





# calculate and plot pseudotime using monocle spline curve

# BiocManager::install("monocle")
library(monocle)


# goi <- stack(list(early = c('ETV5', 'MGA', 'OTX2'), middle = c('GBX2', 'HEY1', 'MEIS2', "HIF1A", "MYCN"), late = c("ZNF423", "TCF7L2", "HESX1")))
# 
# 
# goi = stack(list(early = c("CITED4", "MAFA", "OTX2", "MYCN"), middle = c("YEATS4", "SOX2", "ZIC3", "LMO1"), late = c("SOX1", "ZEB2", "IRX2", "HOXC6", "NKX6-2")))
# colnames(goi) = c("gene_name", "category")


goi <- stack(list(early = c('ATF3', 'ETV1', 'NRIP1', 'RUNX1T1', 'TFAP2C', 'NAB1', 'NFKB2', 'CITED4', 'RNF14', 'JADE2', 'SAMD4A', 'SOX3', 'ERNI', 'TRIM24', 'Z-ZNF462'),
                  middle = c('ATF7IP', 'BRD8', 'ENSGALG00000037457', "FUBP1", "HESX1", "MBD3", 'MID1', 'PTTG2', 'SOX2', 'SP4', 'TCF7L2', 'ZNF821', 'MDFI', 'HMGA2', 'LIN28A', 'LMO1', 'LMX1B', 'PDCD4', 'SIX3'),
                  late = c("CBY1", "E2F3", "FEZF2", "GLI2", "HOXB1", "MAML2", "MEIS1", "NKX6-2", "Z-SMAD2Z", "T", "ZEB2", "SIX1")))

colnames(goi) = c("gene_name", "category")













neural.seurat.out@meta.data$Pseudotime = Embeddings(object = neural.seurat.out[["pca"]])[, 1]


monocle_data <- as.CellDataSet(neural.seurat.out, assay = 'RNA')

monocle_data <- estimateSizeFactors(monocle_data)
monocle_data <- estimateDispersions(monocle_data)

cds_subset <- monocle_data[rownames(monocle_data) %in% goi$gene_name,]

newdata <- data.frame(Pseudotime = seq(min(monocle_data$Pseudotime), max(monocle_data$Pseudotime), length.out = 100)) 

test <- genSmoothCurves(cds_subset, cores=1, trend_formula = '~sm.ns(Pseudotime, df=3)', relative_expr = T, new_data = newdata)

# test <- log(test+1)

ggplot_data <- as.data.frame(test) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  pivot_longer(!gene_name, names_to = "pseudotime") %>%
  transform(pseudotime = as.numeric(pseudotime)) %>%
  dplyr::left_join(goi)


# plot profile and facet genes by time
ggplot(ggplot_data, aes(x = pseudotime, y = value, 
                        group = gene_name, colour = gene_name)) +
  facet_wrap(~ category) +
  
  geom_line() +
  theme_classic()



ggplot_data <- ggplot_data[,-1] %>%
  group_by(pseudotime, category) %>%
  dplyr::summarize(mean = mean(value, na.rm=TRUE), sd=sd(value, na.rm=TRUE))

  
# plot profile and facet genes by time
ggplot(ggplot_data, aes(x = pseudotime, y = mean, colour = category)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - sd,
                  ymax = mean + sd), alpha = 0.2) +
  theme_classic()

# plot all genes and group by timepoint
ggplot(ggplot_data, aes(x = pseudotime, y = value, 
                        group = gene_name, colour = gene_name)) +
  geom_line() +
  theme_classic()

