library(Seurat)
library(ggplot2)
# install.packages('ggbeeswarm')
library(ggbeeswarm)
library(tidyverse)
# BiocManager::install("monocle")
library(monocle)

output_path = "./output/R_results/plots/8_pseudotime/"
dir.create(output_path, recursive = TRUE)

neural.seurat.out <- readRDS('./output/R_results/RDS.files/neural.seurat.out.RDS')

# plot PCA 1 and 2
png(filename = paste0(output_path, "PCA.png"), width = 20, height = 16, units = "cm", res = 400)
DimPlot(neural.seurat.out, reduction = 'pca', group.by = "orig.ident")
graphics.off()

# PCA captures differences across stage
plot_dat <- neural.seurat.out@meta.data[,'orig.ident', drop=F] %>%
  tibble::rownames_to_column(var = "cell_name") %>%
  mutate(pc1= Embeddings(object = neural.seurat.out[["pca"]])[, 1])

png(filename = paste0(output_path, "PC1.png"), width = 20, height = 16, units = "cm", res = 400)
ggplot(plot_dat, aes(x = pc1, y = orig.ident, 
                                             colour = orig.ident)) +
  geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")
graphics.off()

# calculate and plot pseudotime using monocle spline curve
# Initialise monocle object and use PC1 position as measurement for pseudotime
neural.seurat.out@meta.data$Pseudotime = Embeddings(object = neural.seurat.out[["pca"]])[, 1]
monocle_data <- as.CellDataSet(neural.seurat.out, assay = 'RNA')
monocle_data <- estimateSizeFactors(monocle_data)
# monocle_data <- estimateDispersions(monocle_data)

# select genes to plot
goi = stack(list(early = c("GATA2", "ID3", "KLF5", "EPAS1", "TFAP2C", "PDLIM1", "ERNI"), late = c("SOX2", "SOX1", "ZEB2", "HOXC6", "NKX6-2", "HOXB1", "GLI2", "MEIS1", "SIX3")))
colnames(goi) = c("gene_name", "category")

cds_subset <- monocle_data[rownames(monocle_data) %in% goi$gene_name,]
newdata <- data.frame(Pseudotime = seq(min(monocle_data$Pseudotime), max(monocle_data$Pseudotime), length.out = 100)) 
test <- monocle::genSmoothCurves(cds_subset, cores=1, trend_formula = '~sm.ns(Pseudotime, df=3)', relative_expr = T, new_data = newdata)
test <- log(test+1)

ggplot_data <- as.data.frame(test) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  pivot_longer(!gene_name, names_to = "pseudotime") %>%
  transform(pseudotime = as.numeric(pseudotime)) %>%
  dplyr::left_join(goi)
  
png(filename = paste0(output_path, "pseudotime_facet.png"), width = 20, height = 10, units = "cm", res = 400)
ggplot(ggplot_data, aes(x = pseudotime, y = value, 
                 group = gene_name, colour = gene_name)) +
  facet_wrap(~ category) +
  geom_line() +
  theme_classic()
graphics.off()


png(filename = paste0(output_path, "pseudotime_timepoint.png"), width = 15, height = 10, units = "cm", res = 400)
ggplot(ggplot_data, aes(x = pseudotime, y = value, 
                        group = gene_name, colour = category)) +
  geom_line() +
  theme_classic()
graphics.off()

