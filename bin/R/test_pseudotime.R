library(Seurat)
library(ggplot2)
library(ggbeeswarm)
library(tidyverse)

plot_dat <- neural.seurat.out@meta.data[,'orig.ident', drop=F] %>%
  tibble::rownames_to_column(var = "cell_name") %>%
  mutate(rank = Embeddings(object = neural.seurat.out[["pca"]])[, 1]) %>%
  mutate(rank = rank(rank))
  
ggplot(plot_dat, aes(x = pc1, y = orig.ident, 
                                             colour = orig.ident)) +
  geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")




norm.expression <- GetAssayData(neural.seurat.out, slot = 'scale.data')

goi <- stack(list(early = c('ETV5', 'MGA', 'OTX2'), middle = c('GBX2', 'HEY1', 'MEIS2', "HIF1A", "MYCN", "SOX3"), late = c("ZNF423", "TCF7L2", "HESX1")))
colnames(goi) = c("gene_name", "category")

temp <- as.data.frame(t(as.matrix(norm.expression[rownames(norm.expression) %in% goi$gene_name,]))) %>%
  tibble::rownames_to_column(var = "cell_name") %>%
  dplyr::full_join(plot_dat) %>%
  pivot_longer(!c(cell_name, orig.ident, rank), names_to = "gene_name", values_to = "normalised_count") %>%
  dplyr::left_join(goi)


ggplot(temp, aes(x = rank, y = normalised_count, 
                     group = gene_name, colour = category)) +
  geom_smooth(method="auto", se=FALSE, fullrange=FALSE, level=0.95) +
  theme_classic()

  

# calculate and plot pseudotime using monocle spline curve

# BiocManager::install("monocle")
# library(monocle)


goi <- stack(list(early = c('ETV5', 'MGA', 'OTX2'), middle = c('GBX2', 'HEY1', 'MEIS2', "HIF1A", "MYCN"), late = c("ZNF423", "TCF7L2", "HESX1")))


goi = stack(list(early = c("CITED4", "MAFA", "OTX2", "MYCN"), middle = c("YEATS4", "SOX2", "ZIC3", "LMO1"), late = c("SOX1", "ZEB2", "IRX2", "HOXC6", "NKX6-2")))
colnames(goi) = c("gene_name", "category")


neural.seurat.out@meta.data$Pseudotime = Embeddings(object = neural.seurat.out[["pca"]])[, 1]


monocle_data <- as.CellDataSet(neural.seurat.out, assay = 'RNA')

monocle_data <- estimateSizeFactors(monocle_data)
monocle_data <- estimateDispersions(monocle_data)

cds_subset <- monocle_data[rownames(monocle_data) %in% goi$gene_name,]

newdata <- data.frame(Pseudotime = seq(min(monocle_data$Pseudotime), max(monocle_data$Pseudotime), length.out = 100)) 

test <- monocle::genSmoothCurves(cds_subset, cores=1, trend_formula = '~sm.ns(Pseudotime, df=3)', relative_expr = T, new_data = newdata)

test <- log(test+1)

ggplot_data <- as.data.frame(test) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  pivot_longer(!gene_name, names_to = "pseudotime") %>%
  transform(pseudotime = as.numeric(pseudotime)) %>%
  dplyr::left_join(goi)
  

ggplot(ggplot_data, aes(x = pseudotime, y = value, 
                 group = gene_name, colour = gene_name)) +
  facet_wrap(~ category) +
  geom_line() +
  theme_classic()

ggplot(ggplot_data, aes(x = pseudotime, y = value, 
                        group = gene_name, colour = gene_name)) +
  geom_line() +
  theme_classic()

rownames(monocle_data) %in% goi$gene_name




goi$gene_name %in% rownames(monocle_data)
