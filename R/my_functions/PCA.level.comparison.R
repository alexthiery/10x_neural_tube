#### compare PCA levels ####
# loop through nested list containing Seurat objects for different filtering thresholds, developmental stages and PCA levels
PCA.level.comparison <- function(data, PCA.levels = c(10, 15, 20, 30), plot.path = curr.plot.path, cluster_res = 0.5, basename = "UMAP_PCA_comparison.pdf", height = 20, width = 15){
  if(class(data) == "list") {
    data.stage.PCA <- list()
    for (stage in names(data)){
      for (val in PCA.levels){
        name<-paste0(stage, " PCA ", val)
        data.stage.PCA[[name]] <- FindNeighbors(data[[stage]], dims = 1:val, verbose = FALSE)
        data.stage.PCA[[name]] <- FindClusters(data.stage.PCA[[name]], resolution = cluster_res, verbose = FALSE)
        data.stage.PCA[[name]] <- RunUMAP(object = data.stage.PCA[[name]], dims = 1:val, verbose = FALSE)
      }
    }
    plots <- lapply(seq_along(1:length(data.stage.PCA)), function(dat){
      DimPlot(data.stage.PCA[[dat]]) + ggtitle(paste(names(data.stage.PCA[dat])))
    })
    pdf(paste0(plot.path, basename), width = width, height = height)
    gridExtra::grid.arrange(grobs = plots)
    dev.off()
  } else {
    data.PCA <- list()
    for (val in PCA.levels){
      name<-paste0(" PCA ", val)
      data.PCA[[name]] <- FindNeighbors(data, dims = 1:val, verbose = FALSE)
      data.PCA[[name]] <- FindClusters(data.PCA[[name]], resolution = cluster_res, verbose = FALSE)
      data.PCA[[name]] <- RunUMAP(object = data.PCA[[name]], dims = 1:val, verbose = FALSE)
    }
    plots <- lapply(seq_along(1:length(data.PCA)), function(dat){
      DimPlot(data.PCA[[dat]]) + ggtitle(paste(names(data.PCA[dat])))
    })
    pdf(paste0(plot.path, basename), width = width, height = height)
    gridExtra::grid.arrange(grobs = plots)
    dev.off()
  }
}
