# function for plotting gene module heatmap
tenx.pheatmap <- function(data, metadata, primary_ordering = metadata[1], secondary_ordering = NULL, tertiary_ordering = NULL,
                           custom_order = NULL, custom_order_column = primary_ordering, assay = "RNA", slot = "scale.data",
                           selected_genes,
                           main = '', hide_annotation = NULL, show_rownames = TRUE, hclust_rows = FALSE, hclust_cols = FALSE, gaps_col = NULL,
                           cell_subset = NULL, treeheight_row = 0, use_seurat_colours = TRUE,  colour_scheme = c("PRGn", "RdYlBu", "Greys")){
  
  if(!is.null(cell_subset)){
    data <- subset(data, cells = cell_subset)
  } else {}
  
  HM.col <- droplevels(data@meta.data[, metadata])
  
  if(!is.null(custom_order)){
    if(!all(as.factor(custom_order) %in% HM.col[[custom_order_column]])){
      stop("custom_order factors not all found in custom_order_column")
    } else {}
    HM.col[[custom_order_column]] <-  factor(HM.col[[custom_order_column]], levels = custom_order)
    HM.col <- HM.col[order(HM.col[[custom_order_column]]),,drop = FALSE]
    
    if(!all(custom_order %in% HM.col[[primary_ordering]]) & custom_order_column == primary_ordering){
      stop("custom_order_column needs to be specified")
    } else {}
  }
  
  # three ordering levels possible, primary, secondary and tertiary
  if(!is.null(tertiary_ordering)){
    HM.col <- HM.col[order(HM.col[,primary_ordering], HM.col[,secondary_ordering], HM.col[,tertiary_ordering]),,drop = FALSE]
  } else if(!is.null(secondary_ordering)){
    HM.col <- HM.col[order(HM.col[,primary_ordering], HM.col[,secondary_ordering]),,drop = FALSE]
  } else {
    HM.col <- HM.col[order(HM.col[[primary_ordering]]),,drop = FALSE]
  }
  
  # gaps_col specifies a metadata column which column gaps are calculated from
  if(!is.null(gaps_col)) {
    if(class(gaps_col) != "character"){
      stop("gaps_col must be a metadata column name")
    } else {
      gaps_col = cumsum(rle(as.vector(HM.col[[gaps_col]]))[["lengths"]])
    }
  } else {
  }
  
  # hide as many annotations in metadata as desired with hide_annotation
  if(!is.null(hide_annotation)){
    HM.col[,hide_annotation] <- NULL
  } else{}
  
  if(use_seurat_colours == FALSE){
    # set colours for metadata
    ann_colours <- list()
    for(tic in 1:ncol(HM.col)){
      ann_colours[[colnames(HM.col[tic])]] <- setNames(colorRampPalette(brewer.pal(9, colour_scheme[tic])[2:9])(length(unique(HM.col[,tic]))),
                                                       unique(HM.col[,tic]))
    }
  } else {
    # set colours ggplot default colours, as in Seurat::DimPlot
    ann_colours <- list()
    for(column in metadata){
      ann_colours[[column]] <- setNames(ggplotColours(n = length(levels(data@meta.data[,column]))), levels(data@meta.data[,column]))
      
      # change levels of HM col so that heatmap annotations are in the same order as plotted
      ann_colours[[column]] <- ann_colours[[column]][match(levels(HM.col[[column]]), names(ann_colours[[column]]))]
    }
  }
  
  # extract data to plot from seurat object
  new.dat <- t(as.matrix(x = GetAssayData(object = data, assay = assay, slot = slot)[selected_genes, rownames(HM.col), drop = FALSE]))
  if(!is.null(cell_subset)){
    cat("rescaling data as cells have been subset \n")
    new.dat <- t(scale(t(new.dat)))
  } else {}
  new.dat <- replace(new.dat, new.dat >= 2, 2)
  new.dat <- replace(new.dat, new.dat <= -2, -2)
  
  print(pheatmap(t(new.dat), color = PurpleAndYellow(),
                 cluster_rows = hclust_rows, cluster_cols = hclust_cols, show_colnames = F,
                 annotation_col = HM.col[,rev(metadata)], fontsize = 22, fontsize_row = 12, gaps_col = gaps_col,
                 main = main, show_rownames = show_rownames, annotation_colors = ann_colours, treeheight_row = treeheight_row))
}
