
# cbind dataframes by rowname
cbind.rownames <- function(dat.1, dat.2, dat.3, dat.4){
if(rownames(dat.1) %in% rownames(dat.2) && rownames(dat.1) %in% rownames(dat.3) && rownames(dat.1) %in% rownames(dat.4)){
  dat.2 <- dat.2[match(rownames(dat.2), rownames(dat.1)),]
  dat.3 <- dat.3[match(rownames(dat.3), rownames(dat.1)),]
  dat.4 <- dat.4[match(rownames(dat.4), rownames(dat.1)),]
  return(cbind(dat.1, dat.2, dat.3, dat.4))
}
else{cat("check data frame rownames are equal")}
}

# plot quality metrics for each cluster
QC_plot <- function(seurat_obj, group_by = "seurat_clusters", y_elements = c("nCount_RNA", "nFeature_RNA", "percent.mt") ,
                    y_lab = c("UMI Count", "Gene Count", "% MT"), x_lab = "Cluster ID", basename = "cluster.QC", rds.path = rds.path, plot.path = curr.plot.path){
  plots <- list()
  plots$a = box.plot(dat = seurat_obj@meta.data, y_col = y_elements[1], group_by = group_by, y_lab = y_lab[1], x_lab = x_lab)
  plots$b = box.plot(dat = seurat_obj@meta.data, y_col = y_elements[2], group_by = group_by, y_lab = y_lab[2], x_lab = x_lab)
  plots$c = box.plot(dat = seurat_obj@meta.data, y_col = y_elements[3], group_by = group_by, y_lab = y_lab[3], x_lab = x_lab)
  
  pdf(file = paste0(plot.path, basename, ".pdf"), height = 7, width = 18)
  grid.arrange(grobs = plots, ncol = 3)
  dev.off()
}


# function to plot UMAP of cells co-expressing genes put as input to function
coexpressed_UMAP <- function(seurat_data, features){
  GOI <- t(GetAssayData(seurat_data, slot = "counts", assay = "RNA")[features,])
  coexp_cell_IDs <- rownames(GOI[apply(GOI, 1, function(x) !any(x==0)),])
  seurat_data@meta.data$new <- sapply(rownames(seurat_data@meta.data), function(x) ifelse(x %in% coexp_cell_IDs, "co-expressed", NA))
  if(all(is.na(seurat_data@meta.data$new)) == T){
    stop("no cells co-express all features")
  }
  DimPlot(seurat_data, group.by = "new") + 
    ggtitle(paste(c(features, "co-expressed"), collapse = " ")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    NoLegend()
}

# make all combinations of vector inserted
make_combinations <- function(features){
  combis <- vector("list", length(features))
  combi1 <- data.frame(features)
  colnames(combi1) <- 'X1'
  combis[[1]] <- combi1
  combis[2:length(features)] <- lapply(2:length(features), function(n) data.frame(t(combn(features, n))))
  do.call(plyr::rbind.fill, combis)
}



# simple merge column function for when using data from different dataframe sources
merge.col <- function(col1, col2, col3 = NA, col4 = NA, order.by = NA, add.colour.cols = F){
  
  if(add.colour.cols == T){
    myColors <- RColorBrewer::brewer.pal(length(unique(col1[,1])), "Pastel1")
    names(myColors) <- unique(col1[,1])
    col1$col1_hex <- unname(sapply(col1[,1], function(x) myColors[names(myColors) %in% x]))
    
    myColors <- RColorBrewer::brewer.pal(length(unique(col2[,1])), "Pastel2")
    names(myColors) <- unique(col2[,1])
    col2$col2_hex <- unname(sapply(col2[,1], function(x) myColors[names(myColors) %in% x]))
  }else{}
  dat <- merge(col1, col2, by = "row.names", all.x = T)
  rownames(dat) <- dat$Row.names
  dat$Row.names <- NULL
  
  if(!is.na(col3)){
    if(add.colour.cols == T){
      myColors <- RColorBrewer::brewer.pal(length(unique(col3[,1])), "Set3")
      names(myColors) <- unique(col3[,1])
      col3$col3_hex <- unname(sapply(col3[,1], function(x) myColors[names(myColors) %in% x]))
    }else{}
    dat <- merge(dat, col3, by = "row.names", all.x = T)
    rownames(dat) <- dat$Row.names
    dat$Row.names <- NULL
  } else{}
  if(!is.na(col4)){
    if(add.colour.cols == T){
      myColors <- RColorBrewer::brewer.pal(length(unique(col4[,1])), "Accent")
      names(myColors) <- unique(col4[,1])
      col4$col4_hex <- unname(sapply(col4[,1], function(x) myColors[names(myColors) %in% x]))
    }else{}
    dat <- merge(dat, col4, by = "row.names", all.x = T)
    rownames(dat) <- dat$Row.names
    dat$Row.names <- NULL
  } else{}
  if(!is.na(order.by)){
    # order cells by order.by
    dat <- dat[order(dat[[order.by]]),]
  }
  return(dat)
}




#### heatmap plots ####
# merge two dataframes together (clusters and filter condition)
HM.metadata.col <- function(dat, metadata, order.by = metadata[1]){
  if(length(metadata) == 0){
    stop("metadata not input")
  } else {
    HMcol <- dat@meta.data[, metadata[1], drop=F]
  }
  if(length(metadata) > 1){
    HMcol <- dat@meta.data[, metadata[1], drop=F]
    HMcol <- merge(dat@meta.data[, metadata[2], drop=F], HMcol, by="row.names", all.x = T)
    rownames(HMcol) <- HMcol$Row.names
    HMcol$Row.names <- NULL
  } else {}
  if(length(metadata) > 2){
      HMcol <- merge(dat@meta.data[, metadata[3], drop=F], HMcol, by="row.names", all.x = T)
      rownames(HMcol) <- HMcol$Row.names
      HMcol$Row.names <- NULL
    } else {}
    if(length(metadata) == 4){
      HMcol <- merge(dat@meta.data[, metadata[4], drop=F], HMcol, by="row.names", all.x = T)
      rownames(HMcol) <- HMcol$Row.names
      HMcol$Row.names <- NULL
    } else if (length(metadata) > 4){
      stop("Maximum 4 metadata columns")
    } else {}
  # order cells by clusters
  HMcol <- unclass(HMcol[order(HMcol[[order.by]]),,drop=F])
  return(as.data.frame(HMcol, row.names = attributes(HMcol)$row.names))
}







#### box plot for meta data stats ####
box.plot <- function(dat, y_col, group_by, y_lab, x_lab){
  ggplot(dat, aes(x = dat[[group_by]], y = dat[[y_col]], fill = dat[[group_by]])) +
    geom_boxplot() +
    xlab(x_lab) +
    ylab(y_lab) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.background = element_rect(colour = "white", fill = "white"),
      strip.text = element_text(size = 18),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18))
}


#### plot multiple feature plots at once ####
# if seurat object is not a list of objects then change multi.obj.list = F and
# stage = "ss8" or equivalent stage. Can also input cluster resolution that you want to plot the clusters
multi.feature.plot <- function(seurat.obj, stage.name = NULL, gene.list, multi.obj.list = T, plot.clusters = T, plot.stage = F, cluster.col = "seurat_clusters", cluster.res,
                               n.col = 4, plot_path = curr.plot.path, basename = "UMAP_GOI.pdf", label = "UMAP plots for GOI on normalised filtered data", legend.pos = "right", plot.celltype = F,
                               celltype.col = NA){
  if(plot.stage == F & plot.clusters == F){
    tot.len = length(gene.list)
  } else if (plot.stage == T & plot.clusters == T){
    tot.len = length(gene.list) + 2
  } else {
    tot.len = length(gene.list) + 1
  }
  
  if(multi.obj.list == T){
    for(i in names(seurat.obj)){
      if(plot.celltype == T){
        celltype.plot <- list(DimPlot(seurat.obj[[i]], group.by = celltype.col, label = T) + 
                             ggtitle(paste("Cell Types")) +
                             theme(plot.title = element_text(hjust = 0.5)) +
                             NoLegend())
      } else {
        celltype.plot <-NULL
      }
      if(plot.clusters == T){
        clust.plot <- list(DimPlot(seurat.obj[[i]], group.by = cluster.col) + 
                             ggtitle(paste("Clusters")) +
                             theme(plot.title = element_text(hjust = 0.5)))
      }else{
        clust.plot <- NULL
      }
      if(plot.stage == T){
        stage.plot <- list(DimPlot(seurat.obj[[i]], group.by =  "orig.ident") + 
                             ggtitle("Developmental Stage") +
                             theme(plot.title = element_text(hjust = 0.5)))
      }else{
        stage.plot <- NULL
      }
      plots <- lapply(gene.list, function(x) FeaturePlot(seurat.obj[[i]], features = x))
      plots <- c(stage.plot, celltype.plot, clust.plot, plots)
      pdf(paste0(plot_path, i, "/", basename),width = n.col*5, height = 5*ceiling(tot.len/n.col))
      print(gridExtra::grid.arrange(grobs = plots, ncol = n.col, top = textGrob(label = paste0(label, " (Stage = ", i, ")"), gp=gpar(fontsize=20, font = 2))))
      dev.off()
    }
  }else{
    if(plot.clusters == T){
      if(plot.celltype == T){
        celltype.plot <- list(DimPlot(seurat.obj, group.by = celltype.col, label = T) + 
                             ggtitle(paste("Cell Types")) +
                             theme(plot.title = element_text(hjust = 0.5)) +
                             NoLegend())
      } else {
        celltype.plot <-NULL
      }
      if(plot.clusters == T){
        clust.plot <- list(DimPlot(seurat.obj, group.by = cluster.col) + 
                             ggtitle(paste("Clusters")) +
                             theme(plot.title = element_text(hjust = 0.5)))
      }else{
        clust.plot <- NULL
      }
    if(plot.stage == T){
      stage.plot <- list(DimPlot(seurat.obj, group.by =  "orig.ident") + 
                           ggtitle("Developmental Stage") +
                           theme(plot.title = element_text(hjust = 0.5)))
    }else{
      stage.plot <- NULL
    }
    
    plots <- lapply(gene.list, function(x) FeaturePlot(seurat.obj, features = x))
    plots <- c(stage.plot, celltype.plot, clust.plot, plots)
    pdf(paste0(plot_path, "/", basename),width = n.col*5, height = 5*ceiling(tot.len/n.col))
    print(gridExtra::grid.arrange(grobs = plots, ncol = n.col, top = textGrob(label = label, gp=gpar(fontsize=20, font = 2))))
    dev.off()
  }
}
}


# Simple function to plot UMAP and dev.stage
clust.stage.plot <- function(data, cluster.col = "seurat_clusters", basename = "UMAP", plotpath = curr.plot.path){
  clust.plot <- DimPlot(data, group.by = cluster.col) + 
    ggtitle(paste("Clusters")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  stage.plot <- DimPlot(data, group.by =  "orig.ident") + 
    ggtitle("Developmental Stage") +
    theme(plot.title = element_text(hjust = 0.5))
  
  pdf(paste0(plotpath, basename, ".pdf"), width=10, height=5)
  gridExtra::grid.arrange(stage.plot, clust.plot, nrow = 1)
  dev.off()
}



#### clustree clustering ####
clust.res <- function(seurat.obj = progen_seurat, multi.obj.list = F, stage = NULL, plot_path = curr.plot.path, basename = "clustree.pdf", by = 0.1, starting_res = 0){
  plots <- list()
  resolutions <- c(seq(starting_res, starting_res+9*by, by=by))
  if(multi.obj.list == T){
    if(length(seurat.obj[[1]]@reductions) == 0){
      stop("Carry out dimensionality reduction (PCA) before clustering")
    }
    for(stage in names(seurat.obj)){
      seurat.obj[[stage]] <- seurat.obj[[stage]]@meta.data[,!grepl("RNA_snn_res.", colnames(seurat.obj[[stage]]@meta.data))]
      seurat.obj[[stage]] <- FindClusters(seurat.obj[[stage]], resolution = resolutions, verbose = F)
      plots[[stage]][["clustree"]] <- clustree(seurat.obj[[stage]]@meta.data, prefix = "RNA_snn_res.")
      for(res in resolutions[2:length(resolutions)]){
        plots[[stage]][[paste(res)]] <- DimPlot(seurat.obj[[stage]], group.by =  paste0("RNA_snn_res.", res)) +
          ggtitle(paste("resolution = ", res))
      }
      lay <- rbind(c(1,1,1,2,3,4),
                   c(1,1,1,5,6,7),
                   c(1,1,1,8,9,10))
      plots2 <- gridExtra::arrangeGrob(grobs = plots[[stage]], layout_matrix = lay)
      dev.off()
    }
  }else{
    if(length(seurat.obj@reductions) == 0){
      stop("Carry out dimensionality reduction (PCA) before clustering")
    }
    seurat.obj@meta.data <- seurat.obj@meta.data[,!grepl("RNA_snn_res.", colnames(seurat.obj@meta.data))]
    seurat.obj <- FindClusters(seurat.obj, resolution = resolutions, verbose = F)
    plots[["clustree"]] <- clustree(seurat.obj@meta.data, prefix = "RNA_snn_res.")
    for(res in resolutions[2:length(resolutions)]){
      plots[[paste(res)]] <- DimPlot(seurat.obj, group.by =  paste0("RNA_snn_res.", res)) +
        ggtitle(paste("resolution = ", res))
    }
    lay <- rbind(c(1,1,1,2,3,4),
                 c(1,1,1,5,6,7),
                 c(1,1,1,8,9,10))
    plots2 <- gridExtra::arrangeGrob(grobs = plots, layout_matrix = lay)
    dev.off()
  }
  pdf(paste0(plot_path, stage, basename), width= 25, height= 15)
  gridExtra::grid.arrange(plots2)
  dev.off()
  return(seurat.obj)
}




#### compare PCA levels ####
# loop through nested list containing Seurat objects for different filtering thresholds, developmental stages and PCA levels
PCA.level.comparison <- function(data, PCA.levels = c(10, 15, 20, 30), plot_path = curr.plot.path, cluster_res = 0.5, basename = "UMAP_PCA_comparison.pdf", height = 20, width = 15){
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
    pdf(paste0(plot_path, basename), width = width, height = height)
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
    pdf(paste0(plot_path, basename), width = width, height = height)
    gridExtra::grid.arrange(grobs = plots)
    dev.off()
  }
}




# custom function for violin plots in ggplot from Seurat metadata
# only a few arguments are customisable on purpose to keep the plots simple
# y.dat is the data you want to be plotted
# x_split is which feature you want to split the x axis by
# dot.col = factor to colour the dots by
# fill.vio = T/F as to whether to colour fill the violins
vio.plot <-
  function(meta.data,
           y.dat,
           y.lim = NULL,
           y.lab = y.dat,
           x_split = NULL,
           x.lab = NULL,
           dot.col = NULL,
           fill.vio = TRUE,
           group.var = "orig.ident"){
    meta.data <- as.data.frame(meta.data)
    if (fill.vio == TRUE) {
      col1 <- as.character(meta.data[[group.var]])
    } else{
      col1 <- NULL
    }
    p = ggplot(meta.data, aes(x = as.character(meta.data[[group.var]]), y = meta.data[[y.dat]], fill = col1)) +
      geom_violin() +
      xlab(x.lab) +
      ylab(y.lab) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill = NA, colour = NA),
        strip.background = element_rect(colour = "white", fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)
      )
    if (hasArg(dot.col)) {
      p = p + geom_jitter(
        shape = 21,
        position = position_jitter(0.2),
        size = 0.75,
        stroke = 0.05,
        aes(
          colour = relevel(factor(meta.data[[dot.col]]), ref = "Remaining Cells"),
          fill = relevel(factor(meta.data[[dot.col]]), ref = "Remaining Cells")
        )
      ) +
        scale_colour_manual(name = "colour", values = alpha(c("#000000", "#FF8C00", "#ff0000"), 0.5)) +
        scale_fill_manual(name = "fill", values = alpha(c("#000000", "#FF8C00", "#ff0000"), 0.1)) +
        guides(colour = guide_legend(override.aes = list(size = 3)))
    } else{
      p = p + geom_jitter(shape = 16,
                          position = position_jitter(0.2),
                          size = 0.1) +
        theme(legend.position = "none")
    }
    if (hasArg(y.lim)) {
      p = p + coord_cartesian(ylim = 1:y.lim)
    }
    if (hasArg(x_split)) {
      p + facet_wrap(~ meta.data[[x_split]])
    } else{
      p
    }
    return(p)
  }

# Alex gm function for correlation matrix without subsequent looping through function
subset.gm <- function(gm, fav.genes) {
  outlist <- list()
  for(mod in gm){
    if(any(fav.genes %in% mod)){
      gene.match <- str_c(fav.genes[which(fav.genes %in% mod)], collapse = "; ")
      outlist[[gene.match]] <- mod
    }
    else{next}
  }
  return(outlist)
}



# vio plot for filter comparison - doesnt split sample by stage, instead produces individual plot for each filtering parameter

vio.plot.filt <- function(meta.data,
                          filt.param.df,
                          y.dat,
                          y.lim = NULL,
                          y.lab = y.dat,
                          x.lab = NULL) {
  temp.sub <- lapply(rownames(filt.param.df), function(x) {
    meta.data[meta.data[[x]] == "Remaining Cells", y.dat]
  })
  names(temp.sub) <- rownames(filt.param.df)
  temp <- rbind.ident(temp.sub, data.name = y.lab)
  p = ggplot(temp, aes(x = temp[["ident"]], y = temp[[y.lab]], fill = temp[["ident"]])) +
    geom_violin() +
    xlab(x.lab) +
    ylab(y.lab) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.key = element_rect(fill = NA, colour = NA),
      strip.background = element_rect(colour = "white", fill = "white"),
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 18),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18)
    ) +
    geom_jitter(shape = 16,
                position = position_jitter(0.2),
                size = 0.1) +
    theme(legend.position = "none")
  if (hasArg(y.lim)) {
    p = p + coord_cartesian(ylim = 0:y.lim)
  } else {
    p
  }
  return(p)
}





# rbind lists of vectors into a single dataframe with a new column indicating their original identity
rbind.ident <- function(x, data.name){
  temp <- list()
  if(class(x) != "list") stop(
    "data is not of class: list"
  )
  if(is.null(names(x))){
    for(i in 1:length(x)){
      temp[[i]] <- data.frame(data.name = x[[i]], ident = i)
    }
  } else {
    for(i in names(x)){
      temp[[i]] <- data.frame(data.name = x[[i]], ident = i)
    }
  }
  
  temp <- do.call(rbind, temp)
  colnames(temp) <- c(data.name, "ident")
  rownames(temp) <- NULL
  return(temp)
}







# function for plotting gene module heatmap
GM.plot <- function(data, metadata, primary_ordering = metadata[1], secondary_ordering = NULL, tertiary_ordering = NULL,
                    custom_order = NULL, custom_order_column = primary_ordering, assay = "RNA", slot = "scale.data",
                    gene_modules = m$topCorr_DR$genemodules, fav_genes = NULL, colour_scheme = c("RdYlBu", "PRGn", "Greys"),
                    width = 28, height = 17, main = NULL, plotpath = curr.plot.path, hide_annotation = NULL,
                    basename = "HM.top15.DE", file.format = "png", res = 100, show_rownames = T, annotation_colors = NA,
                    hclust_rows = F, hclust_cols = F, gaps_row = T, gaps_col = NULL, gm_row_annotation = T, cell_subset = NULL){
  
  if(!is.null(cell_subset)){
    data <- subset(data, cells = cell_subset)
  } else {}
  
  HM.col <- HM.metadata.col(dat = data, metadata = metadata, order = primary_ordering)
  
  if(!is.null(custom_order)){
    if(!all(as.factor(custom_order) %in% HM.col[[custom_order_column]])){
      stop("custom_order factors not all found in custom_order_column")
    } else {}
    HM.col[[custom_order_column]] <-  factor(HM.col[[custom_order_column]], levels = custom_order)
    if(!all(custom_order %in% HM.col[[primary_ordering]]) & custom_order_column == primary_ordering){
      stop("custom_order_column needs to be specified")
    } else {}
  }
  
  if(!is.null(tertiary_ordering)){
    HM.col <- HM.col[order(HM.col[,primary_ordering], HM.col[,secondary_ordering], HM.col[,tertiary_ordering]),,drop = FALSE]
  } else if(!is.null(secondary_ordering)){
    HM.col <- HM.col[order(HM.col[,primary_ordering], HM.col[,secondary_ordering]),,drop = FALSE]
  } else {
    HM.col <- HM.col[order(HM.col[[primary_ordering]]),,drop = FALSE]
  }
  
  if(!is.null(fav_genes)){
    selected_GM <- subset.gm(gene_modules, fav_genes)
  } else {
    if(is.null(names(gene_modules))){
      names(gene_modules) <- paste0("GM:", 1:length(gene_modules))
      selected_GM <- gene_modules
    } else {}
    selected_GM <- gene_modules
  }
  
  new.dat <- t(as.matrix(x = GetAssayData(object = data, assay = assay, slot = slot)[unlist(selected_GM), rownames(HM.col), drop = FALSE]))
  if(!is.null(cell_subset)){
    cat("rescaling data as cells have been subset \n")
    new.dat <- t(scale(t(new.dat)))
  } else {}
  new.dat <- replace(new.dat, new.dat >= 2, 2)
  new.dat <- replace(new.dat, new.dat <= -2, -2)
  
  if(file.format == "pdf"){
    pdf(paste0(plotpath, basename, ".pdf"), width=width, height=height)
  } else if(file.format == "png"){
    png(paste0(plotpath, basename, ".png"), width=width*2.5, height=height*2.5, units = "cm", res = res)
  } else{
    stop("file format must be png or pdf")
  }
  
  if(gm_row_annotation == T) {
    row_ann <- stack(selected_GM)
    rownames(row_ann) <- row_ann$values
    colnames(row_ann)[2] <- "Gene Modules"
    row_ann$values <- NULL
  } else {
    row_ann = NA
  }
  
  if(gaps_row == T) {
    gaps_row = cumsum(summary(as.factor(row_ann[["Gene Modules"]]), maxsum = max(lengths(lapply(row_ann, unique)))))
  } else {
    gaps_row = NULL
  }
  
  if(!is.null(gaps_col)) {
    gaps_col = cumsum(as.vector(table(HM.col[[gaps_col]])))
  } else {
  }
  
  if(!is.null(hide_annotation)){
    HM.col[,hide_annotation] <- NULL
  } else{}
  
  ann_colours <- list()
  for(tic in 1:ncol(HM.col)){
    ann_colours[[colnames(HM.col[tic])]] <- setNames(colorRampPalette(brewer.pal(9, colour_scheme[tic])[2:9])(length(unique(HM.col[,tic]))),
                                                     unique(HM.col[,tic]))
  }
  
  ann_colours[["Gene Modules"]] <- setNames(colorRampPalette(brewer.pal(9, "Paired"))(length(unique(row_ann$`Gene Modules`))),
                                            unique(row_ann$`Gene Modules`))
  
  pheatmap(t(new.dat), color = PurpleAndYellow(),
           cluster_rows = hclust_rows, cluster_cols = hclust_cols, show_colnames = F,
           annotation_col = HM.col, fontsize = 22, fontsize_row = 12, gaps_col = gaps_col,
           gaps_row = gaps_row,
           main = main, show_rownames = show_rownames, annotation_row = row_ann, annotation_colors = ann_colours)
  dev.off()
}





tenx_pheatmap <- function(data, metadata, order.by = "seurat_clusters", custom_order = NULL, assay = "RNA", slot = "scale.data", used.genes = top15,
                          width = 28, height = 17, main = "Heatmap of top15 DE genes per cluster", plotpath = curr.plot.path,
                          basename = "HM.top15.DE", file.format = "png", res = 200, show_rownames = T, annotation_colors = NA,
                          hclust_rows = F, hclust_cols = F){
  
  HM.col <- HM.metadata.col(dat = data, metadata = metadata, order.by = order.by)
  HM.col <- HM.col[order(HM.col[[order.by]]),]
  if(!is.null(custom_order)){
    HM.col[[order.by]] <- factor(HM.col[[order.by]], levels = custom_order)
    HM.col <- HM.col[order(HM.col[[order.by]]),,drop = FALSE]
  } else {}
  
  new.dat <- t(as.matrix(x = GetAssayData(object = data, assay = assay, slot = slot)[used.genes,, drop = FALSE]))
  new.dat <- new.dat[rownames(HM.col), used.genes] 
  new.dat <- replace(new.dat, new.dat >= 2, 2)
  new.dat <- replace(new.dat, new.dat <= -2, -2)
  
  if(file.format == "pdf"){
    pdf(paste0(plotpath, basename, ".pdf"), width=width, height=height)
  } else if(file.format == "png"){
    png(paste0(plotpath, basename, ".png"), width=width*2.5, height=height*2.5, units = "cm", res = res)
  } else{
    stop("file format must be png or pdf")
  }
  
  pheatmap(t(new.dat), color = PurpleAndYellow(),
           cluster_rows = hclust_rows, cluster_cols = hclust_cols, show_colnames = F,
           annotation_col = HM.col, fontsize = 22, fontsize_row = 12, gaps_col = cumsum(as.vector(table(HM.col[[order.by]]))),
           annotation_colors = annotation_colors,
           main = main, show_rownames = show_rownames)
  dev.off()
}









tenx_pheatmap.old <- function(data, metadata, order.by = "seurat_clusters", custom_order = NULL, assay = "RNA", slot = "scale.data", used.genes = top15,
                              width = 28, height = 17, main = "Heatmap of top15 DE genes per cluster", plotpath = curr.plot.path,
                              basename = "HM.top15.DE", file.format = "png", res = 200, show_rownames = T, annotation_colors = NA,
                              hclust_rows = F, hclust_cols = F){
  
  HM.col <- HM.metadata.col(dat = data, metadata = metadata, order.by = order.by)
  HM.col <- HM.col[order(HM.col[[order.by]]),]
  if(!is.null(custom_order)){
    HM.col <-  HM.col %>%
      rownames_to_column('gene') %>%
      mutate(order.by = factor(order.by, levels = custom_order)) %>%
      arrange(order.by) %>%
      column_to_rownames('gene')
  } else {}
  
  new.dat <- t(as.matrix(x = GetAssayData(object = data, assay = assay, slot = slot)[used.genes,, drop = FALSE]))
  new.dat <- new.dat[rownames(HM.col), used.genes] 
  new.dat <- replace(new.dat, new.dat >= 2, 2)
  new.dat <- replace(new.dat, new.dat <= -2, -2)
  
  if(file.format == "pdf"){
    pdf(paste0(plotpath, basename, ".pdf"), width=width, height=height)
  } else if(file.format == "png"){
    png(paste0(plotpath, basename, ".png"), width=width*2.5, height=height*2.5, units = "cm", res = res)
  } else{
    stop("file format must be png or pdf")
  }
  
  pheatmap(t(new.dat), color = PurpleAndYellow(),
           cluster_rows = hclust_rows, cluster_cols = hclust_cols, show_colnames = F,
           annotation_col = HM.col, fontsize = 22, fontsize_row = 12, gaps_col = cumsum(as.vector(table(HM.col[[order.by]]))),
           annotation_colors = annotation_colors,
           main = main, show_rownames = show_rownames)
  dev.off()
}



