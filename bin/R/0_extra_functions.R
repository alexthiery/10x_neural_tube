
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













