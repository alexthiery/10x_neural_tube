
# custom function to plot heatmap from 10x data with extra user flexibility
tenx.pheatmap <- function(data, metadata, order.by = "seurat_clusters", custom_order = NULL, assay = "RNA", slot = "scale.data", used.genes = top15,
                          width = 28, height = 17, main = "Heatmap of top15 DE genes per cluster", plot.path = curr.plot.path,
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
    pdf(paste0(plot.path, basename, ".pdf"), width=width, height=height)
  } else if(file.format == "png"){
    png(paste0(plot.path, basename, ".png"), width=width*2.5, height=height*2.5, units = "cm", res = res)
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