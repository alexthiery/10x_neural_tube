
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