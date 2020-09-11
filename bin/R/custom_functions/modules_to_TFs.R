
# gene_annotations is a dataframe from the seurat object misc slot with gene ensembl IDs and corresponding names
# GO:0043565 = sequence specific DNA binding and GO:0003700 = DNA-binding transcription factor activity

modules_to_TFs <- function(gm, gene_annotations, GO_terms = c('GO:0003700', 'GO:0043565')){
  
  gene_annotations <- as.data.frame(gene_annotations)
  module_ids <- gene_annotations[gene_annotations$gene_name %in% unlist(gm),]
  
  ensembl = useMart("ensembl",dataset="ggallus_gene_ensembl")
  ensembl_list <- getBM(attributes=c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"), 
                        filters = 'ensembl_gene_id', 
                        values = module_ids$gene_ID, 
                        mart = ensembl)
  
  # subset genes based on transcription factor GO terms
  gene_subset <- list()
  for(go in GO_terms){
    gene_subset[[go]] <- ensembl_list$ensembl_gene_id[ensembl_list$go_id %in% go]
  }
  gene_subset <- unique(unlist(gene_subset))
  
  # get current gene names from ensembl IDs
  gene_sub <- module_ids$gene_name[module_ids$gene_ID %in% gene_subset]
  
  # filter gene modules by transcription factors
  gm = lapply(gm, function(x) x[x %in% gene_sub])
  
  # remove empty elements from list
  return(gm[lapply(gm, length)>0])
}

