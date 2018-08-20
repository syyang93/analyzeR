#' Reorganizes the object to be in the right format, prints number of genes to look at
#' 
#' @export
#' 

reorg <- function(full){
  look_genes <- colnames(full)
  ENSG_indices <- grep('ENSG', look_genes)
  print(length(ENSG_indices))
  not_ENSG_indices <- which(!(1:length(look_genes) %in% ENSG_indices))
  reorder <- full[,c(ENSG_indices, not_ENSG_indices)]
  return(reorder)
}
