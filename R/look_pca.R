#' Performs PCA, and tells you what percent of the variance each PC is. 
#' 
#' 
#' @export
#' 

look_pca <- function(no.outlier){
  no.outlier <- no.outlier[[1]]
  no.outlier$submitted_subject_id <- NULL
  pca <- prcomp(no.outlier, center = T)
  value = 20
  pca.variance.explained = pca$sdev^2 / sum(pca$sdev^2)
  barplot(100*pca.variance.explained[1:value], las=2, xlab='', ylab='% Variance Explained', main = paste0('Variance explained by first ', value,' PCs'))
  # makes a scree plot for your pca!
}

