#' Filters transcripts based on expression levels.  Also filters out individuals who are PCA outliers
#' 
#' @export
#' 
#' @examples

identify_pca_outliers <- function(transcripts_t, name)
{
  no_runs <- dplyr::select(transcripts_t, -submitted_subject_id)
  outliers <- 2
  print('filtering transcripts')
  means <- apply(no_runs, 2, mean)
  keep <- which(means>5)
  exclude <- which(means <= 5)
  no_runs <- no_runs[,keep]
  print(paste0(length(keep), ' transcripts after removing less than 5 avg counts'))
  # remove skewed transcripts
  mmratio <-apply(no_runs,2,mean, na.rm = T)/apply(no_runs,2,median, na.rm = T)
  print(paste0(ncol(no_runs[,mmratio<2]), ' transcripts after skew filtering'))
  no_runs <- no_runs[,mmratio<2]
  while(length(outliers) != 0){
    print('doing pca')
    pca <- prcomp(no_runs, center = T)
    outliers <- find_outliers(pca)
    print(paste0('removed ', length(outliers), ' outliers'))
    remove <- which(rownames(no_runs) %in% outliers)
    if (length(remove != 0)) {
      no_runs <- no_runs[-remove,]
    } 
    #    png(file = paste0('images/', name, '_pca.png'))
    #    explain_variance(pca)
    #    dev.off()
  }
  no_runs$PC1 <- pca$x[,1]
  no_runs$PC2 <- pca$x[,2]
  no_runs$PC3 <- pca$x[,3]
  no_runs$PC4 <- pca$x[,4]
  no_runs$PC5 <- pca$x[,5]
  no_runs$PC6 <- pca$x[,6]
  no_runs$PC7 <- pca$x[,7]
  no_runs$PC8 <- pca$x[,8]
  no_runs$PC9 <- pca$x[,9]
  no_runs$PC10 <- pca$x[,10]
  no_runs$submitted_subject_id <- rownames(no_runs)
  return(no_runs)
}
