#' Filters transcripts based on expression levels.  
#' Also filters out individuals who are PCA outliers until there are no more PCA outliers - using find_outlier
#' 
#' @export
#' 


identify_pca_outliers <- function(transcripts_t, indiv_count)
{
  no_runs <- transcripts_t
  no_runs$submitted_subject_id <- NULL
  outliers <- 2
  print('filtering transcripts')
  means <- apply(no_runs, 2, mean)
  mean(means)
  keep <- which(means>50)
  exclude <- which(means <= 50)
  no_runs <- no_runs[,keep]
  print(paste0(length(keep), ' transcripts after removing less than 50 avg counts'))
  
  # Remove transcripts not expressed in >90% of individuals
  zeroes <- colSums(no_runs == 0)
  cutoff <- nrow(no_runs) * 0.1
  length(zeroes[which(zeroes > cutoff)])
  no_runs2 <- no_runs[,-which(zeroes > cutoff)]
  print(paste0(ncol(no_runs2), ' transcripts after removing transcripts not expressed in > 90% of individuals'))
  no_runs <- no_runs2

  while(length(outliers) != 0){
    print('doing pca')
    pca <- prcomp(no_runs, center = T)
    outliers <- find_outliers(pca)
    print(paste0('removed ', length(outliers), ' outliers'))
    remove <- which(rownames(no_runs) %in% outliers)
    if (length(remove != 0)) {
      no_runs <- no_runs[-remove,]
    } 
  }
  num_removed <- nrow(transcripts_t) - nrow(no_runs)
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
  
  indiv_count <- rbind(indiv_count, c(num_removed, 'PCA Outlier'))
  return(list(no_runs, indiv_count))
}