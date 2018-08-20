#' Remove individuals that are outliers in corrected mito % (Greater or lower than 3 SD from the mean)
#' Greater than 3SD from the mean
#' 
#' @param with.mt.unalign Data frame after filtering? 
#' @param with_gtex data frame with GTEx IDs
#' 
#' @export
#' 
#' @return Transcripts for those individuals
#' 
#' @examples


remove_mito_outlier <- function(with.mt.unalign){
  indiv.count <- with.mt.unalign[[2]]
  with.mt.unalign <- with.mt.unalign[[1]]
  cutoff <- mean(na.omit(with.mt.unalign$mtDNA_adjust_AGE)) + 3*sd(na.omit(with.mt.unalign$mtDNA_adjust_AGE))
  cutoff_low <- mean(na.omit(with.mt.unalign$mtDNA_adjust_AGE)) - 3*sd(na.omit(with.mt.unalign$mtDNA_adjust_AGE))
  
  high_outliers <- which(with.mt.unalign$mtDNA_adjust_AGE > cutoff)
  low_outliers <- which(with.mt.unalign$mtDNA_adjust_AGE < cutoff_low)
  high_outliers <- c(high_outliers, low_outliers)
  
  with.mt.unalign$mtDNA_adjust_AGE[high_outliers]
  print(paste0('removed ', length(with.mt.unalign$mtDNA_adjust_AGE[high_outliers]), ' outliers in mtDNA-CN.'))
  if(length(high_outliers) != 0) {with.mt.mito <- with.mt.unalign[-high_outliers,]} else{with.mt.mito <- with.mt.unalign}
  indiv.count <- rbind(indiv.count, c(length(with.mt.unalign$mtDNA_adjust_AGE[high_outliers]), 'Outlier in mtDNA-CN'))
  return(list(with.mt.mito, indiv.count))
}
