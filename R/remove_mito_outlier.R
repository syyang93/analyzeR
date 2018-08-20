#' Remove individuals that are outliers in corrected mito %
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
  cutoff <- mean(na.omit(with.mt.unalign$mtDNA_adjust_AGE)) + 3*sd(na.omit(with.mt.unalign$mtDNA_adjust_AGE))
  
  cutoff_low <- mean(na.omit(with.mt.unalign$mtDNA_adjust_AGE)) - 3*sd(na.omit(with.mt.unalign$mtDNA_adjust_AGE))
  
  high_outliers <- which(with.mt.unalign$mtDNA_adjust_AGE > cutoff)
  
  with.mt.unalign$mtDNA_adjust_AGE[high_outliers]
  print(paste0('removed ', length(with.mt.unalign$mtDNA_adjust_AGE[high_outliers]), ' outliers in mtDNA-CN.'))
  with.mt.mito <- with.mt.unalign[-high_outliers,]
  return(with.mt.mito)
}
