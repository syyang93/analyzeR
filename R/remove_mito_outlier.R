# remove individuals that are outliers in corrected mito %

remove_mito_outlier <- function(blood.with.mt.unalign){
  cutoff <- mean(na.omit(blood.with.mt.unalign$mtDNA_adjust_AGE)) + 3*sd(na.omit(blood.with.mt.unalign$mtDNA_adjust_AGE))
  
  cutoff_low <- mean(na.omit(blood.with.mt.unalign$mtDNA_adjust_AGE)) - 3*sd(na.omit(blood.with.mt.unalign$mtDNA_adjust_AGE))
  
  high_outliers <- which(blood.with.mt.unalign$mtDNA_adjust_AGE > cutoff)
  
  blood.with.mt.unalign$mtDNA_adjust_AGE[high_outliers]
  print(paste0('removed ', length(blood.with.mt.unalign$mtDNA_adjust_AGE[high_outliers]), ' outliers in mtDNA-CN.'))
  blood.with.mt.mito <- blood.with.mt.unalign[-high_outliers,]
  return(blood.with.mt.mito)
}
