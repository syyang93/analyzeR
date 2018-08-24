#' Function that permutes the dataset and finds an appropriate cutoff value.
#' 
#' @param permutations The number of times you want to permute the data
#' @param lm_full The dataset you want to test
#' 
#' @export
#' 
#' @return A 95% cutoff value, also a data frame with all the p-values from each permutation (to make the overlaid p-val qqplot)
#' 
#' @examples
#' perm.results <- perm.for.cutoff(permutations = 5, subcut_full)


perm.for.cutoff <- function(permutations = 100, lm_full){
  permute_pvals<-as.data.frame(matrix(NA, nrow=permutations, ncol=2))
  colnames(permute_pvals) <- c('Permutation', 'Minimum pval')
  for_test <- lm_full[,1:length(grep('ENSG', colnames(lm_full)))] # the transcripts you're testing as your dependent variables
  
  # this data frame will contain all the pvalues for every permutation
  all_pvals<-as.data.frame(matrix(NA, nrow=ncol(for_test), ncol = permutations))
  for(i in 1:permutations)
  {
    print(paste0('on permutation ', i))
    lm_full$mtDNA_adjust_AGE <- sample(lm_full$mtDNA_adjust_AGE)
    lm_results <- testing_assoc(lm_full)
    
    # keeping minimum p-value in df
    permute_pvals[i, 1] <- i
    pvals <- lm_results$`Pr(>|t|)`
    permute_pvals[i, 2] <- min(pvals)
    all_pvals[,i] <- pvals
  }
  permute_pvals <- permute_pvals[order(permute_pvals$`Minimum pval`, decreasing = F),]
  cutoff_95_percent <- permute_pvals$`Minimum pval`[5]
  print(cutoff_95_percent)
  return(list(cutoff_95_percent, all_pvals))
}
