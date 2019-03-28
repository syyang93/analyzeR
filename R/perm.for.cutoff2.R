#' Function that permutes the dataset and finds an appropriate cutoff value.  (Compatible with Vamsee's run.all.lms)
#' 
#' @param tx_expr Expression matrix in form: [genes x samples]. Will be converted to a list(!) of gene-expr vectors (if not input as list)
#' @param gene.ids Character vector of gene IDs, corresponding to rows in `tx_expr` [genes x samples]
#' @param cov Regression covariates [cov x samples]
#' @param SCORE Main covariate to be permuted (not included in `cov`)
#' @param num.cores The number of cores you would like to use
#' @export
#' 
#' @return A 95% cutoff value, also a data frame with all the p-values from each permutation (to make the overlaid p-val qqplot)
#' 
#' @examples
#' perm.results <- perm.for.cutoff2(permutations = 100, my.list[[1]], my.list[[2]], my.list[[3]], my.list[[4]])

perm.for.cutoff2 <- function(permutations = 100, tx_expr, cov, gene.ids, SCORE, num.cores = 10) {
  suppressPackageStartupMessages(library(parallel))
  suppressPackageStartupMessages(library(stringr))
  print(paste0('Permutations running on: ', ncol(tx_expr),' genes.'))
  
  perm.res <- rep(NA, permutations) # min p-values by permutation
  all.pvals <- as.data.frame(matrix(nrow=ncol(tx_expr), 
                                    ncol=permutations))
  rownames(all.pvals) <- gene.ids
  
  for(i in 1:permutations) {
    print(paste0('Permutation #: ', i, sep=''))
    SCORE.perm <- sample(SCORE)
    lm.res <- run.all.lms(tx_expr, cov, gene.ids, SCORE.perm, omit.outlier = T, num.cores = num.cores)
    
    # Harvest results
    all.pvals[,i] = lm.res$pval
    perm.res[i] <- min(lm.res$pval)
    print(paste0('min p-val: ', perm.res[i], sep=''))
  }
  
  # Gather 95% cutoff
  cutoff_95pct <- sort(perm.res)[ceiling(0.05 * permutations)]
  print(paste0('95% Cutoff: ', formatC(cutoff_95pct, format = 'e', digits = 3), sep=''))
  return(list(cutoff_95pct, all.pvals))
}
