#' Function that will run a regression for your variable of interest across all genes, in parallel.  Adapted from Vamsee (github.com/vkp3/pillalamarRi)
#' 
#' @param tx_expr Expression matrix in form: [genes x samples]. Will be converted to a list(!) of gene-expr vectors (if not input as list)
#' @param gene.ids Character vector of gene IDs, corresponding to rows in `tx_expr` [genes x samples]
#' @param cov Regression covariates [cov x samples]
#' @param SCORE Main covariate to be permuted (not included in `cov`)
#' @param omit.outlier Whether or not you want to omit gene expression outliers
#' @param num.cores The number of cores you would like to use
#' @export 
#' 
#' @return All regression coefficients for lm(gene expression ~ SCORE + cov) for all genes
#' 
#' @examples
#' lm_res.sort <- run.all.lms(my.list[[1]], my.list[[2]], my.list[[3]], my.list[[4]], omit.outlier = T, num.cores = 10)

run.all.lms <- function(tx_expr, cov, gene.ids, SCORE, omit.outlier = T, num.cores = 10)
{
  require(pbapply)
  lm.res <-
    pblapply(tx_expr,            # Expression vector list for `pbapply::pblapply`
             run_lm,             # This function
             cov = cov,           # Covariate matrix, as desribed above
             SCORE = SCORE,       # PRS
             omit.outlier = omit.outlier,
             method = 'default',# Choose between 'default' or 'two-stage'
             cl = num.cores)      # Number of cores to parallelize over
  
  lm.res <- simplify2array(lm.res, higher=F)
  rownames(lm.res) <-
    c('intercept',
      'beta',
      'SE',
      't_value',
      'pval',
      'conf.low',
      'conf.high',
      'corr.rho')
  colnames(lm.res) <- gene.ids
 lm.res <- as.data.frame(t(lm.res))
  return(lm.res)
  # Sort results by p-value
  # No sorting! This will mess up the permutation saving
  # lm_res.sort <- lm.res[order(lm.res$pval), ]
  # return(lm_res.sort)
}
