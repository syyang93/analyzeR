#' Function that will get the variables you need to run all the lms in parallel 
#' 
#' @param lm_full A data frame with all variables.  Column = genes, covariates, etc. Rows = observations (in this case, individuals)
#' @param permutevar The main variable, as a string (for my case, mtDNA-CN)
#' @param covnames Names of the covariates to be included in the regression, a vector of strings.
#' @export
#' 
#' @return a list containing tx_expr, cov, gene.ids, and SCORE, all of which are to be fed into the function "run.all.lms"
#' 
#' @examples
#' my.list <- get.vars.runall.lms(lm_full, permutevar = 'mtDNA_adjust_AGE')

run.all.lms <- function(tx_expr, cov, gene.ids, SCORE, num.cores = 10)
{
  require(pbapply)
  lm.res <-
    pblapply(tx_expr,            # Expression vector list for `pbapply::pblapply`
             run_lm,              # This function
             cov = cov,           # Covariate matrix, as desribed above
             SCORE = SCORE,       # PRS
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
  
  # Sort results by p-value
  lm.res <- as.data.frame(t(lm.res))
  lm_res.sort <- lm.res[order(lm.res$pval), ]
  return(lm_res.sort)
}
