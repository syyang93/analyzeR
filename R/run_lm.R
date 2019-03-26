#' A wrapper for deciding either to run a default lm or a two-stage lm.  I always use the default.  Adapted from Vamsee (github.com/vkp3/pillalamarRi)
#' 
#' Runs linear regression to estimate effect of SCORE on gene expression
#' @description formula: Expr ~ SCORE + Cov, where Cov = covariates
#' 
#' @param expr Expression vector (numeric) with length = #samples
#' @param cov Regression covariates in form [cov x samples]
#' @param SCORE The PRS, not included in `cov`
#' @param omit.outlier Whether or not you want to omit gene expression outliers
#' @param method Choose between 'default' or 'two-stage' for lm() method (see desc. in support functions below)
#' @return A [1 x 8] vector output from an lm() like below:
#'         ['intercept', 'beta', 'SE', 't_value', 'pval', 'beta.conf.low', 'beta.conf.high', 'corr.rho']
#' 
#' Author: Vamsee Pillalamarri
#' @export
 

run_lm <- function(expr, cov, SCORE, omit.outlier = T, method='default') {
  res <- switch (method,
                 "default" = run_lm_default(expr, cov, SCORE, omit.outlier),
                 "two-stage" = run_lm_two_stage(expr, cov, SCORE, omit.outlier)
  )
  return(res)
}
