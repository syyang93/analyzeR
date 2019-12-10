#' Basic function that is used in run.all.lms.  Needs to be called by run_lm.  Adapted from Vamsee (github.com/vkp3/pillalamarRi)
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


run_lm_default <- function(expr, cov, SCORE, omit.outlier = T, outlier_sd = 3) {
  expr <- as.numeric(expr)
  
  # stop scaling the expression values!
  # expr <- scale(expr)
  expr_cov <- cbind(SCORE, expr, cov)
  
  if(omit.outlier == T)
  {
    m <- mean(expr)
    s <- sd(expr)
    outliers <- which(expr > m + outlier_sd*s | expr < m - outlier_sd*s)
    if(length(outliers) != 0){
      expr_cov$expr[outliers] <- NA
      # second iteration
      m <- mean(expr)
      s <- sd(expr)
      outliers <- which(expr > m + outlier_sd*s | expr < m - outlier_sd*s)
      if(length(outliers) != 0){
        # remove more outliers?  if outlier removal has to be done twice, omit that gene from testing?  
        
        
      }
    }
  }
  
  # # Run lm() normal procedure
  lm.fit <- lm(expr ~ ., data = expr_cov)
  lm.fit.summary <- summary(lm.fit)
  
  # Get corr
  cor_expr_score <-
    cor(expr, SCORE)
  
  # Capture p-val, etc.
  lm.res_ <-
    as.data.frame(t(coef(lm.fit.summary)['SCORE',]))
  intercept <- coef(lm.fit)[1]
  names(intercept) <- NULL
  lm.res_ <-
    cbind(data.frame(intercept),
          lm.res_,
          t(confint(lm.fit)['SCORE', ]),
          cor_expr_score)
  
  return(as.matrix(lm.res_))
}