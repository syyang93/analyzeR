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

get.vars <- function(lm_full, permutevar = 'mtDNA_adjust_AGE', covnames = c(paste0('PC', 1:10), paste0('Genotyping.PC', 1:3), 'sex', 'COHORT'))
{
  require(dplyr)
  tx_expr <- lm_full[, grep("ENSG", colnames(lm_full))]
  cov <- dplyr::select(lm_full, covnames)
  gene.ids <- colnames(tx_expr)
  
  # search for the main variable you're testing associations for
  search.term <- paste0('^', permutevar, '$')
  permute.index <- grep(search.term, colnames(lm_full))
  SCORE <- as.numeric(lm_full[, permute.index])
  all.vars <- list(tx_expr, cov, gene.ids, SCORE)
  return(all.vars)
}
