#' Function that tests for associations and returns results
#' This function will do a single linear regression and extract the coefficients from the regression.
#' 
#' @param lm_data Our dependent variables.  Column = dependent variable, rows = observations (in this case, individuals) (numeric vector)
#' @param count Our independent variable (numeric vector)
#' @param covariates Covariates associated with each individual. Columns = covariates, rows = observations. (data.frame)
#' @param to_correct A string containing additional covariates to correct for.  ex: '+as.factor(covariates$DTHHRDY)+as.factor(covariates$RACE)' 
#' 
#' @export
#' 
#' @return Coefficients from the regression
#' 
#' @examples
#' lm_results <- testing_assoc(lm_full)

testing_assoc <- function(lm_full){
  for_test <- lm_full[,1:length(grep('ENSG', colnames(lm_full)))] # the transcripts you're testing as your dependent variables
  covariates <- lm_full[,(length(grep('ENSG', colnames(lm_full)))+1):ncol(lm_full)] # the covariates
  MT_count_index <- grep('mtDNA_adjust_AGE', colnames(lm_full))
  MT_count <- as.numeric(lm_full[,MT_count_index]) # the MT_count
  
  to_correct <- '+ as.factor(covariates$GENDER) + as.numeric(covariates$AGE)+ as.numeric(covariates$smrin) + as.factor(covariates$smcenter) + as.factor(covariates$RACE)+ as.numeric(covariates$TRISCHD) +as.factor(covariates$COHORT) + as.factor(covariates$DTHHRDY) + covariates$PC1 + covariates$PC2 +covariates$PC3+ covariates$PC4 + covariates$PC5'
  lm_results <- apply(for_test, 2, lm_test, count = MT_count, covariates = covariates, correct_for = to_correct)
  lm_results <- as.data.frame(t(lm_results))
  lm_results <- lm_results[order(lm_results$`Pr(>|t|)`),]
  load('R_objects/info.rds')
  lm_results$gene_id <- rownames(lm_results)
  with_genes <- merge(info, lm_results, by = 'gene_id')
  with_genes <- as.data.frame(with_genes)
  with_genes <- with_genes[order(with_genes$`Pr(>|t|)`),]
  return(with_genes)
}
