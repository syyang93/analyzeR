#' Function that tests for associations and returns results (TRANSCRIPT LEVEL).  DOES NOT merge with gene info, for use with tpm datasets.
#' This function will do a single linear regression and extract the coefficients from the regression.
#' 
#' @param lm_data Our dependent variables.  Column = dependent variable, rows = observations (in this case, individuals) (numeric vector)
#' @param lm_full The data frame containing genes, covariates, and mtDNA-CN info
#' @param to_correct A string containing additional covariates to correct for.  ex: '+as.factor(covariates$DTHHRDY)+as.factor(covariates$RACE)' 
#' 
#' @export
#' 
#' @return Coefficients from the regression
#' 
#' @examples
#' lm_results <- testing_assoc(lm_full)
 
testing_assoc.transcript <- function (lm_full, to_correct = "+ as.factor(covariates$SEX) + as.numeric(covariates$AGE)+ as.numeric(covariates$RACE) + covariates$PC1 + covariates$PC2 +covariates$PC3+ covariates$PC4 + covariates$PC5 + covariates$PC6 + covariates$PC7 + covariates$PC8 + covariates$PC9 + covariates$PC10")
{
  for_test <- lm_full[, grep("ENST", colnames(lm_full))]
  covariates <- lm_full[, -grep("ENST", colnames(lm_full))]
  MT_count_index <- grep("mtDNA_adjust_AGE", colnames(lm_full))
  MT_count <- as.numeric(lm_full[, MT_count_index])
  lm_results <- apply(for_test, 2, lm_test, count = MT_count,
                      covariates = covariates, correct_for = to_correct)
  lm_results <- as.data.frame(t(lm_results))
  lm_results <- lm_results[order(lm_results$`Pr(>|t|)`), ]
  return(lm_results)

  # deprecated, for v.6 with recount
  
  # load('R_objects/info.rds')
  # lm_results$gene_id <- rownames(lm_results)
  # with_genes <- merge(info, lm_results, by = 'gene_id')
  # with_genes <- as.data.frame(with_genes)
  # with_genes <- with_genes[order(with_genes$`Pr(>|t|)`),]
  # return(with_genes)
}