#' Gets the adjusted variable plot for mtDNA vs transcripts for people who remain after filtering
#' 
#' @param lm_data lm_full you want to examine the transcript in 
#' @param symbol a string with the gene symbol you want to look at
#' @param gene_key  a gene key with symbol -> transcript information (tpms_gene_key.rds)
#' @param correct_for a string that says what you need to correct transcripts for
#' @param omit.outlier T/F whether you want to omit outliers > 3SD away from the mean
#' @param col A string, what you want to color the adjusted variable plot with
#' 
#' @export
#' 
#' @return Transcripts for those individuals
#' 
#' @examples
#' examine_symbol(blood_full, 'POLG', omit.outlier = T, correct_for = 'as.factor(lm_data$SEX)', col = 'esoph_day')

examine_symbol <- function (lm_data, symbol, gene_key, omit.outlier = T, correct_for = "as.factor(lm_data$SEX) + as.numeric(lm_data$AGE) + as.numeric(lm_data$RACE) +\n                  lm_data$PC1 + lm_data$PC2 + lm_data$PC3 + lm_data$PC4 + lm_data$PC5 + lm_data$PC6 + lm_data$PC7 + lm_data$PC8 + lm_data$PC9 + lm_data$PC10 + lm_data$Genotyping.PC1 + lm_data$Genotyping.PC2 + lm_data$Genotyping.PC3",
                            col = NA)
{
  indiv <- nrow(lm_data)
  require(ggplot2)
  index <- grep(paste0("^", symbol, "$"), gene_key$symbol)
  if (length(index) == 0) {
    print("Symbol not found in dataframe")
  }
  transcript <- gene_key$gene_id[index]
  resids <- as.data.frame(matrix(NA, nrow = indiv, ncol = 2))
  colnames(resids) <- c("transcript", "CN")
  index <- which(colnames(lm_data) == transcript)
  if (omit.outlier == T) {
    outlier_sd <- 3
    m <- mean(lm_data[, index])
    s <- sd(lm_data[, index])
    outliers <- which(lm_data[, index] > m + outlier_sd *
                        s | lm_data[, index] < m - outlier_sd * s)
    if (length(outliers) != 0) {
      lm_data[outliers, index] <- NA
      # add a second iteration to remove more outliers?
      ## iterate!
      
      
    }
  }
  #### added this line to scale things!

  #### removed scaling.
  form <- as.formula(paste0("lm_data[,index] ~ lm_data$mtDNA_adjust_AGE + ",
                            correct_for))
  lm_test <- lm(form, na.action = na.exclude)
  
  print(coef(summary(lm_test))[2, ])
  form2 <- as.formula(paste0("lm_data[,index] ~ ", correct_for))
  lm_test <- lm(form2, na.action = na.exclude)
  
  # don't scale residuals
  resids$transcript <- resid(lm_test)
  resids$CN <- as.numeric(lm_data$mtDNA_adjust_AGE)
  resids$col <- lm_data[[col]]
  resids$esoph_day <- lm_data$esoph_day
  if (is.na(col)) {
    g <- ggplot(resids, aes(x = CN, y = transcript)) + geom_point() +
      geom_smooth(method = "lm", formula = y ~ x) + xlab("mtDNA-CN") +
      ylab(paste0("residuals for ", symbol))
  }
  else {
    g <- ggplot(resids, aes(x = CN, y = transcript, col = col)) +
      geom_point() + geom_smooth(method = "lm", formula = y ~
                                   x) + xlab("mtDNA-CN") + ylab(paste0("residuals for ",
                                                                       symbol))
  }
  return(g)
}
