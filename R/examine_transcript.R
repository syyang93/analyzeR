#' Gets the transcripts for people who remain after filtering
#' 
#' @param counts1 data frame containing the scaled counts from rse_gene (assay(rse_gene))
#' @param with_gtex data frame with GTEx IDs
#' 
#' @export
#' 
#' @return Transcripts for those individuals
#' 
#' @examples
#' examine_transcript(blood_full, 'ENSG00000141527.16')

examine_transcript <- function(lm_data, transcript, indiv = 304, correct_for = 'as.factor(lm_data$SEX) + as.numeric(lm_data$AGE) + as.numeric(lm_data$RACE) +
                  lm_data$PC1 + lm_data$PC2 + lm_data$PC3 + lm_data$PC4 + lm_data$PC5 + lm_data$PC6 + lm_data$PC7 + lm_data$PC8 + lm_data$PC9 + lm_data$PC10'){
  require(ggplot2)
  resids<-as.data.frame(matrix(NA, nrow=indiv, ncol=2))
  colnames(resids) <- c('transcript', 'CN')
  index <- which(colnames(lm_data) == transcript)
  form <- as.formula(paste0('lm_data[,index] ~ ', correct_for))
  lm_test <- lm(form, na.action=na.exclude)
  resids$transcript<-scale(resid(lm_test))
  resids$CN<-lm_data$mtDNA_adjust_AGE
  resids <- na.omit(resids)
  print(summary(lm(resids$transcript~resids$CN)))
  g <- ggplot(resids, aes(x=CN, y=transcript))+geom_point()+geom_smooth(method="lm", formula=y~x)
  print(summary(lm(resids$transcript~resids$CN))$coef['resids$CN',])
  return(g)
}

# examine_transcript(lm_full, 'ENSG00000225972.1', 291) + ylab('MT-ND1P23')
# lm_full$submitted_subject_id[order(lm_full$ENSG00000225972.1)]
# write.csv(lm_full$submitted_subject_id[order(lm_full$ENSG00000225972.1, decreasing = T)], file = '/Users/arkinglab/Desktop/vamsee.csv', quote = F, row.names = F)


# examine_transcript(lm_full, 'ENSG00000256525.6', 308) + ggtitle('POLG2')
# examine_transcript(lm_full, 'ENSG00000140521.12', 308) + ggtitle('POLG')
# examine_transcript(lm_full, "ENSG00000108064.10", 308) + ggtitle('TFAM')
# examine_transcript(lm_full, "ENSG00000186814.12", 308) + ggtitle('ZSCAN30')
# examine_transcript(lm_full, "ENSG00000108064.10", 308) + ggtitle('TFAM')


# examine_transcript(lm_full, 'ENSG00000256525.6', 308) + ggtitle('POLG2')
# examine_transcript(lm_full, 'ENSG00000140521.12', 308) + ggtitle('POLG')
# examine_transcript(lm_full, "ENSG00000108064.10", 308) + ggtitle('TFAM')
# examine_transcript(lm_full, "ENSG00000186814.12", 308) + ggtitle('ZSCAN30')
# examine_transcript(lm_full, "ENSG00000178057.14", 291) + ggtitle('NDUFAF3')g
