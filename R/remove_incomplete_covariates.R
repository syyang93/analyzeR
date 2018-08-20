# take out incomplete covariates

remove_incomplete_covariates <- function(blood.with.mt.mito) {
  covariate_frame <- blood.with.mt.mito[,c(which(colnames(blood.with.mt.mito) == 'COHORT'), which(colnames(blood.with.mt.mito) == 'GENDER'), which(colnames(blood.with.mt.mito) == 'AGE'), which(colnames(blood.with.mt.mito) == 'TRISCHD'), which(colnames(blood.with.mt.mito) == 'DTHHRDY'), which(colnames(blood.with.mt.mito) == 'RACE'), which(colnames(blood.with.mt.mito) == 'submitted_subject_id'))]
  
  length(which(complete.cases(covariate_frame)))
  
  incomplete <- which(!(complete.cases(covariate_frame)))
  omit <- which(blood.with.mt.mito$submitted_subject_id %in% covariate_frame$submitted_subject_id[incomplete])
  print(paste0('incomplete covariates: ', length(omit)))
  if (length(omit) != 0){
    blood.with.mt.complete <- blood.with.mt.mito[-omit,]
  } else{blood.with.mt.complete <- blood.with.mt.mito}
  return(blood.with.mt.complete)
}
