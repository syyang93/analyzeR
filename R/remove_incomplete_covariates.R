#' Filter out individuals with incomplete phenotypes
#' Phenotypes filtered on: cohort, race, hardy, sex, ischemic time, age
#' 
#' @export
#' 
#' @examples

remove_incomplete_covariates <- function(with.mt.mito) {
  indiv_count <- with.mt.mito[[2]]
  with.mt.mito <- with.mt.mito[[1]]
  covariate_frame <- with.mt.mito[,c(which(colnames(with.mt.mito) == 'COHORT'), which(colnames(with.mt.mito) == 'GENDER'), which(colnames(with.mt.mito) == 'AGE'), which(colnames(with.mt.mito) == 'TRISCHD'), which(colnames(with.mt.mito) == 'DTHHRDY'), which(colnames(with.mt.mito) == 'RACE'), which(colnames(with.mt.mito) == 'submitted_subject_id'))]
  
  length(which(complete.cases(covariate_frame)))
  
  incomplete <- which(!(complete.cases(covariate_frame)))
  omit <- which(with.mt.mito$submitted_subject_id %in% covariate_frame$submitted_subject_id[incomplete])
  print(paste0('incomplete covariates: ', length(omit)))
  if (length(omit) != 0){
    with.mt.complete <- with.mt.mito[-omit,]
  } else{with.mt.complete <- with.mt.mito}
  indiv_count <- rbind(indiv_count, c(length(omit), 'Incomplete covariate information'))
  return(list(with.mt.complete, indiv_count))
}