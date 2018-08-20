#' Adds the rnaseq covariates and makes the MT_count column
#' 
#' @export
#' 

add_cov <- function(with.mt, with_gtex){
  with.cov <- merge(with.mt, with_gtex, by = 'submitted_subject_id')
  with.cov$MT_count <- with.cov$mito_percents
  return(with.cov)
}
