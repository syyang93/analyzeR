#' Function that will return GTEx subject IDs when given sample IDs
#' 
#' @param x GTEx sample IDs
#' 
#' @export
#' 
#' @return GTEx subject-level IDs
#' 
#' @examples
#' blood.runs$submitted_subject_id <- make.subjids(blood.runs$sampid)

make.subjids <- function(sampids){
  subjids <- substr(sampids, 6, 10)
  subjids <- gsub('-', '', subjids)
  subjids <- paste0('GTEX-', subjids)
  return(subjids)
}