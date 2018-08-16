#' Adds GTEx ID information to RNAseq data from recount
#' 
#' @param scaled_rse the scaled rse gene counts after it's been subsetted for tissue detail.
#' 
#' @export
#' 
#' @return Returns the scaled rse gene counts with GTEx ID information
#' 
#' @examples
#' pituitary <- subset(as.data.frame(colData(scaled)), smtsd== "Pituitary" & smafrze == 'USE ME')
#' pituitary <- getid(pituitary)


getid <- function(scaled_rse){
  scaled_rse$submitted_subject_id2 <- substr(scaled_rse$sampid, 6, 10)
  scaled_rse$submitted_subject_id2 <- gsub('-', '', scaled_rse$submitted_subject_id2)
  scaled_rse$submitted_subject_id <- paste0('GTEX-', scaled_rse$submitted_subject_id2)
  scaled_rse$submitted_subject_id2 <- NULL
  return(scaled_rse)
}