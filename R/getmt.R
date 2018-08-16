#' Load and merge mtDNA-CN info
#' 
#' @param transcripts_t scaled rse gene counts 
#' @param mt_frame the data frame containing mtDNA-CN and read depth information
#' 
#' @export
#' 
#' @return Returns the data frame with mtDNA-CN info
#' 
#' @examples



getmt <- function(transcripts_t, mt_frame = mtDNA_CN_with_pheno){
  index <- which(colnames(transcripts_t) == 'submitted_subject_id')
  only_id <- transcripts_t[,(index-1):index]
  with.mt <- merge(only_id, mt_frame, by = 'submitted_subject_id')
  return(with.mt)
}
