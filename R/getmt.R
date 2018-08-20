#' Load and merge mtDNA-CN info
#' 
#' @param transcripts_t scaled rse gene counts 
#' @param mt_frame the data frame containing mtDNA-CN and read depth information
#' 
#' @export
#' 
#' @return Returns a data frame with mtDNA-CN info
#' 
#' @examples

getmt <- function(transcripts_t, indiv_count, mt_frame = mtDNA_CN_with_pheno){
  index <- which(colnames(transcripts_t) == 'submitted_subject_id')
  only_id <- transcripts_t[,index]
  with.mt <- mt_frame[which(mt_frame$submitted_subject_id %in% only_id),]
  no_mt_info <- only_id[-which(only_id %in% mt_frame$submitted_subject_id)]
  indiv_count <- rbind(indiv_count, c(length(no_mt_info), 'No mtDNA-CN info'))
  return(list(with.mt, indiv_count))
}
