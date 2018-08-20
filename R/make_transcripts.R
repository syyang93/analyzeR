#' Filter the transcripts to match the individuals that have been filtered out
#' 
#' @export
#' 

make_transcripts <- function(transcripts_t, with.mt.complete){
  kept <- which(transcripts_t$submitted_subject_id %in% with.mt.complete$submitted_subject_id)
  transcripts_t <- transcripts_t[kept,]
  return(transcripts_t)
}