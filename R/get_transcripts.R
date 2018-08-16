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

get_transcripts <- function(counts1, with_gtex) {
  tissue_transcripts <- pituitary_transcripts <- counts1[,colnames(counts1) %in% with_gtex$run]
  tissue_transcripts_t <- as.data.frame(t(tissue_transcripts))
  tissue_transcripts_t$submitted_subject_id <- with_gtex$submitted_subject_id
  rownames(tissue_transcripts_t) <- tissue_transcripts_t$submitted_subject_id
  return(tissue_transcripts_t)
}
