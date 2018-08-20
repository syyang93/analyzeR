#' Filters out people who had incomplete RNAseq downloads from GTEx to recount
#' 
#' @param with_gtex the scaled rse gene counts (along with the GTEx ID info)
#' 
#' @export
#' 
#' @return The data frame, after removing people who had incomplete RNAseq downloads
#' 
#' 

filtering_steps <- function(with_gtex) {
  indiv_count <- with_gtex[[2]]
  with_gtex <- with_gtex[[1]]
  failed_dl <- which(with_gtex$read_count_as_reported_by_sra != with_gtex$reads_downloaded)
  if(length(failed_dl) != 0) {
    print(paste0('number of failed downloads:', length(failed_dl)))
    print(with_gtex$submitted_subject_id[failed_dl])
    with_gtex <- with_gtex[-failed_dl, ]
  }
  indiv_count <- rbind(indiv_count, c(length(failed_dl), 'Failed RNAseq Downloads'))
  return(list(with_gtex, indiv_count))
}