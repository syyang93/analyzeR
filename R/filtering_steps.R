#' Filters out people who had incomplete RNAseq downloads from GTEx to recount
#' 
#' @param with_gtex the scaled rse gene counts (along with the GTEx ID info)
#' 
#' @export
#' 
#' @return The data frame, after removing people who had incomplete RNAseq downloads
#' 
#' @examples
#' 

filtering_steps <- function(with_gtex) {
  failed_dl <- which(with_gtex$read_count_as_reported_by_sra != with_gtex$reads_downloaded)
  if(length(failed_dl) != 0) {
    print(paste0('number of failed downloads:', length(failed_dl)))
    print(with_gtex$submitted_subject_id[failed_dl])
    with_gtex <- with_gtex[-failed_dl, ]
  } else {print('no failed downloads')}
  return(with_gtex)
}
