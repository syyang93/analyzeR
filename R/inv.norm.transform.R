#' Code for an inverse normal transformation, from Ryan
#' 
#' @param mtDNA.resid.15nat The variable that you'd like to inverse normal transform 
#' @export
#' 
#' @return The same value, but after an inverse normal transformation
#' 
#' @examples
#' mtDNA.resid.15nat.inv <- inv.norm.transform(mtDNA.resid.15nat)

inv.norm.transform <- function(mtDNA.resid.15nat){
  mtDNA.resid.15nat=scale(qnorm((rank(mtDNA.resid.15nat, na.last="keep")-0.5) / sum(!is.na(mtDNA.resid.15nat))))
  return(mtDNA.resid.15nat)
}