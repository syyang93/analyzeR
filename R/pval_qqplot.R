#' Function for drawing a p-value qqplot
#' 
#' @param lm_object The results of an lm that you want a p-value QQ-plot for
#' 
#' @export
#' 
#' @return A p-value QQ-plot
#' 
#' @examples
#' pval_qqplot(subcut_lm_results)


pval_qqplot <- function(lm_object, title = 'QQ-plot') {
  pvals <- lm_object$`Pr(>|t|)`
  observed <- sort(pvals)
  observed2 <- c(length(pvals))
  observed2null <- -(log10(observed2 / (length(observed2)+1)))
  pvals <- c(pvals, observed2null)
  observed <- sort(pvals)
  lobs <- -(log10(observed))
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  #creating uniform disn
  m <- title
  plot(c(0,20), c(0,20), col = 'red', lwd = 4, type = 'l', xlab = 'Expected (-logP)', ylab= 'Observed (-logP)', xlim = c(0,7), ylim = c(0,7), las = 1, xaxs = 'i', yaxs = 'i', bty = 'l', main = m)
  points(lexp, lobs, pch=23, cex = 0.5, col = 'black', bg = 'black')
}
