#' Makes forestplots for genes of interest across tissues
#' 
#' @param gene gene you want to look at
#' @param combined dataframe with genes + effect sizes + standard errors for each tissue
#' @param title what you want to name the plot you make
#' 
#' @export
#' 
#' @return writes out a forestplot to a pdf
#' 
#' @examples
#' make.forestplot('YBX1P10', combined)


make.forestplot <- function(gene, combined, title = ''){
  # you will save the plot as gene_title
  gene.title <- paste0(gene, title)
  
  # required packages:
  require(forestplot)
  require(meta)
  
  # get effects for specific gene:	
  gene.only <- subset(combined, symbol == gene) # 49 tissues --> but might not be same for every tissue!
  
  # get means and upper/lower limits
  means <- gene.only$beta
  upper <- gene.only$beta+gene.only$SE
  lower <- gene.only$beta-gene.only$SE
  
  # random effects meta-analysis, inverse variance weighted
  m1 <- metamean(gene.only$samps, gene.only$beta, gene.only$SE, comb.random = T)
  
  # add meta to vectors
  means <- c(NA, NA, NA, means, NA, m1$TE.random)
  upper <- c(NA, NA, NA, upper, NA, m1$TE.random + m1$seTE.random)
  lower <- c(NA, NA, NA, lower, NA, m1$TE.random - m1$seTE.random)
  
  # create table text
  text <-cbind(c(paste0("Forestplot for ", gene), NA, "Tissue", gene.only$tissue, NA, 'Summary (Random)'), 
               c(NA, NA, 'N', gene.only$samps, NA, NA), 
               c(NA, NA,"Beta", formatC(gene.only$beta, format = "e", digits = 2), NA, formatC(m1$TE.random, format = "e", digits = 2)))
  
  # draw forestplot
  pdf(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_code/Cross.tissue.look/forestplots/', gene.title, '.pdf'), width = 13.3, height = 7.5, onefile = F)
  forestplot(text, means, lower, upper, col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), is.summary=c(TRUE, FALSE, TRUE, rep(FALSE, nrow(gene.only)+1), TRUE))
  dev.off()
}