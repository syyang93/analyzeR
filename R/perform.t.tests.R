#' Will perform t-tests and return pvals, confints, and estimates for difference in means
#' 
#' @param gene.set gene sets you would like to query (read using read.gmt())
#' @param with.gene dataframe containing genes and t-values for association with mtdna-cn
#' 
#' @export
#' 
#' @return a dataframe containing pvals/confints/estimates for gene sets
#' 
#' @examples
#' kegg.gene.sets <- perform.t.tests(kegg.sets, with.gene)


perform.t.tests <- function(gene.sets, with.gene){
  signif.gene.sets <- as.data.frame(matrix(nrow = 1, ncol = 7))
  colnames(signif.gene.sets) <- c('Gene.Set.Name', 'T.test.pval', 'Rank.t.test.pval', 'Num.genes.in.set', 'Beta', 'Confint.Upper', 'Confint.Lower')
  
  for(i in 1: length(gene.sets)){
    test.set <- gene.sets[[i]]
    set.name <- names(gene.sets[i])
    selected.indices <- which(with.gene$symbol %in% test.set)
    in.set <- with.gene[selected.indices,]
    out.set <- with.gene[-selected.indices,]
    if(nrow(in.set) < 2){
      print(paste0("Not enough samples for ", set.name))
    } else{
      t.stats <- t.test(in.set$abs.tscore, out.set$abs.tscore)
      if(t.stats$p.value < (0.05/length(gene.sets))){
        print(set.name)
        print(t.stats$p.value)
        ranked.t.stats <- t.test(in.set$ranked.tvals, out.set$ranked.tvals)
        print(ranked.t.stats$p.value)
        sig <- c(set.name, t.stats$p.value, ranked.t.stats$p.value, nrow(in.set), t.stats$estimate[1]-t.stats$estimate[2], t.stats$conf.int[1], t.stats$conf.int[2])
        signif.gene.sets <- rbind(signif.gene.sets, sig)
      }
    }
  }
  
  
  signif.gene.sets <- na.omit(signif.gene.sets)
  signif.gene.sets$Rank.t.test.pval <- as.numeric(signif.gene.sets$Rank.t.test.pval)
  
  signif.gene.sets <- signif.gene.sets[order(signif.gene.sets$Rank.t.test.pval, decreasing = F),]
  return(signif.gene.sets)
}
