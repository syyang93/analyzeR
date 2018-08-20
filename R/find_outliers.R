#' Finds PCA outliers - greater than 4 SD from the mean of the first 10 PCs.  
#' 
#' @export
#' 

find_outliers <- function(pca) {
  outlier<- c()
  all_outliers<- c()
  for(i in 1:10){
    a<-subset(rownames(pca$x), pca$x[,i] > (mean(pca$x[,i])+4*sd(pca$x[,i])))
    b<-subset(rownames(pca$x), pca$x[,i] < (mean(pca$x[,i])-4*sd(pca$x[,i])))
    out<-c(a,b)
    outlier <- c(outlier, out)
    print(paste("outliers in PCA",i,":",sep=""))
    print(outlier)
    all_outliers=c(all_outliers, outlier)
    outlier=c()
  }
  all_outliers
}
