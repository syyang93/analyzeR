% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perform.t.tests.R
\name{perform.t.tests}
\alias{perform.t.tests}
\title{Will perform t-tests and return pvals, confints, and estimates for difference in means}
\usage{
perform.t.tests(gene.sets, with.gene)
}
\arguments{
\item{with.gene}{dataframe containing genes and t-values for association with mtdna-cn}

\item{gene.set}{gene sets you would like to query (read using read.gmt())}
}
\value{
a dataframe containing pvals/confints/estimates for gene sets
}
\description{
Will perform t-tests and return pvals, confints, and estimates for difference in means
}
\examples{
kegg.gene.sets <- perform.t.tests(kegg.sets, with.gene)
}
