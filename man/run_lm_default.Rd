% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_lm_default.R
\name{run_lm_default}
\alias{run_lm_default}
\title{Basic function that is used in run.all.lms.  Needs to be called by run_lm.  Adapted from Vamsee (github.com/vkp3/pillalamarRi)}
\usage{
run_lm_default(expr, cov, SCORE, omit.outlier = T, outlier_sd = 3)
}
\arguments{
\item{expr}{Expression vector (numeric) with length = #samples}

\item{cov}{Regression covariates in form [cov x samples]}

\item{SCORE}{The PRS, not included in `cov`}

\item{omit.outlier}{Whether or not you want to omit gene expression outliers}

\item{method}{Choose between 'default' or 'two-stage' for lm() method (see desc. in support functions below)}
}
\value{
A [1 x 8] vector output from an lm() like below:
        ['intercept', 'beta', 'SE', 't_value', 'pval', 'beta.conf.low', 'beta.conf.high', 'corr.rho']

Author: Vamsee Pillalamarri
}
\description{
formula: Expr ~ SCORE + Cov, where Cov = covariates
}
\details{
Runs linear regression to estimate effect of SCORE on gene expression
}
