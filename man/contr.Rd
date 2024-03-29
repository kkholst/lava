% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contr.R
\name{contr}
\alias{contr}
\alias{parsedesign}
\alias{pairwise.diff}
\title{Create contrast matrix}
\usage{
contr(p, n, diff = TRUE, ...)
}
\arguments{
\item{p}{index of non-zero entries (see example)}

\item{n}{Total number of parameters (if omitted the max number in p will be used)}

\item{diff}{If FALSE all non-zero entries are +1, otherwise the second non-zero element in each row will be -1.}

\item{...}{Additional arguments to lower level functions}
}
\description{
Create contrast matrix typically for use with 'estimate' (Wald tests).
}
\examples{
contr(2,n=5)
contr(as.list(2:4),n=5)
contr(list(1,2,4),n=5)
contr(c(2,3,4),n=5)
contr(list(c(1,3),c(2,4)),n=5)
contr(list(c(1,3),c(2,4),5))

parsedesign(c("aa","b","c"),"?","?",diff=c(FALSE,TRUE))

## All pairs comparisons:
pdiff <- function(n) lava::contr(lapply(seq(n-1), \(x) seq(x, n)))
pdiff(4)
}
