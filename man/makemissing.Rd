\name{makemissing}
\alias{makemissing}
\title{Create random missing data}
\usage{
makemissing(data, p = 0.2, cols = 1:ncol(data), rowwise = FALSE,
  nafun = function(x) x)
}
\arguments{
  \item{data}{data.frame}

  \item{p}{Fraction of missing data in each column}

  \item{cols}{Which columns (name or index) to alter}

  \item{rowwise}{Should missing occur row-wise (either none
  or all selected columns are missing)}

  \item{nafun}{(Optional) function to be applied on
  data.frame before return (e.g. \code{na.omit} to return
  complete-cases only)}
}
\value{
data.frame
}
\description{
Generates missing entries in data.frame/matrix
}
\author{
Klaus K. Holst
}
\keyword{utilities}

