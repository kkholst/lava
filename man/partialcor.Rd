\name{partialcor}
\alias{partialcor}
\title{Calculate partial correlations}
\usage{
partialcor(formula, data, level = 0.95)
}
\arguments{
  \item{formula}{formula speciying the covariates and
  optionally the outcomes to calculate partial correlation
  for}

  \item{data}{data.frame}

  \item{level}{Level of confidence limits}
}
\value{
A coefficient matrix
}
\description{
Calculate partial correlation coefficients and confidence
limits via Fishers z-transform
}
\examples{
m <- lvm(c(y1,y2,y3)~x1+x2)
covariance(m) <- c(y1,y2,y3)~y1+y2+y3
d <- sim(m,500)
partialcor(~x1+x2,d)
}
\author{
Klaus K. Holst
}
\keyword{models}
\keyword{regression}

