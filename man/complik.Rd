% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complik.R
\name{complik}
\alias{complik}
\title{Composite Likelihood for probit latent variable models}
\usage{
complik(
  x,
  data,
  k = 2,
  type = c("all", "nearest"),
  pairlist,
  messages = 0,
  estimator = "normal",
  quick = FALSE,
  ...
)
}
\arguments{
\item{x}{\code{lvm}-object}

\item{data}{data.frame}

\item{k}{Size of composite groups}

\item{type}{Determines number of groups. With \code{type="nearest"} (default)
only neighboring items will be grouped, e.g. for \code{k=2}
(y1,y2),(y2,y3),... With \code{type="all"} all combinations of size \code{k}
are included}

\item{pairlist}{A list of indices specifying the composite groups. Optional
argument which overrides \code{k} and \code{type} but gives complete
flexibility in the specification of the composite likelihood}

\item{messages}{Control amount of messages printed}

\item{estimator}{Model (pseudo-likelihood) to use for the pairs/groups}

\item{quick}{If TRUE the parameter estimates are calculated but all additional
information such as standard errors are skipped}

\item{\dots}{Additional arguments parsed on to lower-level functions}
}
\value{
An object of class \code{estimate.complik} inheriting methods from \code{lvm}
}
\description{
Estimate parameters in a probit latent variable model via a composite
likelihood decomposition.
}
\examples{
m <- lvm(c(y1,y2,y3)~b*x+1*u[0],latent=~u)
ordinal(m,K=2) <- ~y1+y2+y3
d <- sim(m,50,seed=1)
if (requireNamespace("mets", quietly=TRUE)) {
   e1 <- complik(m,d,control=list(trace=1),type="all")
}
}
\seealso{
estimate
}
\author{
Klaus K. Holst
}
\keyword{models}
\keyword{regression}
