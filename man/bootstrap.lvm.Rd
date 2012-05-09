\name{bootstrap.lvm}
\alias{bootstrap.lvm}
\alias{bootstrap.lvmfit}
\title{Calculate bootstrap estimates of a lvm object}
\usage{
  \method{bootstrap}{lvm}(x,R=100,data,fun=NULL,control=list(),
    p, parametric=FALSE,
    constraints=TRUE,sd=FALSE,silent=FALSE,...)

  \method{bootstrap}{lvmfit}(x,R=100,data=model.frame(x),
    control=list(start=coef(x)), p=coef(x),
    parametric=FALSE,
    estimator=x$estimator,weight=Weight(x),...)
}
\arguments{
  \item{x}{\code{lvm}-object.}

  \item{R}{Number of bootstrap samples}

  \item{fun}{Optional function of the (bootstrapped)
  model-fit defining the statistic of interest}

  \item{data}{The data to resample from}

  \item{control}{Options to the optimization routine}

  \item{p}{Parameter vector of the null model for the
  parametric bootstrap}

  \item{parametric}{If TRUE a parametric bootstrap is
  calculated. If FALSE a non-parametric (row-sampling)
  bootstrap is computed.}

  \item{constraints}{Logical indicating whether non-linear
  parameter constraints should be included in the bootstrap
  procedure}

  \item{estimator}{String definining estimator, e.g.
  'gaussian' (see \code{estimator})}

  \item{weight}{Optional weight matrix used by
  \code{estimator}}

  \item{sd}{Logical indicating whether standard error
  estimates should be included in the bootstrap procedure}

  \item{silent}{Suppress messages}

  \item{\dots}{Additional arguments, e.g. choice of
  estimator.}
}
\value{
  A \code{bootstrap.lvm} object.
}
\description{
  Draws non-parametric bootstrap samples
}
\examples{
m <- lvm(y~x)
d <- sim(m,100)
e <- estimate(y~x, d)
\dontrun{
B <- bootstrap(e,R=100)
B
}
}
\author{
  Klaus K. Holst
}
\seealso{
  \code{\link{confint.lvmfit}}
}
\keyword{models}
\keyword{regression}

