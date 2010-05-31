\name{bootstrap}
\alias{bootstrap}
\alias{bootstrap.lvm}
\alias{bootstrap.lvmfit}
\title{Calculate bootstrap estimates of a lvm}
\description{
  Draws samples from a non-parametric bootstrap
}
\usage{
\method{bootstrap}{lvm}(x,R=100,data,control=list(),constraints=TRUE,sd=FALSE,...)

\method{bootstrap}{lvmfit}(x,R=100,data=model.frame(x),control=list(start=coef(x)),...) 
}\arguments{
  \item{x}{\code{lvm}-object.}
  \item{R}{Number of bootstrap samples}
  \item{data}{The data to resample from}
  \item{control}{Options to the optimization routine}
  \item{constraints}{Logical indicating whether non-linear parameter
    constraints should be included in the bootstrap procedure}
  \item{sd}{Logical indicating whether standard error estimates
    should be included in the bootstrap procedure}
  \item{\dots}{Additional arguments, e.g. choice of estimator.}
}
\value{
A \code{bootstrap.lvm} object.
}
\seealso{
  \code{\link{confint.lvmfit}}
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
\author{Klaus K. Holst}
\keyword{models}
\keyword{regression}
