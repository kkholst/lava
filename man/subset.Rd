\name{subset.lvm}
\alias{subset.lvm}
\alias{measurement}
\title{ Extract subset of latent variable model}
\description{
  Extract measurement models or user-specified subset of model 
}
\usage{
\method{subset}{lvm}(x,vars,...)
measurement(x,silent=TRUE,debug=FALSE)
}
\arguments{
  \item{x}{\code{lvm}-object.}
  \item{vars}{Character vector or formula specifying variables to include in subset.}
  \item{silent}{Logical variable which indicates whether messages are turned on/off.}
  \item{debug}{For debugging only}
  \item{\dots}{}
}
\value{
  A \code{lvm}-object.
}
\examples{
m <- lvm(c(y1,y2)~x1+x2)
subset(m,~y1+x1)
}
\author{Klaus K. Holst}
\keyword{models}
\keyword{regression}
