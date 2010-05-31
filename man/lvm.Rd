\name{lvm}
\alias{lvm}
\alias{print.lvm}
\alias{summary.lvm}
\title{Initialize new latent variable model}
\description{
  Function that constructs a new latent variable model object
}
\usage{
lvm(var = NULL, silent=FALSE, ...)
}
\arguments{
  \item{var}{ Vector of variable names. Optional but gives control of
    the sequence of appearance of the variables. The argument can be
    given as a character vector or formula, e.g. \code{~y1+y2} is
    equivalent to \code{c("y1","y2")}. Alternatively the argument can be
    a formula specifying a linear model.
  }
  \item{silent}{Logical variable which indicates whether messages are turned
  on/off.}
  \item{\dots}{}
}
\value{
  Returns an object of class \code{lvm}.
}
\author{Klaus K. Holst}
\seealso{ \code{\link{regression}}, \code{\link{covariance}},
  \code{\link{intercept}}, ... }
\examples{
m <- lvm() # Empty model
m1 <- lvm(y~x) # Simple linear regression
m2 <- lvm(~y1+y2) # Model with two independent variables (argument)
}
\keyword{models}
\keyword{regression}
