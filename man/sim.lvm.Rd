\name{sim}
\alias{sim}
\alias{sim.lvmfit}
\alias{sim.lvm}
\alias{functional}
\alias{functional}
\alias{functional<-}
\alias{functional.lvm}
\alias{functional<-.lvm}
\alias{distribution}
\alias{distribution}
\alias{distribution<-}
\alias{distribution.lvm}
\alias{distribution<-.lvm}
\alias{heavytail}
\alias{heavytail<-}
\alias{weibull.lvm}
\alias{binomial.lvm}
\alias{poisson.lvm}
\alias{uniform.lvm}
\alias{normal.lvm}
\alias{probit.lvm}
\alias{logit.lvm}

\title{ Simulate model }
\description{
  Simulate data from a general SEM model including non-linear effects
  and general link and distribution of variables.
}
\usage{
\method{sim}{lvm}(x, n = 100, p = NULL, normal = FALSE, cond = FALSE, sigma = 1, rho = 0.5, ...)
}
\arguments{
  \item{x}{ Model object }
  \item{n}{ Number of simulated values/individuals  }
  \item{p}{ Parameter value (optional) }
  \item{normal}{ Logical indicating whether to simulate data from a
    multivariate normal distribution conditional on exogenous variables
    hence ignoring functional/distribution definition }
  \item{cond}{for internal use}
  \item{sigma}{Default residual variance (1)}
  \item{rho}{Default covariance parameter (0.5)}
  \item{\dots}{Additional arguments to be passed to the low level functions  }
}
\details{
  E.g. \eqn{ E(y|x) = 2*x^2 } could be specified as
  
  \code{regression(m, "y3", fn=function(x) x^2) <- "x$2"}
}

\author{ Klaus K. Holst }

\keyword{ models }
\keyword{ datagen }
\keyword{ regression }
