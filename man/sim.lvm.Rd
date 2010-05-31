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
  \item{cond}{  }
  \item{sigma}{ }
  \item{rho}{  }
  \item{\dots}{  }
}
\details{
  E.g. \eqn{ E(y|x) = 2*x^2 } could be specified as
  
  \code{regression(m, "y3", fn=function(x) x^2) <- "x$2"}
}

\author{ Klaus K. Holst }

\keyword{ models }
\keyword{ datagen }
\keyword{ regression }
