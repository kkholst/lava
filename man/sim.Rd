\name{sim}
\alias{binomial.lvm}
\alias{distribution}
\alias{distribution<-}
\alias{distribution<-.lvm}
\alias{distribution.lvm}
\alias{functional}
\alias{functional<-}
\alias{functional<-.lvm}
\alias{functional.lvm}
\alias{heavytail}
\alias{heavytail<-}
\alias{logit.lvm}
\alias{normal.lvm}
\alias{poisson.lvm}
\alias{probit.lvm}
\alias{sim}
\alias{sim.lvm}
\alias{sim.lvmfit}
\alias{student.lvm}
\alias{uniform.lvm}
\alias{weibull.lvm}
\title{Simulate model}
\usage{
  \method{sim}{lvm}(x, n = 100, p = NULL, normal = FALSE,
    cond = FALSE, sigma = 1, rho = 0.5, X, unlink=FALSE,
    ...)
}
\arguments{
  \item{x}{Model object}

  \item{n}{Number of simulated values/individuals}

  \item{p}{Parameter value (optional)}

  \item{normal}{Logical indicating whether to simulate data
  from a multivariate normal distribution conditional on
  exogenous variables hence ignoring
  functional/distribution definition}

  \item{cond}{for internal use}

  \item{sigma}{Default residual variance (1)}

  \item{rho}{Default covariance parameter (0.5)}

  \item{X}{Optional matrix of covariates}

  \item{unlink}{Return Inverse link transformed data}

  \item{\dots}{Additional arguments to be passed to the low
  level functions}
}
\description{
  Simulate data from a general SEM model including
  non-linear effects and general link and distribution of
  variables.
}
\details{
  E.g. \eqn{ E(y|x) = 2*x^2 } could be specified as

  \code{regression(m, "y3", fn=function(x) x^2) <- "x$2"}
}
\author{
  Klaus K. Holst
}
\keyword{datagen}
\keyword{models}
\keyword{regression}

