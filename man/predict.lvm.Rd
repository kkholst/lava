\name{predict.lvm}
\alias{predict.lvm}
\alias{predict.lvmfit}
\title{Prediction in structural equation models}
\usage{
  \method{predict}{lvm} (object, x = NULL,
    residual = FALSE, p, data, path = FALSE,
    quick = is.null(x) & !(residual | path), ...)
}
\arguments{
  \item{object}{Model object}

  \item{x}{optional list of (endogenous) variables to
  condition on}

  \item{residual}{If true the residuals are predicted}

  \item{p}{Parameter vector}

  \item{data}{Data to use in prediction}

  \item{path}{Path prediction}

  \item{quick}{If TRUE the conditional mean and variance
  given covariates are returned (and all other calculations
  skipped)}

  \item{\dots}{Additional arguments to lower level
  function}
}
\description{
  Prediction in structural equation models
}
\examples{
m <- lvm()
}

