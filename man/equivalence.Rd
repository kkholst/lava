\name{equivalence}
\alias{equivalence}
\title{Identify candidates of equivalent models}
\usage{
  equivalence(x, rel, tol = 0.001, k = 1, omitrel = TRUE,
    ...)
}
\arguments{
  \item{x}{\code{lvmfit}-object}

  \item{rel}{Formula or character-vector specifying two
  variables to omit from the model and subsequently search
  for possible equivalent models}

  \item{tol}{Define two models as empirical equivalent if
  the absolute difference in score test is less than
  \code{tol}}

  \item{k}{Number of parameters to test simultaneously. For
  \code{equivalence} the number of additional associations
  to be added instead of \code{rel}.}

  \item{omitrel}{if \code{k} greater than 1, this boolean
  defines wether to omit candidates containing \code{rel}
  from the output}

  \item{\dots}{Additional arguments to be passed to the low
  level functions}
}
\description{
  Identifies candidates of equivalent models
}
\author{
  Klaus K. Holst
}
\seealso{
  \code{\link{compare}}, \code{\link{modelsearch}}
}

