\name{estimate}
\alias{estimate}
\alias{estimate.lvm}
\alias{estimate.list}
\alias{print.lvmfit}
\alias{summary.lvmfit}
\title{Estimation of parameters in a Latent Variable Model (lvm)}
\description{
  Estimate parameters. MLE, IV or user-defined estimator.
}
\usage{
\method{estimate}{lvm}(x, data, estimator = "gaussian", control = list(),
                 missing = FALSE, weight, fix, index = TRUE, graph =
                 FALSE, debug = FALSE, silent = lava.options()$silent,
                 quick = FALSE, ...)
}
\arguments{
  \item{x}{ \code{lvm}-object }
  \item{data}{ \code{data.frame} }
  \item{estimator}{ String defining the estimator (see details below) }
  \item{control}{ control/optimization parameters (see details below)}
  \item{weight}{ Optional weights to used by the chosen estimator. }
  \item{missing}{ Logical variable indiciating how to treat missing
    data. Setting to FALSE leads to complete case analysis. In the
    other case likelihood based inference is obtained by integrating out
  the missing data under assumption the assumption that data is missing
  at random (MAR). }
  \item{index}{ For internal use only }
  \item{graph}{ For internal use only }
  \item{fix}{ Logical variable indicating whether parameter restriction
    automatically should be imposed (e.g. intercepts of latent variables
    set to 0 and at least one regression parameter of each measurement
    model fixed to ensure identifiability.)}
  \item{debug}{ for debugging only }
  \item{quick}{ If TRUE the parameter estimates are calculated but all
    additional information such as standard errors are skipped }
  \item{silent}{ Logical argument indicating whether information should
    be printed during estimation }
  \item{\dots}{Additional arguments to be passed to the low level functions}
}
\value{
  A \code{lvmfit}-object.
}
\details{
A list of parameters controlling the estimation and optimization
procedures is parsed via the \code{control} argument.
By default Maximum Likelihood is used assuming
multivariate normal distributed measurement errors. A list with one or
more of the following elements is expected:
\describe{
  \item{start:}{Starting value. The order of the parameters can
    be shown by calling \code{coef} (with \code{mean=TRUE}) on the
    \code{lvm}-object or with \code{plot(..., labels=TRUE)}. Note that
    this requires a check that it is actual the model being estimated,
    as \code{estimate} might add additional restriction to the model,
    e.g. through the \code{fix} and \code{exo.fix} arguments. The
    \code{lvm}-object of a fitted model can be extracted with the
    \code{Model}-function.}
  \item{starterfun:}{Starter-function with syntax
    \code{function(lvm, S, mu)}.  Three builtin functions are available:
    \code{startvalues}, \code{startvalues2},
    \code{startvalues3}.}
  \item{estimator:}{ String defining which estimator to use (Defaults to ``\code{gaussian}'')}
    \item{meanstructure}{Logical variable indicating whether to fit
    model with meanstructure.}
  \item{method:}{ String pointing to alternative optimizer
    (e.g. \code{optim} to use simulated annealing).}
  \item{control:}{ Parameters passed to the optimizer (default \code{stats::nlminb}).}

  \item{tol:}{ Tolerance of optimization constraints on lower
    limit of variance parameters.  }
}}

\author{Klaus K. Holst}

\seealso{ \code{\link{score}}, \code{\link{information}}, ... }
\examples{
m <- lvm()
regression(m) <- y~v1+v2
regression(m) <- c(v1,v2) ~ x
\dontrun{
plot(m)
}
dd <- sim(m,10000) ## Simulate 10000 observations from model
e <- estimate(m, dd) ## Estimate parameters
e
}
\keyword{models}
\keyword{regression}
