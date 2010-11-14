\name{gof}
\alias{gof}
\alias{gof.lvmfit}
\alias{moments}
\alias{moments.lvm}
\alias{information}
\alias{information.lvmfit}
\alias{score}
\alias{score.lvmfit}
\alias{logLik.lvmfit}
\title{ Extract model summaries and GOF statistics for model object }
\description{
  Calculates various GOF statistics for model object including global
  chi-squared test statistic and AIC. Extract model-specific mean and
  variance structure, residuals and various predicitions.
}
\usage{
gof(object, ...)
\method{gof}{lvmfit}(object, chisq=FALSE, ...)
moments(x,...)
\method{moments}{lvm}(x, p, debug=FALSE, conditional=FALSE, data=NULL, ...) 
\method{logLik}{lvmfit}(object, p=coef(object),
                      data=model.frame(object),
                      model=object$estimator,
                      weight=Weight(object),
                          ...)
\method{score}{lvmfit}(x, data=model.frame(x), p=pars(x), model=x$estimator, weight=Weight(x), ...)
\method{information}{lvmfit}(x,p=pars(x),n=x$data$n,data=model.frame(x),model=x$estimator,weight=Weight(x),...)
}
\arguments{
  \item{object}{ Model object }
  \item{x}{ Model object }
  \item{p}{ Parameter vector used to calculate statistics }
  \item{data}{ Data.frame to use }
  \item{weight}{ Optional weight matrix }
  \item{n}{ Number of observations }
  \item{conditional}{If TRUE the conditional moments given the covariates are
    calculated. Otherwise the joint moments are calculated }
  \item{model}{ String defining estimator, e.g. "gaussian" (see \code{estimate})}
  \item{debug}{Debugging only}
  \item{chisq}{ Boolean indicating whether to calculate chi-squared
    goodness-of-fit (always TRUE for estimator='gaussian')}  
  \item{\dots}{Additional arguments to be passed to the low level functions}
}
\value{
  A \code{htest}-object.
}
\author{ Klaus K. Holst }

\keyword{ methods }
\keyword{ models }
