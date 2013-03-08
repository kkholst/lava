\name{curereg}
\alias{curereg}
\title{Regression model for binomial data with unkown group of immortals}
\usage{
  curereg(formula, cureformula = ~1, data,
    fam = binomial(), offset = NULL, start, var = "robust",
    ...)
}
\arguments{
  \item{formula}{Formula specifying}

  \item{cureformula}{Formula for model of disease
  prevalence}

  \item{data}{data frame}

  \item{fam}{Distribution family (see the help page
  \code{family})}

  \item{offset}{Optional offset}

  \item{start}{Optional starting values}

  \item{var}{Type of variance (robust, expected, hessian,
  outer)}

  \item{...}{Additional arguments to lower level functions}
}
\description{
  Regression model for binomial data with unkown group of
  immortals
}
\author{
  Klaus K. Holst
}

