\name{PD}
\alias{PD}
\title{Dose response calculation for binomial regression models}
\usage{
PD(model, intercept = 1, slope = 2, prob = NULL, x, level = 0.5,
  ci.level = 0.95, vcov, family, EB = NULL)
}
\arguments{
  \item{model}{Model object or vector of parameter
  estimates}

  \item{intercept}{Index of intercept parameters}

  \item{slope}{Index of intercept parameters}

  \item{prob}{Index of mixture parameters (only relevant
  for \code{curereg} models)}

  \item{x}{Optional weights
  length(x)=length(intercept)+length(slope)+length(prob)}

  \item{level}{Probability at which level to calculate
  dose}

  \item{ci.level}{Level of confidence limits}

  \item{vcov}{Optional estimate of variance matrix of
  parameter estimates}

  \item{family}{Optional distributional family argument}

  \item{EB}{Optional ratio of treatment effect and adverse
  effects used to find optimal dose (regret-function
  argument)}
}
\description{
Dose response calculation for binomial regression models
}
\author{
Klaus K. Holst
}

