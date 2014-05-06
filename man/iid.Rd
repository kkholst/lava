\name{iid}
\alias{iid}
\alias{iid.default}
\title{Extract i.i.d. decomposition (influence function) from model object}
\usage{
  iid(x,...)

  \method{iid}{default}(x,score.deriv,id,...)
}
\arguments{
  \item{x}{model object}

  \item{id}{id/cluster variable (optional)}

  \item{score.deriv}{(optional) derivative of mean score
  function}

  \item{...}{additional arguments}
}
\description{
  Extract i.i.d. decomposition (influence function) from
  model object
}
\examples{
m <- lvm(y~x+z)
distribution(m, ~y+z) <- binomial.lvm("logit")
d <- sim(m,1e3)
g <- glm(y~x+z,data=d,family=binomial)
crossprod(iid(g))
}

