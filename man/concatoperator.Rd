\name{concatoperator}
\alias{concatoperator}
\title{Concatenation operator}
\usage{
  %+%(x, y)
}
\arguments{
  \item{x}{First object}

  \item{y}{Second object}
}
\description{
  Concatenation operator
}
\details{
  For matrices a block-diagonal matrix is created. For
  latent variable models ('lvm' objects) \code{merge.lvm}
  is called. For all other data types he operator is a
  wrapper of \code{paste}.
}
\author{
  Klaus K. Holst
}

