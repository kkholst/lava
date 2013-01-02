\name{blockdiag}
\alias{blockdiag}
\title{Combine matrices to block diagonal structure}
\usage{
  blockdiag(x, ..., pad = 0)
}
\arguments{
  \item{x}{Matrix}

  \item{...}{Additional matrices}

  \item{pad}{Value outside block-diagonal}
}
\description{
  Combine matrices to block diagonal structure
}
\examples{
A <- diag(3)+1
A\%+\%A
blockdiag(A,A,A,pad=NA)
}
\author{
  Klaus K. Holst
}

