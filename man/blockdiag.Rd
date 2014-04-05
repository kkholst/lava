\name{blockdiag}
\alias{blockdiag}
\title{Combine matrices to block diagonal structure}
\usage{
blockdiag(x, ..., pad = 0)
}
\arguments{
  \item{x}{Matrix}

  \item{\dots}{Additional matrices}

  \item{pad}{Vyalue outside block-diagonal}
}
\description{
Combine matrices to block diagonal structure
}
\examples{
A <- diag(3)+1
blockdiag(A,A,A,pad=NA)
}
\author{
Klaus K. Holst
}

