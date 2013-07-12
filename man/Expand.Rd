\name{Expand}
\alias{Expand}
\title{Create a Data Frame from All Combinations of Factors}
\usage{
  Expand(x, ...)
}
\arguments{
  \item{x}{Data.frame}

  \item{...}{vectors, factors or a list containing these}
}
\description{
  Create a Data Frame from All Combinations of Factors
}
\details{
  Simple wrapper of the 'expand.grid' function
}
\examples{
dd <- Expand(iris, Sepal.Length=2:8, Species="virginica")
str(dd)
}
\author{
  Klaus K. Holst
}

