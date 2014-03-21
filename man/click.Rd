\name{click}
\alias{click}
\alias{idplot}
\title{Identify points on plot}
\usage{
  \method{click}{default}(x, y=NULL, label=TRUE,
    n=length(x), pch=19, col="orange", cex=3, ...)

  idplot(x,y,...,id=list())
}
\arguments{
  \item{object}{}

  \item{...}{}

  \item{#x}{X coordinates}

  \item{#y}{Y coordinates}

  \item{#label}{Should labels be added?}

  \item{#n}{Max number of inputs to expect}

  \item{#pch}{Symbol}

  \item{#col}{Color}

  \item{#cex}{Size}

  \item{#id}{List of arguments parsed to \code{Id}
  function}

  \item{#\dots}{Additional arguments parsed to \code{plot}
  function}
}
\description{
  Extension of the \code{identify} function
}
\details{
  For the usual 'X11' device the identification process is
  terminated by pressing any mouse button other than the
  first. For the 'quartz' device the process is terminated
  by pressing either the pop-up menu equivalent (usually
  second mouse button or 'Ctrl'-click) or the 'ESC' key.
}
\examples{
\donttest{
n <- 10; x <- seq(1:n); y <- runif(n)
plot(y ~ x); click(x,y)

data(iris)
l <- lm(Sepal.Length ~ Sepal.Width*Species,iris)
res <- plotConf(l,var2="Species")## ylim=c(6,8), xlim=c(2.5,3.3))
with(res, click(x,y))

with(iris, idplot(Sepal.Length,Petal.Length))
}
}
\author{
  Klaus K. Holst
}
\seealso{
  \code{\link{idplot}}, \code{identify}
}
\keyword{iplot}

