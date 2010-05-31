\name{kill}
\alias{kill}
\alias{kill<-}
\alias{kill.lvm}
\alias{cancel}
\alias{cancel}
\alias{cancel<-}
\alias{cancel.lvm}
\alias{addvar<-}

\title{ Add or remove variables from latent variables. Remove associations from latent variable model }
\description{
   Add or remove variables from a latent variable model, or remove
   previously specified associations among variables.
}
\usage{
kill(x, ...) <- value
cancel(x, ...) <- value
addvar(x, ...) <- value
}
\arguments{
  \item{x}{\code{lvm}-object}
  \item{value}{Vector of variables or formula specifying which nodes to
    remove/add (or associations between nodes).}
  \item{\dots}{...}
}
\value{
  A \code{lvm}-object
}
\author{Klaus K. Holst}
\examples{
m <- lvm()
addvar(m) <- ~y1+y2+x
covariance(m) <- y1~y2
regression(m) <- c(y1,y2) ~ x
### Cancel the covariance between the residuals of y1 and y2
cancel(m) <- y1~y2
### Remove y2 from the model 
kill(m) <- ~y2
}
\keyword{models}
\keyword{regression}
