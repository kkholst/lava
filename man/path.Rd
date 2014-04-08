\name{path}
\alias{effects}
\alias{effects.lvmfit}
\alias{path}
\alias{path.lvm}
\alias{totaleffects}
\title{Extract pathways in model graph}
\usage{
  \method{path}{lvm} (object, to = NULL, from, ...)
    \method{effects}{lvmfit} (object, to, from,
    silent=FALSE, ...)
}
\arguments{
  \item{object}{Model object (\code{lvm})}

  \item{to}{Outcome variable (string). Alternatively a
  formula specifying response and predictor in which case
  the argument \code{from} is ignored.}

  \item{from}{Response variable (string), not necessarily
  directly affected by \code{to}.}

  \item{silent}{Logical variable which indicates whether
  messages are turned on/off.}

  \item{\dots}{Additional arguments to be passed to the low
  level functions}
}
\value{
  If \code{object} is of class \code{lvmfit} a list with
  the following elements is returned \item{idx}{ A list
  where each element defines a possible pathway via a
  integer vector indicating the index of the visited nodes.
  } \item{V }{ A List of covariance matrices for each path.
  } \item{coef }{A list of parameters estimates for each
  path} \item{path }{A list where each element defines a
  possible pathway via a character vector naming the
  visited nodes in order.  } \item{edges }{Description of
  'comp2'}

  If \code{object} is of class \code{lvm} only the
  \code{path} element will be returned.

  The \code{effects} method returns an object of class
  \code{effects}.
}
\description{
  Extract all possible paths from one variable to another
  connected component in a latent variable model. In an
  estimated model the effect size is decomposed into
  direct, indirect and total effects including approximate
  standard errors.
}
\note{
  For a \code{lvmfit}-object the parameters estimates and
  their corresponding covariance matrix are also returned.
  The \code{effects}-function additionally calculates the
  total and indirect effects with approximate standard
  errors
}
\examples{
m <- lvm(c(y1,y2,y3)~eta)
regression(m) <- y2~x1
latent(m) <- ~eta
regression(m) <- eta~x1+x2
d <- sim(m,500)
e <- estimate(m,d)

path(Model(e),y2~x1)
parents(Model(e), ~y2)
children(Model(e), ~x2)
children(Model(e), ~x2+eta)
effects(e,y2~x1)
}
\author{
  Klaus K. Holst
}
\seealso{
  \code{children}, \code{parents}
}
\keyword{graphs}
\keyword{methods}
\keyword{models}

