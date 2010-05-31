\name{vars}
\alias{vars}
\alias{vars.lvm}
\alias{vars.lvmfit}
\alias{latent}
\alias{latent<-}
\alias{latent.lvm}
\alias{latent<-.lvm}
\alias{latent.lvmfit}
\alias{latent.multigroup}
\alias{manifest}
\alias{manifest.lvm}
\alias{manifest.lvmfit}
\alias{manifest.multigroup}
\alias{exogenous}
\alias{exogenous.lvm}
\alias{exogenous<-.lvm}
\alias{exogenous.lvmfit}
\alias{exogenous.multigroup}
\alias{endogenous}
\alias{endogenous.lvm}
\alias{endogenous.lvmfit}
\alias{endogenous.multigroup}

\title{Extract variable names from latent variable model}
\description{
  Extract exogenous variables (predictors), endogenous variables
  (outcomes), latent variables (random effects), manifest (observed) variables from a
  \code{lvm} object.
}
\usage{
vars(x,...)
endogenous(x,...)
exogenous(x,...)
manifest(x,...)
latent(x,...)
\method{exogenous}{lvm}(x,...) <- value
\method{exogenous}{lvm}(x,latent=FALSE,index=TRUE,...)
\method{latent}{lvm}(x,clear=FALSE,...) <- value
}
\arguments{
  \item{x}{\code{lvm}-object}
  \item{latent}{Logical defining whether latent variables without parents
  should be included in the result}
  \item{index}{For internal use only}
  \item{clear}{Logical indicating whether to add or remove latent
    variable status}
  \item{value}{Formula or character vector of variable names.}
  \item{\dots}{}
}
\details{
  \code{vars} returns all variables of the \code{lvm}-object including
  manifest and latent variables. Similarily \code{manifest} and
  \code{latent} returns the observered resp. latent variables of the model.
  \code{exogenous} returns all manifest
  variables without parents, e.g. covariates in the model, however the argument
  \code{latent=TRUE} can be used to also include latent variables
  without parents in the result. Pr. default \code{lava} will not
  include the parameters of the exogenous variables in the optimisation
  routine during estimation (likelihood of the remaining observered
  variables conditional on the covariates), however this behaviour can
  be altered via the assignment function \code{exogenous<-} telling
  \code{lava} which subset of (valid) variables to condition on.
  Finally \code{latent} returns a vector with the names of the latent
  variables in \code{x}. The assigment function \code{latent<-} can be
  used to change the latent status of variables in the model.
}
\value{
  Vector of variable names.
}
\author{Klaus K. Holst}

\seealso{\code{\link{endogenous}}, \code{\link{manifest}},
  \code{\link{latent}}, \code{\link{exogenous}}, \code{\link{vars}}
}
\examples{
g <- lvm(eta1 ~ x1+x2)
regression(g) <- c(y1,y2,y3) ~ eta1
latent(g) <- ~eta
endogenous(g)
exogenous(g)
identical(latent(g), setdiff(manifest(g), vars(g)))
}
\keyword{models}
\keyword{regression}
