\name{covariance}
\alias{covariance}
\alias{covariance.lvm}
\alias{covariance<-}
\alias{covariance<-.lvm}
\alias{covfix}
\alias{covfix.lvm}
\alias{covfix<-}
\alias{covfix<-.lvm}
\title{Add covariance structure to Latent Variable Model}
\usage{
  \method{covariance}{lvm}(object, var1=NULL, var2=NULL,
    constrain=FALSE, pairwise=FALSE,...) <- value
}
\arguments{
  \item{object}{\code{lvm}-object}

  \item{var1}{Vector of variables names (or formula)}

  \item{var2}{Vector of variables names (or formula)
  defining pairwise covariance between \code{var1} and
  \code{var2})}

  \item{constrain}{Define non-linear parameter constraints
  to ensure positive definite structure}

  \item{pairwise}{If TRUE and \code{var2} is omitted then
  pairwise correlation is added between all variables in
  \code{var1}}

  \item{\dots}{Additional arguments to be passed to the low
  level functions}

  \item{value}{List of parameter values or (if \code{var1}
  is unspecified)}
}
\value{
  A \code{lvm}-object
}
\description{
  Define covariances between residual terms in a
  \code{lvm}-object.
}
\details{
  The \code{covariance} function is used to specify
  correlation structure between residual terms of a latent
  variable model, using a formula syntax.

  For instance, a multivariate model with three response
  variables,

  \deqn{Y_1 = \mu_1 + \epsilon_1}

  \deqn{Y_2 = \mu_2 + \epsilon_2}

  \deqn{Y_3 = \mu_3 + \epsilon_3}

  can be specified as

  \code{m <- lvm(~y1+y2+y3)}

  Pr. default the two variables are assumed to be
  independent. To add a covariance parameter \eqn{r =
  cov(\epsilon_1,\epsilon_2)}, we execute the following
  code

  \code{covariance(m) <- y1 ~ f(y2,r)}

  The special function \code{f} and its second argument
  could be omitted thus assigning an unique parameter the
  covariance between \code{y1} and \code{y2}.

  Similarily the marginal variance of the two response
  variables can be fixed to be identical (\eqn{var(Y_i)=v})
  via

  \code{covariance(m) <- c(y1,y2,y3) ~ f(v)}

  To specify a completely unstructured covariance
  structure, we can call

  \code{covariance(m) <- ~y1+y2+y3}

  All the parameter values of the linear constraints can be
  given as the right handside expression of the assigment
  function \code{covariance<-} if the first (and possibly
  second) argument is defined as well. E.g:

  \code{covariance(m,y1~y1+y2) <- list("a1","b1")}

  \code{covariance(m,~y2+y3) <- list("a2",2)}

  Defines

  \deqn{var(\epsilon_1) = a1}

  \deqn{var(\epsilon_2) = a2}

  \deqn{var(\epsilon_3) = 2}

  \deqn{cov(\epsilon_1,\epsilon_2) = b1}

  Parameter constraints can be cleared by fixing the
  relevant parameters to \code{NA} (see also the
  \code{regression} method).

  The function \code{covariance} (called without additional
  arguments) can be used to inspect the covariance
  constraints of a \code{lvm}-object.
}
\examples{
m <- lvm()
### Define covariance between residuals terms of y1 and y2
covariance(m) <- y1~y2
covariance(m) <- c(y1,y2)~f(v) ## Same marginal variance
covariance(m) ## Examine covariance structure
}
\author{
  Klaus K. Holst
}
\seealso{
  \code{\link{regression<-}}, \code{\link{intercept<-}},
  \code{\link{constrain<-}} \code{\link{parameter<-}},
  \code{\link{latent<-}}, \code{\link{cancel<-}},
  \code{\link{kill<-}}
}
\keyword{models}
\keyword{regression}

