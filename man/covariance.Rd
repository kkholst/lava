\name{covariance}
\alias{covariance}
\alias{covariance<-}
\alias{covariance.lvm}
\alias{covariance<-.lvm}
\alias{covfix}
\alias{covfix<-}
\alias{covfix.lvm}
\title{ Add covariance structure to Latent Variable Model }
\description{
   Define covariances between residual terms in a \code{lvm}-object.
}
\usage{
\method{covariance}{lvm}(object, var1=NULL, var2=NULL, ...) <- value
\method{covfix}{lvm}(object, ...)

}
\arguments{
  \item{object}{\code{lvm}-object}
  \item{var1}{Vector of variables names
    between (or formula)}
  \item{var2}{Vector of variables names (or formula) defining pairwise covariance between \code{var1} and \code{var2})}
  \item{\dots}{}
  \item{value}{List of parameter values or (if \code{var1} is
    unspecified) a }
}
\details{
  The \code{covariance} function is used to specify correlation
  structure between residual terms of a latent variable model, using
  a formula syntax.
  
  For instance, a multivariate model with two response variables,
  \deqn{Y_1 = \mu_1 + \epsilon_1}{y1 = mu1 + epsilon1}
  \deqn{Y_2 = \mu_2 + \epsilon_2}{y2 = mu2 + epsilon2}
  \deqn{Y_3 = \mu_3 + \epsilon_3}{y3 = mu3 + epsilon3}
  can be specified as

  \code{m <- lvm(~y1+y2+y3)}

  Pr. default the two variables are assumed to be independent. To add a
  covariance parameter \eqn{r = cov(\epsilon_1,\epsilon_2)}, we execute
  the following code
  
  \code{covariance(m) <- y1 ~ f(y2,r)}

  The special function \code{f} and its second argument could be omitted
  thus assigning an unique parameter the covariance between \code{y1}
  and \code{y2}.

  Similarily the marginal variance of the two response variables can be
  fixed to be identical (\eqn{var(Y_i)=v}) via

  \code{covariance(m) <- c(y1,y2,y3) ~ f(v)}

  To specify a completely unstructured covariance structure, we can call

  \code{covariance(m) <- ~y1+y2+y3}

  All the parameter values of the linear constraints can be
  given as the right handside expression of the assigment
  function \code{covariance<-} if the first (and
  possibly second) argument is defined as well. E.g:

  \code{covariance(m,y1~y1+y2) <- list("a1","b1")}
  \code{covariance(m,~y2+y3) <- list("a2",2)}

  Defines
  \deqn{var(\epsilon_1) = a1}
  \deqn{var(\epsilon_2) = a2}
  \deqn{var(\epsilon_3) = 2}
  \deqn{cov(\epsilon_1,\epsilon_2) = b1}  
  
  The assigment function \code{covfix<-} can like \code{covariance<-} be
  used to specify covariance constraints, but conveniently allows the
  parameter values to be given as a list of values (numeric and
  characters).

  Parameter constraints can be cleared by fixing the relevant parameters
  to \code{NA} (see also the \code{regression} method).
  
  The function \code{covfix} (called without
  additional arguments) can be used to inspect the covariance
  constraints of a \code{lvm}-object. 
  
  }
\value{
  A \code{lvm}-object
}
\author{Klaus K. Holst}
\seealso{
  \code{\link{regression<-}},
  \code{\link{intercept<-}},
  \code{\link{constrain<-}}
  \code{\link{parameter<-}},
  \code{\link{latent<-}},
  \code{\link{cancel<-}},
  \code{\link{kill<-}}  
}
\examples{
m <- lvm()
### Define covariance between residuals terms of y1 and y2
covariance(m) <- y1~y2 
covariance(m) <- c(y1,y2)~f(v) ## Same marginal variance
covfix(m) ## Examine covariance structure
}
\keyword{models}
\keyword{regression}
