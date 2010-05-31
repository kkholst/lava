\name{confint.lvmfit}
\alias{confint.lvmfit}
\title{ Calculate confidence limits for parameters }
\description{
  Calculate Wald og Likelihood based (profile likelihood) confidence intervals 
}
\usage{
\method{confint}{lvmfit}(object,parm,level=0.95,profile=FALSE,...) 
}
\arguments{
  \item{object}{\code{lvm}-object.}
  \item{parm}{Index of which parameters to calculate confidence limits
    for.}
  \item{level}{Confidence level}
  \item{profile}{Logical expression defining whether to calculate
    confidence limits via the profile log likelihood}
  \item{\dots}{}
}
\value{
  A 2xp matrix with columns of lower and upper confidence limits
}
\details{
  Calculates either Wald confidence limits:
  \deqn{\hat{\theta} \pm z_{\alpha/2}*\hat\sigma_{\hat\theta}}
  or profile likelihood confidence limits, defined as the set of value
  \eqn{\tau}:
  \deqn{logLik(\hat\theta_{\tau},\tau)-logLik(\hat\theta)< q_{\alpha}/2}

  where \eqn{q_{\alpha}} is the \eqn{\alpha} fractile of the
  \eqn{\chi^2_1} distribution, and \eqn{\hat\theta_{\tau}} are obtained by
  maximizing the log-likelihood with tau being fixed.  
}
\seealso{
  \code{\link{bootstrap.lvm}}  
}
\examples{
m <- lvm(y~x)
d <- sim(m,100)
e <- estimate(y~x, d)
confint(e,3,profile=TRUE)
confint(e,3)
\dontrun{
B <- bootstrap(e,R=100)
B
}
}
\author{Klaus K. Holst}
\keyword{models}
\keyword{regression}
