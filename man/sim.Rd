\name{sim}
\alias{binomial.lvm}
\alias{coxExponential.lvm}
\alias{coxGompertz.lvm}
\alias{coxWeibull.lvm}
\alias{distribution}
\alias{distribution.lvm}
\alias{distribution<-}
\alias{distribution<-.lvm}
\alias{functional}
\alias{functional.lvm}
\alias{functional<-}
\alias{functional<-.lvm}
\alias{Gamma.lvm}
\alias{gaussian.lvm}
\alias{heavytail}
\alias{heavytail<-}
\alias{logit.lvm}
\alias{normal.lvm}
\alias{poisson.lvm}
\alias{probit.lvm}
\alias{sim}
\alias{sim.lvm}
\alias{sim.lvmfit}
\alias{simulate.lvm}
\alias{simulate.lvmfit}
\alias{student.lvm}
\alias{transform.lvm}
\alias{transform<-}
\alias{transform<-.lvm}
\alias{uniform.lvm}
\alias{weibull.lvm}
\title{Simulate model}
\usage{
  \method{sim}{lvm}(x, n = 100, p = NULL, normal = FALSE,
    cond = FALSE, sigma = 1, rho = 0.5, X, unlink=FALSE,
    ...)
}
\arguments{
  \item{x}{Model object}

  \item{n}{Number of simulated values/individuals}

  \item{p}{Parameter value (optional)}

  \item{normal}{Logical indicating whether to simulate data
  from a multivariate normal distribution conditional on
  exogenous variables hence ignoring
  functional/distribution definition}

  \item{cond}{for internal use}

  \item{sigma}{Default residual variance (1)}

  \item{rho}{Default covariance parameter (0.5)}

  \item{X}{Optional matrix of covariates}

  \item{unlink}{Return Inverse link transformed data}

  \item{\dots}{Additional arguments to be passed to the low
  level functions}
}
\description{
  Simulate data from a general SEM model including
  non-linear effects and general link and distribution of
  variables.
}
\details{
  E.g. \eqn{ E(y|x) = 2*x^2 } could be specified as

  \code{regression(m, "y3", fn=function(x) x^2) <- "x$2"}
}
\examples{
##################################################
## Logistic regression
##################################################
m <- lvm(y~x+z)
regression(m) <- x~z
distribution(m,~y+z) <- binomial.lvm("logit")
d <- sim(m,1e3)
head(d)

e <- estimate(m,d,estimator="glm")
e
## Simulate a few observation from estimated model
sim(e,n=5)

##################################################
## Poisson
##################################################
distribution(m,~y) <- poisson.lvm()
d <- sim(m,1e4,p=c(y=-1,"y<-x"=2,z=1))
head(d)
estimate(m,d,estimator="glm")
mean(d$z); lava:::expit(1)

summary(lm(y~x,sim(lvm(y[1:2]~4*x),1e3)))

##################################################
### Gamma distribution
##################################################
m <- lvm(y~x)
distribution(m,~y+x) <- list(Gamma.lvm(shape=2),binomial.lvm())
intercept(m,~y) <- 0.5
d <- sim(m,1e4)
summary(g <- glm(y~x,family=Gamma(),data=d))
\dontrun{MASS::gamma.shape(g)}

args(lava::Gamma.lvm)
distribution(m,~y) <- Gamma.lvm(shape=2,log=TRUE)
sim(m,10,p=c(y=0.5))[,"y"]

##################################################
### Transform
##################################################

m <- lvm()
transform(m,xz~x+z) <- function(x) x[1]*(x[2]>0)
regression(m) <- y~x+z+xz
d <- sim(m,1e3)
summary(lm(y~x+z + x*I(z>0),d))
}
\author{
  Klaus K. Holst
}
\keyword{datagen}
\keyword{models}
\keyword{regression}

