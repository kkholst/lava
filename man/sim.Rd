\name{sim}
\alias{Gamma.lvm}
\alias{aalenExponential.lvm}
\alias{binomial.lvm}
\alias{categorical}
\alias{categorical<-}
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
\alias{gamma.lvm}
\alias{gaussian.lvm}
\alias{heavytail}
\alias{heavytail<-}
\alias{loggamma.lvm}
\alias{logit.lvm}
\alias{lognormal.lvm}
\alias{normal.lvm}
\alias{ones.lvm}
\alias{poisson.lvm}
\alias{probit.lvm}
\alias{sequence.lvm}
\alias{sim}
\alias{sim.lvm}
\alias{sim.lvmfit}
\alias{simulate.lvm}
\alias{simulate.lvmfit}
\alias{student.lvm}
\alias{threshold.lvm}
\alias{transform.lvm}
\alias{transform<-}
\alias{transform<-.lvm}
\alias{uniform.lvm}
\alias{weibull.lvm}
\title{Simulate model}
\usage{
\method{sim}{lvm}(x, n = 100, p = NULL, normal = FALSE, cond = FALSE,
sigma = 1, rho = 0.5, X, unlink=FALSE, ...)
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
Simulate data from a general SEM model including non-linear
effects and general link and distribution of variables.
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
d <- sim(m,1e4,p=c(y=-1,"y~x"=2,z=1))
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


##################################################
### Non-random variables
##################################################
m <- lvm()
distribution(m,~x+z+v) <- list(sequence.lvm(0,5),## Seq. 0 to 5 by 1/n
                               ones.lvm(),       ## Vector of ones
                               ones.lvm(0.5))    ##  0.8n 0, 0.2n 1
sim(m,10)


##################################################
### Cox model
### piecewise constant hazard
################################################

m <- lvm(t~x)
rates <- c(1,0.5); cuts <- c(0,5)
## Constant rate: 1 in [0,5), 0.5 in [5,Inf)
distribution(m,~t) <- coxExponential.lvm(rate=rates,timecut=cuts)


\dontrun{
    d <- sim(m,2e4,p=c("t~x"=0.1)); d$status <- TRUE
    plot(timereg::aalen(survival::Surv(t,status)~x,data=d,
                        resample.iid=0,robust=0),spec=1)
    L <- approxfun(c(cuts,max(d$t)),f=1,
                   cumsum(c(0,rates*diff(c(cuts,max(d$t))))),
                   method="linear")
    curve(L,0,100,add=TRUE,col="blue")
}


##################################################
### Cox model
### piecewise constant hazard, gamma frailty
##################################################

m <- lvm(y~x+z)
rates <- c(0.3,0.5); cuts <- c(0,5)
distribution(m,~y+z) <- list(coxExponential.lvm(rate=rates,timecut=cuts),
                             loggamma.lvm(rate=1,shape=1))
\dontrun{
    d <- sim(m,2e4,p=c("y~x"=0,"y~z"=0)); d$status <- TRUE
    plot(timereg::aalen(survival::Surv(y,status)~x,data=d,
                        resample.iid=0,robust=0),spec=1)
    L <- approxfun(c(cuts,max(d$y)),f=1,
                   cumsum(c(0,rates*diff(c(cuts,max(d$y))))),
                   method="linear")
    curve(L,0,100,add=TRUE,col="blue")
}

## Equivalent via transform (here with Aalens additive hazard model)
m <- lvm(y~x)
distribution(m,~y) <- aalenExponential.lvm(rate=rates,timecut=cuts)
distribution(m,~z) <- Gamma.lvm(rate=1,shape=1)
transform(m,t~y+z) <- prod
sim(m,10)

## Shared frailty
m <- lvm(c(t1,t2)~x+z)
rates <- c(1,0.5); cuts <- c(0,5)
distribution(m,~y) <- aalenExponential.lvm(rate=rates,timecut=cuts)
distribution(m,~z) <- loggamma.lvm(rate=1,shape=1)
\dontrun{
    mets::fast.reshape(sim(m,100),varying="t")
}

##################################################
### Categorical predictor
##################################################

##library(mets)
m <- lvm()
m <- categorical(m,y~x,K=4)
d <- sim(m,1000,p=c('y~x:2'=4))
lm(y~-1+factor(x),data=d)


m <- lvm()
regression(m) <- z~x
m <- categorical(m,y~x,K=4,p=c(0.1,0.2,0.3))
with(sim(m,1e4),table(x))/1e4

m <- categorical(m,y~z,K=3)
d <- sim(m,1e4,p=c('y~z:0'=3,'y~x:2'=4))
lm(y~factor(z)+factor(x),data=d)
with(d,table(x))/nrow(d)
}
\author{
Klaus K. Holst
}
\keyword{datagen}
\keyword{models}
\keyword{regression}

