% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measurement.error.R
\name{measurement.error}
\alias{measurement.error}
\title{Two-stage (non-linear) measurement error}
\usage{
measurement.error(
  model1,
  formula,
  data = parent.frame(),
  predictfun = function(mu, var, data, ...) mu[, 1]^2 + var[1],
  id1,
  id2,
  ...
)
}
\arguments{
\item{model1}{Stage 1 model}

\item{formula}{Formula specifying observed covariates in stage 2 model}

\item{data}{data.frame}

\item{predictfun}{Predictions to be used in stage 2}

\item{id1}{Optional id-vector of stage 1}

\item{id2}{Optional id-vector of stage 2}

\item{...}{Additional arguments to lower level functions}
}
\description{
Two-stage measurement error
}
\examples{
m <- lvm(c(y1,y2,y3)~u,c(y3,y4,y5)~v,u~~v,c(u,v)~x)
transform(m,u2~u) <- function(x) x^2
transform(m,uv~u+v) <- prod
regression(m) <- z~u2+u+v+uv+x
set.seed(1)
d <- sim(m,1000,p=c("u,u"=1))

## Stage 1
m1 <- lvm(c(y1[0:s],y2[0:s],y3[0:s])~1*u,c(y3[0:s],y4[0:s],y5[0:s])~1*v,u~b*x,u~~v)
latent(m1) <- ~u+v
e1 <- estimate(m1,d)

pp <- function(mu,var,data,...) {
    cbind(u=mu[,"u"],u2=mu[,"u"]^2+var["u","u"],v=mu[,"v"],uv=mu[,"u"]*mu[,"v"]+var["u","v"])
}
(e <- measurement.error(e1, z~1+x, data=d, predictfun=pp))

## uu <- seq(-1,1,length.out=100)
## pp <- estimate(e,function(p,...) p["(Intercept)"]+p["u"]*uu+p["u2"]*uu^2)$coefmat
if (interactive()) {
    plot(e,intercept=TRUE,line=0)

    f <- function(p) p[1]+p["u"]*u+p["u2"]*u^2
    u <- seq(-1,1,length.out=100)
    plot(e, f, data=data.frame(u), ylim=c(-.5,2.5))
}
}
\seealso{
stack.estimate
}
