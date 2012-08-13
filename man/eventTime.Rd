\name{eventTime}
\alias{eventTime}
\title{Add an observed event time outcome to a latent variable model.}
\usage{
  eventTime(object, formula, eventName, ...)
}
\arguments{
  \item{object}{Model object}

  \item{formula}{Formula (see details)}

  \item{eventName}{Event names}

  \item{...}{Additional arguments to lower levels
  functions}
}
\description{
  For example, if the model 'm' includes latent event time
  variables are called 'T1' and 'T2' and 'C' is the end of
  follow-up (right censored), then one can specify
}
\details{
  \code{eventTime(object=m,formula=ObsTime~min(T1=a,T2=b,C=0,"ObsEvent"))}

  when data are simulated from the model one gets 2 new
  columns:

  - "ObsTime": the smallest of T1, T2 and C - "ObsEvent":
  'a' if T1 is smallest, 'b' if T2 is smallest and '0' if C
  is smallest

  Note that "ObsEvent" and "ObsTime" are names specified by
  the user.
}
\examples{
m <- lvm()
distribution(m,~X2) <- binomial.lvm()
regression(m) <- T1~f(X1,-.5)+f(X2,0.3)
regression(m) <- T2~f(X2,0.6)
distribution(m,~T1) <- coxWeibull.lvm(scale=1/100)
distribution(m,~T2) <- coxWeibull.lvm(scale=1/100)
distribution(m,~C) <- coxWeibull.lvm(scale=1/100)
m <- eventTime(m,time~min(T1=1,T2=2,C=0),"event")
d <- sim(m,10000)


m <- eventTime(m,cens.otime~min(T1,T2=E2,C=0),"event")
sim(m,10)
}
\author{
  Thomas A. Gerds
}
\keyword{models}
\keyword{regression}
\keyword{survival}

