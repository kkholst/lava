\name{compare}
\alias{compare}
\alias{contrmat}
\title{Statistical tests

Performs Likelihood-ratio, Wald and score tests}
\usage{
compare(object, ...)
}
\arguments{
  \item{object}{\code{lvmfit}-object}

  \item{\dots}{Additional arguments to low-level functions}
}
\value{
Matrix of test-statistics and p-values
}
\description{
Statistical tests

Performs Likelihood-ratio, Wald and score tests
}
\examples{
m <- lvm();
regression(m) <- c(y1,y2,y3) ~ eta; latent(m) <- ~eta
regression(m) <- eta ~ x
m2 <- regression(m, c(y3,eta) ~ x)
set.seed(1)
d <- sim(m,1000)
e <- estimate(m,d)
e2 <- estimate(m2,d)

compare(e)

compare(e,e2) ## LRT, H0: y3<-x=0
compare(e,scoretest=y3~x) ## Score-test, H0: y3~x=0
compare(e2,par=c("y3~x")) ## Wald-test, H0: y3~x=0

B <- diag(2); colnames(B) <- c("y2~eta","y3~eta")
compare(e2,contrast=B,null=c(1,1))

B <- rep(0,length(coef(e2))); B[1:3] <- 1
compare(e2,contrast=B)

compare(e,scoretest=list(y3~x,y2~x))
}
\author{
Klaus K. Holst
}
\seealso{
\code{\link{modelsearch}}, \code{\link{equivalence}}
}
\keyword{htest}

