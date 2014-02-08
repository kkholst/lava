\name{multinomial}
\alias{gkgamma}
\alias{kappa.multinomial}
\alias{kappa.table}
\alias{multinomial}
\title{Estimate probabilities in contingency table}
\usage{
  multinomial(x, ...)
}
\arguments{
  \item{x}{Matrix or data.frame with observations (1 or 2
  columns)}

  \item{...}{Additional arguments to lower-level functions}
}
\description{
  Estimate probabilities in contingency table
}
\examples{
set.seed(1)
breaks <- c(-Inf,-1,0,Inf)
m <- lvm(); covariance(m,pairwise=TRUE) <- ~y1+y2+y3+y4
d <- transform(sim(m,5e2),
              z1=cut(y1,breaks=breaks),
              z2=cut(y2,breaks=breaks),
              z3=cut(y3,breaks=breaks),
              z4=cut(y4,breaks=breaks))


multinomial(d[,5])
(a1 <- multinomial(d[,5:6]))
(K1 <- kappa(a1)) ## Cohen's kappa

K2 <- kappa(d[,7:8])
## Testing difference K1-K2:
estimate(merge(K1,K2),diff)

sqrt(vcov(K1)+vcov(K2)) ## Wrong std.err ignoring dependence

## Goodman-Kruskal's gamma
m2 <- lvm(); covariance(m2) <- y1~y2
breaks1 <- c(-Inf,-1,0,Inf)
breaks2 <- c(-Inf,0,Inf)
d2 <- transform(sim(m2,5e2),
              z1=cut(y1,breaks=breaks1),
              z2=cut(y2,breaks=breaks2))

(g1 <- gkgamma(d2[,3:4]))
## same as
\dontrun{
gkgamma(table(d2[,3:4]))
gkgamma(multinomial(d2[,3:4]))
}
}
\author{
  Klaus Kähler Holst
}

