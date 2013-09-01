\name{estimate.default}
\alias{estimate.default}
\title{Estimation of functional of parameters}
\usage{
  \method{estimate}{default} (x, f, data = model.frame(x),
    id, subset, score.deriv, level = 0.95, iid = FALSE,
    contrast, null, vcov, coef, ...)
}
\arguments{
  \item{x}{model object (\code{glm}, \code{lvmfit}, ...)}

  \item{f}{transformation of model parameters and
  (optionally) data}

  \item{data}{\code{data.frame}}

  \item{id}{(optional) id-variable corresponding to iid
  decomposition of model parameters}

  \item{subset}{(optional) subset of data.frame on which to
  condition (logical expression or variable name)}

  \item{score.deriv}{(optional) derivative of mean score
  function}

  \item{level}{Level of confidence limits}

  \item{iid}{If TRUE the iid decompositions are returned
  instead of variance estimates}

  \item{contrast}{(optional) Contrast matrix for final Wald
  test}

  \item{null}{(optional) Null hypothesis to test}

  \item{vcov}{(optional) covariance matrix of parameter
  estimates (e.g. Wald-test)}

  \item{coef}{(optional) parameter coefficient}

  \item{...}{additional arguments to lower level functions}
}
\description{
  Estimation of functional of parameters. Wald tests,
  robust standard errors, cluster robust standard errors,
  LRT (when \code{f} is not a function)...
}
\examples{
## Simulation from logistic regression model
m <- lvm(y~x+z);
distribution(m,y~x) <- binomial.lvm("logit")
d <- sim(m,1000)
g <- glm(y~z+x,data=d,family=binomial())
g0 <- glm(y~1,data=d,family=binomial())

## LRT
estimate(g,g0)

## Plain estimates (robust standard errors)
estimate(g)

## Testing contrasts
estimate(g,null=0)
estimate(g,contrast=rbind(c(1,1,0),c(1,0,2)))
estimate(g,contrast=rbind(c(1,1,0),c(1,0,2)),null=c(1,2))
estimate(g,contrast=2:3) ## same as rbind(c(0,1,0),c(0,0,1))

## Transformations
estimate(g,function(p) p[1]+p[2])

## Multiple parameters
e <- estimate(g,function(p) c(p[1]+p[2],p[1]*p[2]))
e
vcov(e)

## Label new parameters
estimate(g,function(p) list("a1"=p[1]+p[2],"b1"=p[1]*p[2]))

## Marginalize
f <- function(p,data)
  list(p0=lava:::expit(p[1] + p[3]*data[,"z"]),
       p1=lava:::expit(p[1] + p[2] + p[3]*data[,"z"]))
e <- estimate(g, f)
e
estimate(e,diff)
estimate(e,contrast=cbind(1,1))
}

