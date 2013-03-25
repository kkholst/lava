\name{estimate.default}
\alias{estimate.default}
\title{Estimation of functional of parameters}
\usage{
  estimate.default(x, fun, data = model.frame(x), id1, id2,
    vcov = TRUE, level = 0.95, iid = FALSE, contrast, null,
    ...)
}
\arguments{
  \item{x}{model object (\code{glm}, \code{lvmfit}, ...)}

  \item{fun}{transformation of model parameters and
  optionally data.}

  \item{data}{\code{data.frame}}

  \item{id1}{(optional) id-variable corresponding to iid
  decomposition of model parameters}

  \item{id2}{(optional) id-variable of data.frame}

  \item{vcov}{If TRUE the inverse of \code{vcov(x)} is used
  to calculate ---}

  \item{level}{Level of confidence limits}

  \item{iid}{If TRUE the iid decompositions are returned
  instead of variance estimates}

  \item{contrast}{(optional) Contrast matrix for final Wald
  test}

  \item{null}{(optional) Null hypothesis to test}

  \item{...}{additional arguments to lower level functions}
}
\description{
  Estimation of functional of parameters
}
\examples{
## Simulation from logistic regression model
m <- lvm(y~x+z);
distribution(m,y~x) <- binomial.lvm("logit")
d <- sim(m,1000)
g <- glm(y~z+x,data=d,family=binomial())

## Plain estimates (robust standard errors)
estimate(g)
## Testing contrasts
estimate(g,null=0)
estimate(g,contrast=rbind(c(1,1,0),c(1,0,2)))
estimate(g,contrast=rbind(c(1,1,0),c(1,0,2)),null=c(1,2))

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

## Look at lava:::TN.curereg for example where gradient is supplied as attribute.
}

