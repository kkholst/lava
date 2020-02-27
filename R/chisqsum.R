rchisqsum <- function(n,lambda) {
    p <- length(lambda)
    X2 <- matrix(rnorm(n*p)^2,ncol=p) ## Chi-squared (df=1)
    res <- numeric(n)
    for (i in seq(p)) {
        res <- res + X2[,i]*lambda[i]
    }
    return(res)
}

pchisqsum <- function(x, lambda=1, B=1e6, seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    y <- rchisqsum(B,lambda)
    mean(y<=x)
}
