pzmax <- function(alpha,S) {
    ##P(Zmax > z) Family wise error rate, Zmax = max |Z_i|
    k <- nrow(S)
    z <- qnorm(1-alpha/2)    
    1-mets::pmvn(lower=rep(-z,k),upper=rep(z,k),sigma=cov2cor(S))
}

##' @export
p.correct <- function(object,idx,alpha=0.05) {
    S <- vcov(object); if (!missing(idx)) S <- S[idx,idx,drop=FALSE]
    f <- function(a) pzmax(a,S)-alpha
    uniroot(f,lower=0,upper=0.05)$root
}

##' @export
closed.testing <- function(object,idx=seq_along(coef(object)),...) {    
    B <- diag(nrow=length(idx))
    e <- estimate(object,keep=idx)
    combs <- pvals <- c()
    for (i in seq_along(idx)) {        
        co <- combn(length(idx),i)
        pp <- numeric(ncol(co))
        for (j in seq_along(pp)) {
            pp[j] <- compare(e,contrast=B[co[,j],,drop=FALSE])$p.value
        }
        combs <- c(combs,list(co))
        pvals <- c(pvals,list(pp))
    }
    pmax <- c()
    for (k in seq_along(idx)) {
        pk <- c()
        for (i in seq_along(idx)) {
            cols <- apply(combs[[i]],2,function(x) k%in%x)
            pk <- c(pk,pvals[[i]][which(cols)])
        }
        pmax <- c(pmax,max(pk))
    }
    return(pmax)
}
