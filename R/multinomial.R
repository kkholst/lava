
##' Estimate probabilities in contingency table
##'
##' @title Estimate probabilities in contingency table
##' @aliases multinomial kappa.multinomial kappa.table gkgamma
##' @param x Matrix or data.frame with observations (1 or 2 columns)
##' @param marginal If TRUE the marginals are estimated
##' @param transform Optional transformation of parameters (e.g., logit)
##' @param ... Additional arguments to lower-level functions
##' @export
##' @examples
##' set.seed(1)
##' breaks <- c(-Inf,-1,0,Inf)
##' m <- lvm(); covariance(m,pairwise=TRUE) <- ~y1+y2+y3+y4
##' d <- transform(sim(m,5e2),
##'               z1=cut(y1,breaks=breaks),
##'               z2=cut(y2,breaks=breaks),
##'               z3=cut(y3,breaks=breaks),
##'               z4=cut(y4,breaks=breaks))
##' 
##' multinomial(d[,5])
##' (a1 <- multinomial(d[,5:6]))
##' (K1 <- kappa(a1)) ## Cohen's kappa
##' 
##' K2 <- kappa(d[,7:8])
##' ## Testing difference K1-K2:
##' estimate(merge(K1,K2,id=TRUE),diff)
##' 
##' estimate(merge(K1,K2,id=FALSE),diff) ## Wrong std.err ignoring dependence
##' sqrt(vcov(K1)+vcov(K2))
##' 
##' ## Average of the two kappas:
##' estimate(merge(K1,K2,id=TRUE),function(x) mean(x))
##' estimate(merge(K1,K2,id=FALSE),function(x) mean(x)) ## Independence
##' ##'
##' ## Goodman-Kruskal's gamma
##' m2 <- lvm(); covariance(m2) <- y1~y2
##' breaks1 <- c(-Inf,-1,0,Inf)
##' breaks2 <- c(-Inf,0,Inf)
##' d2 <- transform(sim(m2,5e2),
##'               z1=cut(y1,breaks=breaks1),
##'               z2=cut(y2,breaks=breaks2))
##' 
##' (g1 <- gkgamma(d2[,3:4]))
##' ## same as
##' \dontrun{
##' gkgamma(table(d2[,3:4]))
##' gkgamma(multinomial(d2[,3:4]))
##' }
##' 
##' ##partial gamma
##' d2$x <- rbinom(nrow(d2),2,0.5)
##' gkgamma(z1~z2|x,data=d2)
##' @author Klaus K. Holst
multinomial <- function(x,marginal=FALSE,transform,...) {
    if (is.table(x)) x <- lava::Expand(x)
    if (NCOL(x)==1) {
        x <- as.factor(x)
        lev <- levels(x)
        k <- length(lev)
        n <- length(x)
        P <- table(x)/n
        iid <- matrix(0,n,k)
        for (i in seq(k)) {
            iid[,i] <- (1*(x==lev[i])-P[i])/n
        };
        coefs <- as.vector(P); names(coefs) <- paste0("p",seq(k))
        res <- list(coef=coefs,P=P,vcov=crossprod(iid),iid=iid,position=seq(k),levels=list(lev),data=x)
        class(res) <- "multinomial"
        return(res)
    }
    if (NCOL(x)!=2L) stop("Matrix or data.frame with one or two columns expected")
    x <- as.data.frame(x)
    x[,1] <- as.factor(x[,1])
    x[,2] <- as.factor(x[,2])
    lev1 <- levels(x[,1])
    lev2 <- levels(x[,2])
    k1 <- length(lev1)
    k2 <- length(lev2)
    M <- table(x)
    n <- sum(M)
    Pos <- P <- M/n
    iid <- matrix(0,n,k1*k2)
    for (j in seq(k2)) {
        for (i in seq(k1)) {
            pos <- (j-1)*k1+i
            iid[,pos] <- (x[,1]==lev1[i])*(x[,2]==lev2[j])-P[i,j]
            Pos[i,j] <- pos
        }
    }; iid <- iid/n
    coefs <- as.vector(P);
    names(coefs) <-  as.vector(outer(seq(k1),seq(k2),function(...) paste0("p",...)))
    position1 <- position2 <- NULL
    if (marginal) {
        p1 <- rowSums(P)
        p2 <- colSums(P)
        iid1 <- apply(Pos,1,function(x) rowSums(iid[,x]))
        iid2 <- apply(Pos,2,function(x) rowSums(iid[,x]))
        iid <- cbind(iid,iid1,iid2)
        names(p1) <- paste0("p",seq(k1),".")
        names(p2) <- paste0("p",".",seq(k2))
        coefs <- c(coefs,p1,p2)
        colnames(iid) <- names(coefs)
        position1 <- length(P)+seq(k1)
        position2 <- length(P)+k1+seq(k2)
    }
    if (!missing(transform)) {
        f <- function(p) do.call(transform,list(p))
        D <- diag(numDeriv::grad(f,coefs),ncol=length(coefs))
        coefs <- f(coefs)
        iid <- iid%*%t(D)
    }
    res <- list(coef=coefs,P=P,vcov=crossprod(iid),iid=iid, position=Pos, call=match.call(), levels=list(lev1,lev2), data=x,
                position1=position1,position2=position2 ## Position of marginals)
                )
    class(res) <- "multinomial"
    if (length(list(...))>0) {
        res <- structure(estimate(res,...),class=c("multinomial","estimate"))
    }
    res
}

##' @export
model.frame.multinomial <- function(formula,...) {
    formula$data
}

##' @export
iid.multinomial <- function(x,...) {
    x$iid
}

##' @export
coef.multinomial <- function(object,...) {
    object$coef
}

##' @export
vcov.multinomial <- function(object,...) {
    object$vcov
}

##' @export
print.multinomial <- function(x,...) {
    cat("Call: "); print(x$call)
    cat("\nEstimates:\n")
    print(x$P,quote=FALSE)
    cat("\n")
    print(estimate(NULL,coef=coef(x),vcov=vcov(x)))
    ## stderr <- diag(vcov(x))^.5
    ## StdErr <- x$position
    ## StdErr[] <- stderr[StdErr]
    ## cat("\nStd.Err:\n")
    ## print(StdErr,quote=FALSE)
    ## cat("\nPosition:\n")
    ## print(x$position,quote=FALSE)
}


