##' Estimate probabilities in contingency table
##'
##' @title Estimate probabilities in contingency table
##' @aliases multinomial kappa.multinomial
##' @param x Matrix or data.frame with observations
##' @param ... Additional arguments to lower-level functions
##' @export
##' @examples
##' set.seed(123)
##' breaks <- c(-Inf,-1,0,Inf)
##' m <- lvm(); covariance(m,pairwise=TRUE) <- ~y1+y2+y3+y4
##' d <- transform(sim(m,5e2),
##'               z1=cut(y1,breaks=breaks),
##'               z2=cut(y2,breaks=breaks),
##'               z3=cut(y3,breaks=breaks),
##'               z4=cut(y4,breaks=breaks))
##' 
##' 
##' multinomial(d[,5])
##' a1 <- multinomial(d[,5:6])
##' (K1 <- kappa(a1))
##' ## irr::kappa(d[,5:6])
##' 
##' K2 <- kappa(multinomial(d[,7:8]))
##' crossprod(iid(K1)-iid(K2))^.5
##' sqrt(vcov(K1)+vcov(K2)) ## Wrong
##' 
##' @author Klaus Kähler Holst
multinomial <- function(x,...) {
    if (NCOL(x)==1) {
        x <- as.factor(x)
        lev <- levels(x)
        p <- length(lev)
        n <- length(x)
        P <- table(x)/n
        iid <- matrix(0,n,p)        
        for (i in seq(p)) {
            iid[,i] <- (1*(x==lev[i])-P[i])/n
        };
        res <- list(coef=P,vcov=crossprod(iid),iid=iid,position=seq(p),levels=lev,data=x)
        class(res) <- "multinomial"
        return(res)
    }
    if (NCOL(x)!=2L) stop("Matrix or data.frame with two columns expected")
    x <- as.data.frame(x)
    x[,1] <- as.factor(x[,1])
    x[,2] <- as.factor(x[,2])
    lev1 <- levels(x[,1])
    lev2 <- levels(x[,1])
    if (!identical(lev1,lev2)) stop("Levels should be identical")
    p <- length(lev1)
    M <- table(x)
    n <- sum(M)
    Pos <- P <- M/n
    iid <- matrix(0,n,p^2)
    for (i in seq(p)) {
        for (j in seq(p)) {
            pos <- (i-1)*p+j
            iid[,pos] <- (x[,1]==lev1[j])*(x[,2]==lev1[i])-P[j,i]
            Pos[j,i] <- pos
        }
    }; iid <- iid/n
    res <- list(coef=P,vcov=crossprod(iid),iid=iid, position=Pos, call=match.call(), levels=lev1, data=x)
    class(res) <- "multinomial"
    res
}

##' @S3method model.frame multinomial
model.frame.multinomial <- function(formula,...) {
    formula$data
}

##' @S3method iid multinomial
iid.multinomial <- function(x,...) {
    x$iid
}

##' @S3method coef multinomial
coef.multinomial <- function(object,...) {
    as.vector(object$coef)
}

##' @S3method vcov multinomial
vcov.multinomial <- function(object,...) {
    object$vcov
}

##' @S3method print multinomial
print.multinomial <- function(x,...) {
    cat("Call: "); print(x$call)
    cat("\nEstimates:\n")
    print(x$coef,quote=FALSE)
    stderr <- diag(vcov(x))^.5
    StdErr <- x$position
    StdErr[] <- stderr[StdErr]
    cat("\nStd.Err:\n")
    print(StdErr,quote=FALSE)
    ## cat("\nPosition:\n")
    ## print(x$position,quote=FALSE)
}

##' @S3method kappa multinomial
kappa.multinomial <- function(z,...) {
    pp <- length(coef(z))
    k <- length(z$levels)
    zeros <- rbind(rep(0,pp))
    A0 <- zeros; A0[diag(z$position)] <- 1
    A <- matrix(0,ncol=pp,nrow=2*k)
    for (i in seq(k)) A[i,z$position[i,]] <- 1
    for (i in seq(k)) A[i+k,z$position[,i]] <- 1    
    b <- estimate(z,function(p) as.vector(rbind(A0,A)%*%p),iid=TRUE)    
    b2 <- estimate(b,function(p) c(p[1],sum(p[seq(k)+1]*p[seq(k)+k+1])),iid=TRUE)
    estimate(b2,function(p) (p[1]-p[2])/(1-p[2]),iid=TRUE)
}
