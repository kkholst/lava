##' Estimate probabilities in contingency table
##'
##' @title Estimate probabilities in contingency table
##' @aliases multinomial kappa.multinomial kappa.table gkgamma
##' @param x Matrix or data.frame with observations (1 or 2 columns)
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
##' 
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
##' @author Klaus K. Holst
multinomial <- function(x,...) {
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
        coefs <- as.vector(P); names(coefs) <- paste("p",seq(k))
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
    names(coefs) <-  as.vector(outer(seq(k1),seq(k2),function(...) paste("p",...,sep="")))
    res <- list(coef=coefs,P=P,vcov=crossprod(iid),iid=iid, position=Pos, call=match.call(), levels=list(lev1,lev2), data=x)
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
    object$coef
}

##' @S3method vcov multinomial
vcov.multinomial <- function(object,...) {
    object$vcov
}

##' @S3method print multinomial
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

##' @S3method kappa multinomial
kappa.multinomial <- function(z,...) {    
    pp <- length(coef(z))
    if ((length(z$levels)!=2) || !(identical(z$levels[[1]],z$levels[[2]])))
        stop("Expected square table and same factor levels in rows and columns")
    k <- length(z$levels[[1]])
    zeros <- rbind(rep(0,pp))
    A0 <- zeros; A0[diag(z$position)] <- 1
    A <- matrix(0,ncol=pp,nrow=2*k)
    for (i in seq(k)) A[i,z$position[i,]] <- 1
    for (i in seq(k)) A[i+k,z$position[,i]] <- 1    
    b <- estimate(z,function(p) as.vector(rbind(A0,A)%*%p),iid=TRUE)    
    b2 <- estimate(b,function(p) c(p[1],sum(p[seq(k)+1]*p[seq(k)+k+1])),iid=TRUE)
    estimate(b2,function(p) list(kappa=(p[1]-p[2])/(1-p[2])),iid=TRUE)
}

##' @S3method kappa table
kappa.table <- function(z,...) {
    kappa(multinomial(Expand(z)),...)
}

##' @S3method kappa data.frame
kappa.data.frame <- function(z,...) {
    kappa(multinomial(z),...)
}

goodmankruskal_gamma <- function(P,...) {
    nr <- nrow(P); nc <- ncol(P)
    Pconc <- 0
    for (i in seq_len(nr-1)) {
        h <- seq(i+1,nr)
        for (j in seq_len(nc-1)) {
                k <- seq(j+1,nc)
                Pconc <- Pconc+2*P[i,j]*sum(P[h,k])
            }
    }
    Pdisc <- 0
    for (i in seq_len(nr-1)) {
        h <- seq(i+1,nr)
        for (j in (seq_len(nc-1)+1)) {
            k <- seq(1,j-1)
            Pdisc <- Pdisc+2*P[i,j]*sum(P[h,k])
        }
    }
    (Pconc-Pdisc)/(Pconc+Pdisc)    
}

##' @export
gkgamma <- function(x,...) {
    if (is.table(x) || is.data.frame(x) || is.matrix(x)) {
        x <- multinomial(x)
    }
    if (!inherits(x,"multinomial")) stop("Expected table, data.frame or multinomial object")

    P <- x$position
    estimate(x,function(p) { P[] <- p[P]
                             list(gamma=goodmankruskal_gamma(P)) },iid=TRUE)    
}




independence <- function(x,...) {
    if (is.table(x) || is.data.frame(x) || is.matrix(x)) {
        x <- multinomial(x)
    }
    if (!inherits(x,"multinomial")) stop("Expected table, data.frame or multinomial object")
    if (length(x$levels)!=2) stop("Data from two categorical variables expected")    
    f <- function(p) {
        P <- x$position; P[] <- p[x$position]
        n <- nrow(x$iid)
        k1 <- length(x$levels[[1]])
        k2 <- length(x$levels[[2]])        
        A1 <- matrix(0,ncol=length(p),nrow=k1)
        for (i in seq(k1)) A1[i,x$position[i,]] <- 1
        A2 <- matrix(0,ncol=length(p),nrow=k2)        
        for (i in seq(k2)) A2[i,x$position[,i]] <- 1
        P1 <- A1%*%p
        P2 <- A2%*%p
        I <- P1%*%t(P2)
        Q <- P-I
#        Q <- sum(n*P*(log(I[1,1])-P1
        sum((P-I)^2)
        ##V <- sqrt(sum((P*n-I*n)^2/I/n) /(n*(min(k1,k2)-1)))        
        V <- sqrt(sum((P-I)^2/I)   / ((min(k1,k2)-1)))
        return(V)
        sum(n*Q^2/I)^0.25
        return((sum((P-I)^2))^.5)
##        V
    }

    ## M <- P*n
    ## O2 <- colSums(M)
    ## O1 <- rowSums(M)
    ## M[1,1]-O1[1]*O2[1]/200
    ## M[2,2]-O1[2]*O2[2]/200
    ## sum((M-I*n)^2/(I*n))
    ## sum((P*n-I*n)^2/I/n)   
    ## sum(Q)
    ## sum(Q^2)
    ## M <- P
    ## chisq.test(M,correct=FALSE)

    return(estimate(x,function(p) list(cramersV=f(p)),iid=TRUE))
    
    
    e <- estimate(x,f,iid=TRUE,print=function(x,...) {
        cat("\tTest for independence\n\n")
        cat("Test statistc:\t ", formatC(x$coefmat[1]/x$coefmat[2]),
            "\nP-value:\t ", x$coefmat[5],"\n\n")
        print(estimate(x))
    })
    return(list(p.value=e$coefmat[5]))
    ## Q <- sum((a$coefmat[1,1]/a$coefmat[1,2]))
    ## df <- nrow(a$coefmat)
    ## res <- list(##data.name=hypothesis,
    ##             statistic = Q, parameter = df,
    ##             p.value=1-pchisq(Q,df=1),
    ##             method = "Test for independence")
    ## class(res) <- "htest"
    ## res
}

## independence(x)
## chisq.test(table(dd))



