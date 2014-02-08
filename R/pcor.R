pcor <- function(x,y,X,Z,start,...) {
    if (is.numeric(x) && is.numeric(y)) {
        e <- estimate(covariance(lvm(),x~y))
        return(estimate(e,function(p) list(rho=p[5]/(p[3]*p[4])^.5),iid=TRUE))
    }
    
    n1 <- 1+seq(nlevels(x)-1)
    n2 <- n1[length(n1)]+seq(nlevels(y)-1)
    browser()
    if (missing(start)) {
        f <- as.formula(ifelse(missing(X),"~1","~X"))
        start <- c(0.5,
                   attr(lava:::ordreg(update(f,x~.),fast=TRUE,family=binomial("probit")),"threshold"),
                   attr(lava:::ordreg(update(f,y~.),fast=TRUE,family=binomial("probit")),"threshold"))
    } 
    browser()
    nn <- table(x,y)
    ff <- function(theta) {
        -sum(as.vector(nn)*log(polycor0(theta[1],theta[n1],theta[n2])))
    }
    gg <- function(theta) {
        pp <- polycor0(theta[1],theta[n1],theta[n2],onlyP=FALSE)
        np <- as.vector(nn)/as.vector(pp$p)
        -colSums(apply(pp$dp,2,function(x) np*x))
    }

    nn0 <- nn; nn[nn==0] <- .5
    p0 <- as.vector(nn)/sum(nn)
    logL0 <- sum(as.vector(nn)*log(p0))
    t0 <- system.time(op <- nlminb(start,ff,gg))
    cc <- op$par
    names(cc) <- c("rho",paste(rownames(nn),"x",sep=".")[-1], paste(colnames(nn),"y",sep=".")[-1])
    V <- solve(numDeriv::jacobian(function(p) gg(p), op$par))
    
    res <- list(coef=cc, vcov=V, tab=nn, logLik0=logL0, logLik=-ff(op$par), n1=n1, n2=n2)
    structure(res,class="pcor")
}

##' @S3method coef pcor
coef.pcor <- function(object,...) object$coef

##' @S3method vcov pcor
vcov.pcor <- function(object,...) object$vcov

##' @S3method logLik pcor
logLik.pcor <- function(object,p=coef(object),...) {
    u <- polycor0(p[1],p[object$n1],p[object$n2],onlyP=TRUE)
    np <- sum(as.vector(object$tab)*log(as.vector(u)))
    nobs <- sum(ps$tab)/2
    structure(np,nall=nobs,nobs=nobs,df=length(p),class="logLik")
}

##' @S3method print pcor
print.pcor <- function(x,...) {
    res <- cbind(coef(x),diag(vcov(x))^0.5)
    colnames(res) <- c("Estimate","Std.Err")
    print(res)
    df <- length(x$tab)-nrow(res)
    q <- with(x,2*(logLik0-logLik))
    cat("\nDeviance = ", q, ", df = ",df,"\n")
}

##' @S3method score pcor
score.pcor <- function(x,p=coef(x),indiv=FALSE,...) {
    u <- polycor0(p[1],p[x$n1],p[x$n2],onlyP=FALSE)
    if (!indiv) {
        np <- as.vector(x$tab)/as.vector(u$p)
        return(colSums(apply(u$dp,2,function(x) np*x)))
    }
    U <- u$dp;
    U <- apply(u$dp,2,function(x) x/as.vector(u$p))
    ii <- unlist(apply(cbind(seq(length(x$tab)),as.vector(x$tab)),1,function(x) rep(x[1],x[2])))    
    return(U[ii,])
}
    


polycor0 <- function(rho,a0,b0,onlyP=TRUE,...) {
    k1 <- length(a0); k2 <- length(b0)
    S <- diag(c(1-rho,1-rho))+rho
    P <- matrix(0,nrow=k1,ncol=k2)
    P1 <- pnorm(a0,sd=1)
    P2 <- pnorm(b0,sd=1)
    for (i in seq(k1))
        for (j in seq(k2)) P[i,j] <- pmvnorm(lower=c(-Inf,-Inf),upper=c(a0[i],b0[j]),sigma=S)

    PP <- Drho <- matrix(0,nrow=k1+1,ncol=k2+1)
    pmvn <- function(i,j,sigma=S) {
        if (i==0 | j==0) return(0)
        if (i==(k1+1) & j==(k2+1)) return(1)
        if (i==(k1+1)) return(P2[j])
        if (j==(k2+1)) return(P1[i])
        P[i,j]
    }

    dpmvn <- function(i,j,type=1,k) {
        if (i==0 | j==0) return(0)
        if (i==(k1+1) & j==(k2+1)) return(0)
        if (i==(k1+1)) {
            if (type==3 && k==j) return(dnorm(b0[j]))
                return(0) 
        }
        if (j==(k2+1)) {
            if (type==2 && k==i) return(dnorm(a0[i]))
            return(0)
        }
        if (type==1) ## rho
            return(dmvnorm(c(a0[i],b0[j]),sigma=S))
        if (type==2) { ## threshold a
            if (k!=i) return(0)
            return(dnorm(a0[i])*pnorm((b0[j]-rho*a0[i])/sqrt(1-rho^2)))
        } ## threshold b
        if (k!=j) return(0)
        dnorm(b0[j])*pnorm((a0[i]-rho*b0[j])/sqrt(1-rho^2))
    }
    
    for (i in seq(k1+1))
        for (j in seq(k2+1)) {
            PP[i,j] <- pmvn(i,j) + pmvn(i-1,j-1) -
                pmvn(i-1,j) - pmvn(i,j-1)
            Drho[i,j] <- dpmvn(i,j) + dpmvn(i-1,j-1) -
                dpmvn(i-1,j) - dpmvn(i,j-1)
        }
    if (onlyP) return(PP)
    
    Da <- matrix(0,length(PP),k1)
    for (k in seq(k1))
        for (i in seq(k1+1))
            for (j in seq(k2+1)) {
                pos <- i + (k1+1)*(j-1)
                Da[pos,k] <- dpmvn(i,j,type=2,k=k) + dpmvn(i-1,j-1,type=2,k=k) -
                    dpmvn(i-1,j,type=2,k=k) - dpmvn(i,j-1,type=2,k=k)
            }
    
    Db <- matrix(0,length(PP),k2)
    for (k in seq(k2))
        for (i in seq(k1+1))
            for (j in seq(k2+1)) {
                pos <- i + (k1+1)*(j-1)
                Db[pos,k] <- dpmvn(i,j,type=3,k=k) + dpmvn(i-1,j-1,type=3,k=k) -
                    dpmvn(i-1,j,type=3,k=k) - dpmvn(i,j-1,type=3,k=k)
            }
    
    list(p=PP,dp=cbind(as.vector(Drho),Da,Db))
    
}

