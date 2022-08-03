normal.threshold <- function(object,p=coef(object),...) {
    M <- moments(object,p=p)
    ord <- ordinal(Model(object))
    K <- attributes(ord)$K
    cK <- c(0,cumsum(K-1))
    breaks.orig <- list()
    for (i in seq(K)) {
        breaks.orig <- c(breaks.orig,list(M$e[seq(K[i]-1)+cK[i]]))
    }
    breaks <- lapply(breaks.orig, ordreg_threshold)
    names(breaks) <- names(K)
    ii <- match(names(K),vars(object))
    sigma <- M$Cfull[ii,ii]
    list(breaks=breaks,sigma=sigma,mean=M$v[ii],K=K)
}

prob.normal <- function(sigma,breaks,breaks2=breaks) {
    if (!requireNamespace("mets",quietly=TRUE)) stop("'mets' package required")
    if (ncol(sigma)!=2 || missing(breaks)) stop("Wrong input")
    P <- matrix(ncol=length(breaks2)-1, nrow=length(breaks)-1)
    for (i in seq(length(breaks)-1))
        for (j in seq(length(breaks2)-1))
            P[i,j] <- mets::pmvn(lower=c(breaks[i],breaks2[j]),upper=c(breaks[i+1],breaks2[j+1]),sigma=sigma)
    return(P)
}

assoc <- function(P,sigma,breaks,...) {
    if (missing(P)) P <- prob.normal(sigma,breaks,...)
    Agree <- sum(diag(P))
    marg.row <- rowSums(P)
    marg.col <- colSums(P)
    Chance <- sum(marg.row*marg.col)
    kap <- (Agree-Chance)/(1-Chance)
    gam <- goodmankruskal_gamma(P)$gamma
    inf <- information_assoc(P)
    res <- c(list(kappa=kap,gamma=gam),inf)
    if (!missing(sigma)) res <- c(res,rho=sigma[1,2])
    return(res)
}



##################################################
### Risk comparison
##################################################

## or:= riskcomp(x,scale=odds) // OR
##' @export
riskcomp <- function(x,...,scale,op="/",type=1,struct=FALSE) {
    val <- c(x,unlist(list(...)))
    if (!missing(scale)) val <- do.call(scale,list(val))
    if (!struct && length(val)==2) {
        if (type==2) {
            return(do.call(op,list(val[2],val[1])))
        } else if (type==1) {
            return(do.call(op,list(val[1],val[2])))
        }
        return(c(do.call(op,list(val[2],val[1])),
                 do.call(op,list(val[1],val[2]))))
    }
    outer(val,val,op)
    offdiag(outer(val,val,op) ,type=type)
}

##' @export
Ratio <- function(x,...) riskcomp(x,...,op="/")

##' @export
Diff <- function(x,...) riskcomp(x,...,op="-")


##################################################
## Odds ratio
##################################################

##' @export
odds <- function(x) x/(1-x)

logor <- function(x) {
    c(log(prod(diag(x))/prod(revdiag(x))),sum(1/x)^.5)
}


##' @export
OR <- function(x,tabulate=FALSE,log=FALSE,...) {
    if (!inherits(x,c("multinomial","table"))) {
        val <- riskcomp(x,...,scale=odds)
        if (log) val <- base::log(val)
        return(val)
    }
    if (inherits(x,"multinomial")) {
        M <- x
    } else {
        M <- multinomial(x)
    }
    pos <- M$position
    if (ncol(pos)!=2 & ncol(pos)!=2) stop("Only for 2x2 tables")
    orfun <- function(p,...) {
        list(logOR=sum(log(p[diag(pos)]))-sum(log(p[revdiag(pos)])))
    }
    estimate(M,orfun,back.transform=exp)
}



##################################################
## Information theoretical measures
##################################################


information_assoc <- function(P,base=exp(1),...) {
    P.row <- rowSums(P)
    P.col <- colSums(P)
    H.row <- H.col <- H <- 0
    for (j in seq_along(P.col))
        if (P.col[j]>0) H.col <- H.col - P.col[j]*log(P.col[j]+(P.col[j]==0),base=base)
    for (i in seq_along(P.row)) {
        if (P.row[i]>0) H.row <- H.row - P.row[i]*log(P.row[i]+(P.row[i]==0),base=base)
        for (j in seq_along(P.col)) {
            if (P[i,j]>0) H <- H - P[i,j]*log(P[i,j],base=base)
        }
    }
    I <- H.row+H.col-H
    return(list(MI=I,H=H,H.row=H.row,H.col=H.col,
                U.row=I/H.row,U.col=I/H.col,U.sym=2*I/(H.row+H.col)))
}


##' @export
information.data.frame <- function(x,...) {
    information(multinomial(x,marginal=TRUE),...)
}

##' @export
information.table <- function(x,...) {
    information(multinomial(x,marginal=TRUE),...)
}

##' @export
information.multinomial <- function(x,...) {
    estimate(x,function(p,object,...) {
        P <- object$position; P[] <- p[object$position]
        information_assoc(P)},...)
}


##################################################
## Independence tests
##################################################

independence <- function(x,...) {
    if (is.table(x) || is.data.frame(x) || is.matrix(x)) {
        x <- multinomial(x)
    }
    if (!inherits(x,"multinomial")) stop("Expected table, data.frame or multinomial object")
    if (length(x$levels)!=2) stop("Data from two categorical variables expected")
    f <- function(p) {
        P <- x$position; P[] <- p[x$position]
        n <- nrow(x$IC)
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

    return(estimate(x,function(p) list(cramersV=f(p)),IC=TRUE,...))


    e <- estimate(x,f,IC=TRUE,print=function(x,...) {
        cat("\tTest for independence\n\n")
        cat("Test statistc:\t ", formatC(x$coefmat[1]/x$coefmat[2]),
            "\nP-value:\t ", x$coefmat[5],"\n\n")
        print(estimate(x))
    },...)
    return(list(p.value=e$coefmat[5]))
    ## Q <- sum((a$coefmat[1,1]/a$coefmat[1,2]))
    ## df <- nrow(a$coefmat)
    ## res <- list(##data.name=hypothesis,
    ##             statistic = Q, parameter = df,
    ##             p.value=pchisq(Q,df=1,lower.tail=FALSE),
    ##             method = "Test for independence")
    ## class(res) <- "htest"
    ## res
}

## independence(x)
## chisq.test(table(dd))
