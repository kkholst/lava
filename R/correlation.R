##' Generic correlation method
##'
##' @title Generic method for extracting correlation coefficients of model object
##' @param x Object
##' @param \dots Additional arguments
##' @author Klaus K. Holst
##' @export
"correlation" <- function(x,...) UseMethod("correlation")

##' @export
correlation.lvmfit <- function(x,z=TRUE, IC=FALSE, back.transform=TRUE,...) {
    pp <- matrices2(Model(x), with(index(x),seq_len(npar+npar.mean+npar.ex)))$P
    pos <- pp[lower.tri(pp)][(index(x)$P0)[lower.tri(pp)]==1]
    if (length(pos)<1) return(NULL)
    pp0 <- pp
    pp0[index(x)$P0!=1] <- 0; pp0[lower.tri(pp0)] <- 0
    mynames <- vars(x)
    n <- nrow(pp0)
    ff <-  function(p) {
        res <- numeric(length(pos))
        nn <- character(length(pos))
        for (i in seq_along(pos)) {
            p0 <- pos[i]
            idx <- which(pp0==p0)
            rowpos <- (idx-1)%%n + 1
            colpos <- ceiling(idx/n)
            coefpos <- c(p0,pp0[rbind(c(rowpos,rowpos),c(colpos,colpos))])
            pval <- pp[rbind(c(rowpos,rowpos),c(colpos,colpos))]
            phi.v1.v2 <- numeric(3);
            newval <- p[coefpos]
            phi.v1.v2[coefpos!=0] <- newval
            phi.v1.v2[coefpos==0] <- pval[tail(coefpos==0,2)]
            rho <- atanh(phi.v1.v2[1]/sqrt(prod(phi.v1.v2[-1])))
            res[i] <- rho
            nn[i] <- paste(mynames[c(rowpos,colpos)],collapse="~")
        }
        structure(res,names=nn)
    }
    V <- NULL
    if (!IC) V <- vcov(x)
    if (back.transform) {
        back.transform <- tanh
    } else {
        back.transform <- NULL
    }
    estimate(x,coef=coef(x), vcov=V, f=ff, back.transform=back.transform,
             IC=IC, ...)
}


##' @export
correlation.matrix <- function(x,z=TRUE,back.transform=TRUE,mreg=FALSE,return.all=FALSE,...) {
    if (mreg) {
        m <- lvm()
        covariance(m,pairwise=TRUE) <- colnames(x)
        try(e <- estimate(m,as.data.frame(x),...),silent=TRUE)
        res <- correlation(e,...)
        if (return.all) {
            return(list(model=m,estimate=e,correlation=res))
        }
        return(res)
    }
    if (ncol(x)==2) {
        ii <- IC(x)
        ee <- estimate(coef=attributes(ii)$coef[3:5], IC=ii[,3:5])
        if (z) {
            if (back.transform) {
                ee <- estimate(ee, function(x) atanh(x[2]/sqrt(x[1]*x[3])), back.transform=tanh)
            } else {
                ee <- estimate(ee, function(x) atanh(x[2]/sqrt(x[1]*x[3])))
            }
        } else {
            ee <-  estimate(ee, function(x) x[2]/sqrt(x[1]*x[3]))
        }
        return(ee)
    }

    e <- c()
    R <- diag(nrow=ncol(x))
    dimnames(R) <- list(colnames(x),colnames(x))
    for (i in seq(ncol(x)-1))
        for (j in seq(i+1,ncol(x))) {
            e <- c(e,list(correlation(x[,c(i,j)],z=z,back.transform=FALSE,...)))
            R[j,i] <- coef(e[[length(e)]])
            if (z) R[j,i] <- tanh(R[j,i])
        }
    R <- R[-1,-ncol(R),drop=FALSE]
    res <- do.call(merge, c(e, paired=TRUE))
    if (z && back.transform) {
        res <- estimate(res,back.transform=tanh, print=function(x,digits=1,...) {
            print(x$coefmat[,-2,drop=FALSE],...)
            cat("\n")
            print(offdiag(R,type=4),digits=digits,...)
        })
    }
    return(res)
}

##' @export
correlation.data.frame <- function(x,...) {
    correlation(as.matrix(x),...)
}

