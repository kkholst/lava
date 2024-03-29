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
    list(C=Pconc,D=Pdisc,gamma=(Pconc-Pdisc)/(Pconc+Pdisc))
}


##' @export
gkgamma <- function(x,data=parent.frame(),strata=NULL,all=FALSE,IC=TRUE,...) {
    if (inherits(x,"formula")) {
        xf <- getoutcome(x,sep="|")
        xx <- attr(xf,"x")
        if (length(xx)==0) stop("Not a valid formula")
        yx <- update(as.formula(paste0(xf,"~.")),xx[[1]])
        if (length(xx)>1) {
            strata <- interaction(model.frame(xx[[2]],data=data))
            x <- yx
        } else {
            x <- model.frame(yx,data=data)
        }
    }
    if (!is.null(strata)) {
        dd <- split(data,strata)
        gam <- lapply(dd,function(d,...) gkgamma(x,data=d,...),
                      ...,
                      IC=TRUE,
                      keep=1:2)
        mgam <- Reduce(function(x,y,...) merge(x,y,...),gam)
        ps <- estimate(multinomial(strata),data=data,...)
        mgam <- merge(mgam,ps)
        psi <- 2*length(gam)+seq(length(coef(ps)))
        res <- estimate(mgam,function(p,...) {
            k <- length(p)/3
            cd <- lapply(seq(k),function(x) p[(1:2)+2*(x-1)])
            dif <- unlist(lapply(cd,function(x) x[1]-x[2]))
            tot <- unlist(lapply(cd,function(x) x[1]+x[2]))
            gam <- dif/tot ## Conditional gammas given Z=z
            px2 <- p[psi]^2
            pgamma <- sum(dif*px2)/sum(tot*px2)
            #weights <- px2*tot/sum(px2*tot)
            #pgamma <- sum(weights*gam)
            c(gam,pgamma=pgamma)
        },labels=c(paste0("\u03b3:",names(dd)),"pgamma"),
        IC=IC)
        if (!IC) {
            for (i in seq_along(gam))
                gam[[i]][c("IC","id")] <- NULL
        }
        homtest <- estimate(res,lava::contr(seq_along(gam),length(gam)+1),IC=FALSE)
        attributes(res) <- c(attributes(res),
                             list(class=c("gkgamma","estimate"),
                                  cl=match.call(),
                                  strata=gam,
                                  homtest=homtest))
        return(res)
    }
    if (is.table(x) || is.data.frame(x) || is.matrix(x)) {
        x <- multinomial(x)
    }
    if (!inherits(x,"multinomial")) stop("Expected table, data.frame or multinomial object")
    structure(estimate(x,function(p) {
        P <- x$position; P[] <- p[x$position]
        goodmankruskal_gamma(P)
    },IC=IC,data=data,...),
    cl=match.call(),
    class=c("gkgamma","estimate"))
}

##' @export
print.gkgamma <- function(x,call=TRUE,...) {
    if (call) {
        cat("Call: ")
        print(attr(x,"cl"))
        print(cli::rule())
    }
    n <- x$n

    if (!is.null(attr(x,"strata"))) {
        cat("Strata:\n\n")
        for (i in seq_along(attr(x,"strata"))) {
            with(attributes(x), cat(paste0(names(strata)[i],
                                           " (n=",strata[[i]]$n,
                                           if (strata[[i]]$ncluster<strata[[i]]$n) paste0(",clusters=",strata[[i]]$ncluster),
                                           "):\n",sep="")))
            e <- attr(x,"strata")[[i]]
            print.estimate(e,type=0)
            cat("\n")
        }
        print(cli::rule())
        cat("\n")
        n <- sum(unlist(lapply(attr(x,"strata"),"[[","n")))
    }
    k <- x$ncluster
    if (!is.null(n) && !is.null(k) && k<n) {
        cat("n = ",n,", clusters = ",k,"\n\n",sep="")
    } else {
        if (!is.null(n)) {
            cat("n = ",n,"\n\n",sep="")
        } else if (!is.null(k)) {
            cat("n = ",k,"\n\n",sep="")
        }
    }
    if (!is.null(attr(x,"strata"))) {
        cat("Gamma coefficient:\n\n")
    }
    class(x) <- "estimate"
    print(x)
    ## if (!is.null(attr(x,"homtest"))) {
    ##     print(cli::rule())
    ##     cat("Homogeneity test:\n\n")
    ##     with(attr(x,"homtest")$compare,
    ##          cat("\u03c7\u00b2 = ",statistic,
    ##              ", df = ",parameter,
    ##              ", p-value = ",p.value,"\n",sep=""))
    ## }
    invisible(x)
}



