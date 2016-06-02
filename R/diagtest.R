##' Calculate prevalence, sensitivity, specificity, and positive and
##' negative predictive values
##'
##' @title Calculate diagnostic tests for 2x2 table
##' @param table Table or (matrix/data.frame with two columns)
##' @param positive Switch reference
##' @param exact If TRUE exact binomial proportions CI/test will be used
##' @param p0 Optional null hypothesis (test prevalenc, sensitivity, ...)
##' @param confint Type of confidence limits
##' @param ... Additional arguments to lower level functions
##' @author Klaus Holst
##' @details Table should be in the format with outcome in columns and
##'     test in rows.  Data.frame should be with test in the first
##'     column and outcome in the second column.
##' @examples
##' M <- as.table(matrix(c(42,12,
##'                        35,28),ncol=2,byrow=TRUE,
##'                      dimnames=list(rater=c("no","yes"),gold=c("no","yes"))))
##' diagtest(M,exact=TRUE)
##' @export
diagtest <- function(table,positive=2,exact=FALSE,p0=NA,confint=c("logit","arcsin","pseudoscore","exact"),...) {
    if (!inherits(table,c("table","data.frame","matrix","multinomial")))
        stop("Expecting a table or data.frame.")
    if (is.table(table)) {
        lev <- dimnames(table)[[2]]
        }
    if (inherits(table,"multinomial")) {
        lev <- dimnames(table$P)[[2]]
    }
    if (!is.table(table) & (is.matrix(table) || is.data.frame(table))) {
        if (is.factor(table[,2])) {
            lev <- levels(table[,2])
        } else
            lev <- unique(table[,2])
    }
    if (is.character(positive)) {
        positive <- match(positive,lev)
    }
    if (!(positive%in%c(1,2))) stop("Expecting and index of 1 or 2.")
    negative <- positive%%2+1L
    if (!is.null(confint) && confint[1]=="exact") exact <- TRUE
    if (exact) {
        if (!is.table(table) && (is.matrix(table) || is.data.frame(table))) {
            table <- base::table(table[,c(1,2),drop=FALSE])
            ##names(dimnames(table)) <- colnames(table)[1:2]
        }
        if (!is.table(table) || nrow(table)!=2 || ncol(table)!=2) stop("2x2 table expected")
        n <- sum(table)
        nc <- colSums(table)
        nr <- rowSums(table)
        test <- TRUE
        if (is.na(p0)) {
            test <- FALSE
            p0 <- 0.5
        }
        ## Prevalence
        p1 <- with(stats::binom.test(nc[positive],n,p=p0),c(estimate,conf.int,p.value))
        ## Test marginal
        p2 <- with(stats::binom.test(nr[positive],n,p=p0),c(estimate,conf.int,p.value))
        ## Sensitivity/Specificity
        sens <- with(stats::binom.test(table[positive,positive],nc[positive],p=p0),c(estimate,conf.int,p.value))
        spec <- with(stats::binom.test(table[negative,negative],nc[negative],p=p0),c(estimate,conf.int,p.value))
        ## PPV,NPV
        ppv <- with(stats::binom.test(table[positive,positive],nr[positive],p=p0),c(estimate,conf.int,p.value))
        npv <- with(stats::binom.test(table[negative,negative],nr[negative],p=p0),c(estimate,conf.int,p.value))
        ## Accuracy
        acc <- with(stats::binom.test(table[positive,positive]+table[negative,negative],n,p=p0),c(estimate,conf.int,p.value))
        ## Symmetry (McNemar):
        ##   number of discordant pairs under null: b~bin(b+c,0.5)
        sym <- with(stats::binom.test(table[positive,negative],table[positive,negative]+table[negative,positive],p=0.5),c(estimate,conf.int,p.value))
        coefmat <- rbind(Prevalence=p1,
                         Test=p2,
                         Sensitivity=sens,
                         Specificity=spec,
                         PositivePredictiveValue=ppv,
                         NegativePredictiveValue=npv,
                         Accuracy=acc,
                         Homogeneity=sym)
        if (!test) coefmat[seq(nrow(coefmat)-1),4] <- NA
        coefmat <- cbind(coefmat[,1,drop=FALSE],NA,coefmat[,-1,drop=FALSE])
        colnames(coefmat) <- c("Estimate","Std.Err","2.5%","97.5%","P-value")
        res <- list(table=table, prop.table=table/sum(table),
                    coefmat=coefmat)
    } else {
        if (inherits(table,"table"))
            M <- multinomial(M)
        else {
            if (inherits(table,"multinomial")) {
            } else {
                M <- multinomial(table[,1:2],...)
                table <- base::table(table)
            }
        }
        calc_diag <- function(p,...) {
            P <- matrix(p[1:4],2)
            p1 <- sum(P[,positive])
            p2 <- sum(P[positive,])
            res <- c(Prevalence=p1,  ##(p[1]+p[2]),
                     Test=p2,        ##(p[1]+p[3]),
                     Sensitivity=P[positive,positive]/p1,     ## p[1]/(p[1]+p[2]), # Prob test + | given (true) disease (True positive rate)
                     Specificity=P[negative,negative]/(1-p1), ## p[4]/(1-p[1]-p[2]), # Prob test - | given no disease (True negative rate)
                     PositivePredictiveValue=P[positive,positive]/p2,     ## p[1]/(p[1]+p[3]), # Prob disease | test +
                     NegativePredictiveValue=P[negative,negative]/(1-p2), ## p[4]/(1-p[1]-p[3]), # Prob disease free | test -
                     Accuracy=(P[1,1]+P[2,2])/sum(P),
                     Homogeneity=P[negative,positive]-P[positive,negative]
                     )
            if (!is.null(confint)) {
                if (tolower(confint[1])=="logit") {
                    res[seq(length(res)-1)] <- logit(res[seq(length(res)-1)])
                } else {
                    res[seq(length(res)-1)] <- asin(sqrt(res[seq(length(res)-1)]))
                }
            }
            return(res)
        }

        names(dimnames(table)) <- paste0(c("Test:","Outcome:"),names(dimnames(table)))
        prfun <- function(x,...) {
            printCoefmat(x$coefmat[,c(-2)],na.print="",...)
            printline()
            cat("\n")
            cat("Prevalence:				Prob( outcome+ )\n")
            cat("Test:					Prob( test+ )\n")
            cat("Sensitivity (True positive rate):	Prob( test+ | outcome+ )\n")
            cat("Specificity (True negative rate):	Prob( test- | outcome- )\n")
            cat("Positive predictive value (Precision):	Prob( outcome+ | test+ )\n")
            cat("Negative predictive value:		Prob( outcome- | test- )\n")
            cat("Accuracy:				Prob( correct classification )\n")
            cat("Homogeneity/Symmetry:			Prob( outcome+ ) - Prob( test+ )\n")
        }

        btransform <- NULL
        if (!is.null(confint)) {
            if (tolower(confint[1])=="logit") {
                btransform <- function(x) {
                    rbind(expit(x[seq(nrow(x)-1),,drop=FALSE]),
                          x[nrow(x),,drop=FALSE])
                }
            } else if (tolower(confint[1])=="pseudoscore") {
                ## TODO, agresti-ryu, biometrika 2010
            } else if (tolower(confint[1])=="arcsin")  {
                btransform <- function(x) {
                    rbind(sin(x[seq(nrow(x)-1),,drop=FALSE])^2,
                          x[nrow(x),,drop=FALSE])
                }
            }
        }
        res <- estimate(M,calc_diag,print=prfun,null=c(rep(p0,7),0),transform.ci=btransform,...)
    }

    CI <- confint[1]
    if (exact) CI <- "exact"
    if (is.null(CI)) CI <- "wald"
    res <- structure(c(res,
                       list(table=table, prop.table=table/sum(table),
                            confint=CI,
                            positive=positive,
                            negative=negative,
                            levels=dimnames(table)
                            )),
                     class=c("diagtest","estimate"))
    res$call <- match.call()
    rownames(res$coefmat) <- gsub("\\[|\\]","",rownames(res$coefmat))
    names(res$coef) <- rownames(res$coefmat)
    return(res)
}

print.diagtest <- function(x,...) {
    cat("Call: "); print(x$call)
    cat("Confidence limits: ", x$confint,"\n",sep="")
    printline()
    printmany(x$table,x$prop.table,nspace=2,...)
    cat("\nPositive outcome: '", x$levels[[2]][x$positive],"'\n",sep="")
    ##cat("\tNegative outcome: '", x$levels[[2]][x$positive%%2+1],"'\n",sep="")
    printline()
    printCoefmat(x$coefmat[,c(-2)],na.print="",...)
    printline()
    cat("\n")
    cat("Prevalence:				Prob( outcome+ )\n")
    cat("Test:					Prob( test+ )\n")
    cat("Sensitivity (True positive rate):	Prob( test+ | outcome+ )\n")
    cat("Specificity (True negative rate):	Prob( test- | outcome- )\n")
    cat("Positive predictive value (Precision):	Prob( outcome+ | test+ )\n")
    cat("Negative predictive value:		Prob( outcome- | test- )\n")
    cat("Accuracy:				Prob( correct classification )\n")
    if (x$confint=="exact") {
        cat("Homogeneity/Symmetry:			Prob( outcome+, test- | discordant ), H0: p=0.5 \n")
    } else {
        cat("Homogeneity/Symmetry:			H0: Prob( outcome+ ) - Prob( test+ ), H0: p=0\n")
    }
    cat("\n")
}

summary.diagtest <- function(x,...) {
    x[c("iid","print","id","compare")] <- NULL
    return(x)
}
