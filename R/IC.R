
##' Extract i.i.d. decomposition (influence function) from model object
##'
##' Extract i.i.d. decomposition (influence function) from model object
##' @export
##' @usage
##'
##' IC(x,...)
##'
##' \method{IC}{default}(x, bread, id=NULL, folds=0, maxsize=(folds>0)*1e6,...)
##'
##' @aliases IC.default var_ic
##' @param x model object
##' @param id (optional) id/cluster variable
##' @param bread (optional) Inverse of derivative of mean score function
##' @param folds (optional) Calculate aggregated iid decomposition (0:=disabled)
##' @param maxsize (optional) Data is split in groups of size up to 'maxsize' (0:=disabled)
##' @param ... additional arguments
##' @examples
##' m <- lvm(y~x+z)
##' distribution(m, ~y+z) <- binomial.lvm("logit")
##' d <- sim(m,1e3)
##' g <- glm(y~x+z,data=d,family=binomial)
##' var_ic(IC(g))
##'
IC <- function(x,...) UseMethod("IC")

##' @export
influence.estimate <- function(model, ...)
  IC(model, ...)

##' @export
IC.default <- function(x, bread, id=NULL,
                        folds=0, maxsize=(folds>0)*1e6, ...) {

    if (any(paste("iid",class(x),sep=".") %in% methods("iid"))) {
      ## 'iid' method exists for the specific class.
      ## This is a scaled version of the influence function, hence
      ## we need to rescale.
      cl <- match.call()
      cl[[1]] <- substitute(iid)
      ii <- eval.parent(cl)
      if (!is.null(attr(ii, "bread"))) {
        attr(res, "bread") <- attr(res, "bread")*NROW(res)
      }
      ii <- ii*NROW(ii)
      return(ii)
    }
    if (!any(paste("score",class(x),sep=".") %in% methods("score"))) {
        warning("Not available for this class")
        return(NULL)
    }

    if (folds>0 || maxsize>0 || (!missing(id) && lava.options()$cluster.index)) {
        if (!requireNamespace("mets",quietly=TRUE)) stop("Requires 'mets'")
    }

    if (folds>0) {
        U <- Reduce("rbind",mets::divide.conquer(function(data) score(x,data=data,...),
                                                 id=id,
                                                 data=data,size=round(nrow(data)/folds)))
    } else {
        U <- score(x,indiv=TRUE,...)
    }
    pp <- pars(x)
    if (!missing(bread) && is.null(bread)) {
      bread <- suppressWarnings(vcov(x)*NROW(U))
    }
    if (missing(bread)) bread <- attributes(U)$bread
    if (is.null(bread)) {
        bread <- attributes(x)$bread
        if (is.null(bread)) bread <- x$bread
        if (is.null(bread)) {
            if (maxsize>0) {
                ff <- function(p) colSums(Reduce("rbind",mets::divide.conquer(function(data) score(x,data=data,p=p,...),
                                                                       data=data,size=maxsize)))
                I <- -numDeriv::jacobian(ff,pp,method=lava.options()$Dmethod)
            } else {
                I <- -numDeriv::jacobian(function(p) score(x,p=p,indiv=FALSE,...),pp,method=lava.options()$Dmethod)
            }
            bread <- Inverse(I)*NROW(U)
        }
    }
    ic0 <- U%*%bread
    if (!missing(id)) {
        N <- nrow(ic0)
        if (!lava.options()$cluster.index) {
            ic0 <- matrix(unlist(by(ic0,id,colSums)),byrow=TRUE,ncol=ncol(bread))
        } else {
            ic0 <- mets::cluster.index(id,mat=ic0,return.all=FALSE)
        }
        attributes(ic0)$N <- N
    }
    colnames(ic0) <- colnames(U)
    return(structure(ic0,bread=bread))
}


##' @export
IC.multigroupfit <- function(x,...) IC.default(x, combine=TRUE, ...)

##' @export
IC.matrix <- function(x,...) {
    p <- NCOL(x)
    n <- NROW(x)
    mu <- colMeans(x,na.rm=TRUE); S <- var(x,use="pairwise.complete.obs")*(n-1)/n
    ic1 <- t(t(x)-mu)
    ic2 <- matrix(ncol=(p+1)*p/2,nrow=n)
    pos <- 0
    nn <- c()
    cc <- mu
    for (i in seq(p))
        for (j in seq(i,p)) {
            pos <- pos+1
            cc <- c(cc,S[i,j])
            ic2[,pos] <- (ic1[,i]*ic1[,j])-cc[length(cc)]
            nn <- c(nn,paste(colnames(x)[c(i,j)],collapse=lava.options()$symbols[2]))
        }
    colnames(ic1) <- colnames(x); colnames(ic2) <- nn
    names(cc) <- c(colnames(ic1),colnames(ic2))
    structure(cbind(ic1,ic2),
              coef=cc,
              mean=mu, var=S)
}

##' @export
IC.numeric <- function(x,...) {
    n <- length(x)
    mu <- mean(x); S <- var(x)*(n-1)/n
    ic1 <- t(t(x)-mu)
    structure(cbind(mean=ic1, var=(ic1^2-S)),
              coef=c(mean=mu, var=S), mean=mu, var=S)
}

##' @export
IC.data.frame <- function(x,...) {
    if (!all(apply(x[1,,drop=FALSE],2,function(x) inherits(x,c("numeric","integer")))))
        stop("Don't know how to handle data.frames of this type")
    IC(as.matrix(x))
}
