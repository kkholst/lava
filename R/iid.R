##' Extract i.i.d. decomposition (influence function) from model object
##'
##' Extract i.i.d. decomposition (influence function) from model object 
##' @export
##' @usage
##' 
##' iid(x,...)
##' 
##' \method{iid}{default}(x,score.deriv,id,...)
##' 
##' @aliases iid.default
##' @param x model object
##' @param id id/cluster variable (optional)
##' @param score.deriv (optional) derivative of mean score function 
##' @param ... additional arguments
##' @examples
##' m <- lvm(y~x+z)
##' distribution(m, ~y+z) <- binomial.lvm("logit")
##' d <- sim(m,1e3)
##' g <- glm(y~x+z,data=d,family=binomial)
##' crossprod(iid(g))
##' 
iid <- function(x,...) UseMethod("iid")

##' @S3method iid default
iid.default <- function(x,score.deriv,id,...) {
    if (!any(paste("score",class(x),sep=".") %in% methods("score"))) {
        warning("Not available for this class")
        return(NULL)
    }    
    U <- score(x,indiv=TRUE,...)
    n <- NROW(U)
    pp <- pars(x)
    if (missing(score.deriv)) {
        iI <- vcov(x) 
    } else {
        if (!is.null(score.deriv) && is.function(score.deriv)) {
            score.deriv <- score.deriv(x,p=pp,...)
        } else if (!is.null(score.deriv) && is.matrix(score.deriv)) {
            iI <- Inverse(score.deriv)
        } else {
            score.deriv <- -numDeriv::jacobian(function(p) score(x,p=p,...),pp,method=lava.options()$Dmethod)    
        }
    }
    iid0 <- U%*%iI
    if (!missing(id)) {
        if (!lava.options()$cluster.index) {
            iid0 <- matrix(unlist(by(iid0,id,colSums)),byrow=TRUE,ncol=ncol(iI))
        } else {
            iid0 <- mets::cluster.index(id,mat=iid0,return.all=FALSE)
        }
    }
  return(structure(iid0,iI=iI))
}


##' @S3method iid matrix
iid.matrix <- function(x,...) {
    p <- ncol(x); n <- nrow(x)
    mu <- colMeans(x,na.rm=TRUE); S <- var(x,use="pairwise.complete.obs")*(n-1)/n
    iid1 <- t(t(x)-mu)
    iid2 <- matrix(ncol=(p+1)*p/2,nrow=n)
    pos <- 0
    nn <- c()
    cc <- mu
    for (i in seq(p))
        for (j in seq(i,p)) {
            pos <- pos+1
            cc <- c(cc,S[i,j])
            iid2[,pos] <- (iid1[,i]*iid1[,j])-cc[length(cc)]
            nn <- c(nn,paste(colnames(x)[c(i,j)],collapse=lava.options()$symbols[2]))
        }
    colnames(iid1) <- colnames(x); colnames(iid2) <- nn
    names(cc) <- c(colnames(iid1),colnames(iid2))
    iid1[is.na(iid1)] <- 0
    iid2[is.na(iid2)] <- 0    
    structure(cbind(iid1/n,iid2/n),
              coef=cc,
              mean=mu, var=S)
}

##' @S3method iid numeric
iid.numeric <- function(x,...) {
    n <- length(x)
    mu <- mean(x); S <- var(x)*(n-1)/n
    iid1 <- t(t(x)-mu)
    structure(cbind(mean=iid1/n,var=(iid1^2-S)/n),coef=c(mean=mu,var=S),mean=mu,var=S)
}


##' @S3method iid data.frame
iid.data.frame <- function(x,...) {
    if (!all(apply(x[1,,drop=FALSE],2,function(x) inherits(x,c("numeric","integer")))))
        stop("Don't know how to handle data.frames of this type")
    iid(as.matrix(x))
}




