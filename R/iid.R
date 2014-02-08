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
    if (is.null(score.deriv)) {
      score.deriv <- -numDeriv::jacobian(function(p) score(x,p=p,...),pp)
    }
    if (is.function(score.deriv)) {
      score.deriv <- score.deriv(x,p=pp,...)
    }
    if (is.matrix(score.deriv)) {
      iI <- Inverse(score.deriv)
    }
  }
  iid0 <- U%*%iI
  if (!missing(id)) {
      iid0 <- matrix(unlist(by(iid0,id,colSums)),byrow=TRUE,ncol=ncol(iI))
  }
  return(structure(iid0,iI=iI))
}

##' @S3method iid cox.aalen
iid.cox.aalen <- function(x,time.idx,...) {
    if (missing(time.idx)) return(x$gamma.iid)
    if (!all(time.idx)%in%seq(NCOL(x$cum)-1)) stop("Wrong index")
    if (time.idx==0) return(x$B.iid)
    res <- lapply(time.idx, 
                  function(i) Reduce("rbind",lapply(x$B.iid,function(x) x[,i])))
    if (length(idx)==1) return(res[[1]])
    return(res)
}
