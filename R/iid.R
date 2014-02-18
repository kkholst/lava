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
  U <- score(x,indiv=TRUE,combine=TRUE,...)
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
      if (!lava.options()$cluster.index) {
          iid0 <- matrix(unlist(by(iid0,id,colSums)),byrow=TRUE,ncol=ncol(iI))
      } else {
          iid0 <- mets::cluster.index(id,mat=iid0,return.all=FALSE)
      }
  }
  return(structure(iid0,iI=iI))
}

