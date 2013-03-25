##' Estimation of functional of parameters 
##'
##' Estimation of functional of parameters 
##' @param x model object (\code{glm}, \code{lvmfit}, ...)
##' @param fun transformation of model parameters and optionally data.
##' @param data \code{data.frame}
##' @param id1 (optional) id-variable corresponding to iid decomposition of model parameters
##' @param id2 (optional) id-variable of data.frame
##' @param vcov If TRUE the inverse of \code{vcov(x)} is used to calculate ---
##' @param level Level of confidence limits
##' @param iid If TRUE the iid decompositions are returned instead of variance estimates
##' @param contrast (optional) Contrast matrix for final Wald test
##' @param null (optional) Null hypothesis to test 
##' @param ... additional arguments to lower level functions
##' @examples
##' ## Simulation from logistic regression model
##' m <- lvm(y~x+z);
##' distribution(m,y~x) <- binomial.lvm("logit")
##' d <- sim(m,1000)
##' g <- glm(y~z+x,data=d,family=binomial())
##'
##' ## Plain estimates (robust standard errors)
##' estimate(g)
##' ## Testing contrasts
##' estimate(g,null=0)
##' estimate(g,contrast=rbind(c(1,1,0),c(1,0,2)))
##' estimate(g,contrast=rbind(c(1,1,0),c(1,0,2)),null=c(1,2))
##'
##' ## Transformations
##' estimate(g,function(p) p[1]+p[2])
##' ## Multiple parameters
##' e <- estimate(g,function(p) c(p[1]+p[2],p[1]*p[2]))
##' e
##' vcov(e)
##'
##' ## Label new parameters
##' estimate(g,function(p) list("a1"=p[1]+p[2],"b1"=p[1]*p[2]))
##'
##' ## Marginalize
##' f <- function(p,data) 
##'   list(p0=lava:::expit(p[1] + p[3]*data[,"z"]),
##'        p1=lava:::expit(p[1] + p[2] + p[3]*data[,"z"]))
##' e <- estimate(g, f)
##' e
##' estimate(e,diff)
##' estimate(e,contrast=cbind(1,1))  
##'
##' ## Look at lava:::TN.curereg for example where gradient is supplied as attribute.
##' 
##' @S3method estimate default
estimate.default <- function(x,fun,data=model.frame(x),id1,id2,vcov=TRUE,level=0.95,iid=FALSE,contrast,null,...) {
  alpha <- 1-level
  alpha.str <- paste(c(alpha/2,1-alpha/2)*100,"",sep="%")
  if (!any(paste("score",class(x),sep=".") %in% methods("score"))) {
    if (!missing(fun))
      return(constrain(x,fun,data=data))
    if (!missing(contrast) | !(missing(null))) {
      p <- length(pars(x))
      if (missing(contrast)) contrast <- diag(p)
      if (missing(null)) null <- 0
      return(compare(x,contrast=contrast,null=null))
    }
    stop("transformation ('fun') or contrast matrix/null hypothesis ('contrast'/'null') needed")  
  }
  nn <- NULL
  U <- score(x,indiv=TRUE)  
  pp <- pars(x)
  n <- nrow(U)
  if (vcov) {
    iI <- vcov(x)*n
  } else {
    iI <- Inverse(information(e,type="hessian"))*n
  }
  iid0 <- U%*%iI
  if (missing(fun)) {
    if (iid) return(iid0)
    V <- crossprod(iid0/n)
  } else {
  form <- names(formals(fun))
    dots <- ("..."%in%names(form))
    form0 <- setdiff(form,"...")
    if (length(form0)==1 && !(form0%in%c("object","data"))) {
      names(formals(fun))[1] <- "p"
    }
    arglist <- c(list(object=x,data=data,p=pp),list(...))
    if (!dots) {
      arglist <- arglist[form0]
    }
    val <- do.call("fun",arglist)
    newfun <- NULL
    if (is.list(val)) {
      nn <- names(val)
      val <- do.call("cbind",val)
      newfun <- function(...) do.call("cbind",fun(...))
    }
    k <- NCOL(val)
    N <- NROW(val)
    D <- attributes(val)$grad
    if (is.null(D)) {
      D <- jacobian(function(p,...) {
        arglist$p <- p
        if (is.null(newfun))
          return(do.call("fun",arglist))
      return(do.call("newfun",arglist)) }, pp)      
    }
    if (NROW(val)<NROW(data)) { ## Apparently transformation not depending on data
      pp <- as.vector(val)
      iid2 <- iid0%*%t(D)
      V <- crossprod(iid2/n)
    } else {      
      if (k>1) { ## More than one parameter (and depends on data)
        D0 <- matrix(nrow=k,ncol=length(pp))
        for (i in seq_len(k)) {
          D0[i,] <- colMeans(D[seq(N)+(i-1)*N,,drop=FALSE])
        }
        D <- D0
        iid2 <- iid0%*%t(D)        
      } else {
        D <- colMeans(rbind(D))
        iid2 <- iid0%*%D
      }
      if (iid) return(list(iid1,iid2))
      pp <- as.vector(colMeans(cbind(val)))
      iid1 <- t(rbind(apply(cbind(val),1,function(x) x-pp)))
      if (N!=n) {
        if (missing(id2)) {          
          message("Assuming independence between model iid decomposition and new data frame")
          V <- crossprod(iid1/n) + crossprod(iid2/n)
        } else {
          iid1[id1,,drop=FALSE]+iid2[id2,,drop=FALSE]
        }
      } else {
        ##if (!missing(id))
        V <- crossprod((iid1+iid2)/n)
      }
    }
  }
  res <- cbind(pp,diag(V)^0.5)
  res <- cbind(res,res[,1]-qnorm(1-alpha/2)*res[,2],res[,1]+qnorm(1-alpha/2)*res[,2],(1-pnorm(abs(res[,1]/res[,2])))*2)
  colnames(res) <- c("Estimate","Std.Err",alpha.str,"P-value")
  if (!is.null(nn)) {
    rownames(res) <- nn
  } else {
   nn <- attributes(res)$varnames
   if (!is.null(nn)) rownames(res) <- nn
   if (is.null(rownames(res))) rownames(res) <- paste("p",seq(nrow(res)),sep="")
  }
  res <- structure(list(coef=res[,1],coefmat=res,vcov=V),class="estimate")
  if (!missing(contrast) |  !missing(null)) {
    p <- length(res$coef)
    if (missing(contrast)) contrast <- diag(p)
    if (missing(null)) null <- 0
    cc <- compare(res,contrast=contrast,null=null,vcov=V)
    res <- structure(c(res, list(compare=cc)),class="estimate")
  }
  return(res)  
}

estimate.glm <- function(x,...) {
  estimate.default(x,...)
}


##' @S3method print estimate
print.estimate <- function(x,...) {
  cat("\n")
  print(x$coefmat,...)
  if (!is.null(x$compare)) print(x$compare)    
}

##' @S3method vcov estimate
vcov.estimate <- function(x,...) {
  x$vcov
}

##' @S3method coef estimate
coef.estimate <- function(x,...) {
  x$coef
}
