##' @S3method estimate default
estimate.default <- function(x,fun,data=model.frame(x),vcov=TRUE,level=0.95,iid=FALSE,...) {
  alpha <- 1-level
  alpha.str <- paste(c(alpha/2,1-alpha/2)*100,"",sep="%")
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
    k <- NCOL(val)
    N <- NROW(val)
    D <- attributes(val)$grad
    if (is.null(D)) {
      D <- jacobian(function(p,...) {
        arglist$p <- p
        do.call("fun",arglist)  }, pp)
    }
    if (NROW(val)<NROW(data)) { ## Apparently transformation not depending on data
      pp <- as.vector(val)
      iid <- iid0%*%t(D)
    } else {      
      if (k>1) { ## More than one parameter (and depends on data)
        D0 <- matrix(nrow=k,ncol=length(pp))
        for (i in seq_len(k)) {
          D0[i,] <- colMeans(D[seq(N)+(i-1)*N,,drop=FALSE])
        }
        D <- D0
        iid <- iid0%*%t(D)
        
      } else {
        D <- colMeans(rbind(D))
        iid <- iid0%*%D
      }
      pp <- as.vector(colMeans(cbind(val)))
      iid2 <- t(rbind(apply(cbind(val),1,function(x)-pp)))
      iid <- iid2+iid
    }    
    V <- crossprod(iid/n)
  }
  res <- cbind(pp,diag(V)^0.5)
  res <- cbind(res,res[,1]-qnorm(1-alpha/2)*res[,2],res[,1]+qnorm(1-alpha/2)*res[,2],(1-pnorm(abs(res[,1]/res[,2])))*2)
  colnames(res) <- c("Estimate","Std.Err",alpha.str,"P-value")
  nn <- attributes(res)$varnames
  if (is.null(nn)) nn <- rep("",nrow(res))
  rownames(res) <- nn
  res <- structure(list(coef=res,vcov=V),class="estimate")
  return(res)  
}

##' @S3method print estimate
print.estimate <- function(x,...) {
  cat("\n")
  print(x$coef,...)
}

##' @S3method vcov estimate
vcov.estimate <- function(x,...) {
  x$vcov
}

##' @S3method coef estimate
coef.estimate <- function(x,...) {
  x$coef
}
