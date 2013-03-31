###{{{ summary.lvm

##' @S3method summary lvm
`summary.lvm` <-
function(object,...) {
  k <- length(vars(object))
  ## cat("Latent Variable Model \n\twith: ", k, " variables.\n", sep="");
  print(object)
  if (k==0)
    return()
##  cat("Npar=", index(object)$npar, "+", index(object)$npar.mean, "\n", sep="")
  cat("\n")
  print(regression(object))
  print(covariance(object))
  print(intercept(object))
  if (length(object$exfix)>0) {
    cat("Additional parameters:\n")    
    val <- unlist(object$exfix)
    M <- rbind(val); colnames(M) <- names(val)
    rownames(M) <- "   "
    print(M,quote=FALSE)
  }
  if (length(constrain(object))>0) {
    cat("Non-linear constraints:\n")
    print(constrain(object),quote=FALSE)
  }
    
  ## printmany(object$cov, printmany(object$covpar, object$covfix, name1="Labels:", name2="Fixed:", print=FALSE), name1="covariance:")

  cat("\n")
}

###}}} summary.lvm

###{{{ summary.lvmfit

##' @S3method summary lvmfit
`summary.lvmfit` <-
function(object,std="xy", level=9, labels=2, ...) {
  cc <- CoefMat(object,labels=labels,std=std,level=level,...)
  mycoef <- coef(object,level=9)
  nlincon <- attributes(mycoef)$nlincon
  nonexo <- setdiff(vars(object),index(Model(object))$exogenous)
  attributes(mycoef) <- attributes(mycoef)[1:2]
  mygof <- object$opt$summary.message
  if (is.null(mygof)) {
    mygof <- gof
  }
  if (class(object)[1]=="lvm.missing") {
    nn <- unlist(lapply(object$multigroup$data, nrow))
    nc <- nn[object$cc]
    if (length(nc)==0) nc <- 0
    ngroup <- object$multigroup$ngroup
    res <- list(object=object, coef=mycoef, coefmat=cc, nlincon=nlincon, gof=mygof(object), n=sum(nn), nc=nc, ngroup=ngroup, varmat=modelVar(object)$P[nonexo,nonexo], latent=latent(object), opt=object$opt, vcov=vcov(object), estimator=object$estimator)
  } else {
    n <- nrow(model.frame(object))
    if (is.null(n)) n <- model.frame(object)$n
    res <- list(coef=mycoef, coefmat=cc, nlincon=nlincon, gof=mygof(object), n=n, nc=n, latent=latent(object), opt=object$opt, vcov=vcov(object), estimator=object$estimator)##, varmat=modelVar(object)$P[nonexo,nonexo])
  }
  class(res) <- "summary.lvmfit"
  res
}

##' @S3method print summary.lvmfit
print.summary.lvmfit <- function(x,varmat=TRUE,...) {
  if (!is.null(x$control$method)) {
    l2D <- sum(x$opt$grad^2)
    rnkV <- qr(x$vcov)$rank
    if (l2D>1e-2) warning("Possible problems with convergence!")    
    cat("||score||^2=",l2D,"\n",sep="")
    np <- nrow(x$vcov)
    if (rnkV<np) warning("Possible problems with identification (rank(informaion)=",rnkV,"<",np,"!")
  }
  cat("Latent variables:", x$latent, "\n")
  cat("Number of rows in data=",x$n,sep="")
  if (x$nc!=x$n) {
    cat(" (",x$nc," complete cases, ", x$ngroup, " groups)",sep="")    
  }; cat("\n")
  cat(rep("-", 50), "\n", sep="");
  print(x$coefmat,quote=FALSE,right=TRUE)
##  if (varmat) {
##    cat("\nResidual covariance matrix:\n")
##    print(x$varmat)
##  }
  if (!is.null(x$nlincon)) {
    cat("\nNon-linear constraints:\n")
    printCoefmat(x$nlincon,signif.stars=FALSE)
  }
  cat(rep("-", 50), "\n", sep="");
  cat("Estimator:",x$estimator,"\n")
  cat(rep("-", 50), "\n", sep="");
  if (!is.null(x$gof)) {
    if (class(x$gof)[1]=="list") {
      for (i in x$gof) {
        print(i)
      }
    } else {
      print(x$gof,optim=FALSE)
    }
    cat(rep("-", 50), "\n", sep="");
  }
  invisible(x)
}

##' @S3method coef summary.lvmfit
coef.summary.lvmfit <- function(object,...) object$coef

###}}} summary.lvmfit

###{{{ summary.multigroupfit

##' @S3method summary multigroupfit
summary.multigroupfit <- function(object,groups=NULL,...) {
  if (is.null(groups) | length(groups)==0) {
    if (object$model$missing) {
      groups <- object$model$complete
      if (length(groups)==0)
        groups <- seq_len(object$model0$ngroup)
    } else {
      groups <- seq_len(object$model$ngroup)
    }    
  }
  cc <- CoefMat.multigroupfit(object,groups=groups,...)
  res <- list(coef=coef(object,level=2,groups=groups,...), object=object, coefmat=cc, gof=gof(object), object=object, opt=object$opt, latent=object$latent, estimator=object$estimator)
  class(res) <- "summary.multigroupfit"
  res
}

##' @S3method print summary.multigroupfit
print.summary.multigroupfit <- function(x,...) {
  l2D <- sum(x$opt$grad^2)
  if (l2D>1e-2) warning("Possible problems with convergence!")
  cat("||score||^2=",l2D,"\n")
  cat("Latent variables:", x$latent, "\n")
  print(x$object,...)  
##  cat(rep("-", 50), "\n\n", sep="");
  ##print(x$coefmat,quote=FALSE,right=TRUE)
  cat(rep("-", 50), "\n", sep="");
  if (!is.null(attributes(x$coefmat)$nlincon)) {
    cat("Non-linear constraints:\n")
    print(attributes(x$coefmat)$nlincon)
    cat(rep("-", 50), "\n", sep="");
  }
  cat("Estimator:",x$estimator,"\n")  
  cat(rep("-", 50), "\n", sep="");
  if (!is.null(x$gof)) {
    print(x$gof)
    cat(rep("-", 50), "\n", sep="");
  }
  invisible(x)  
}

###}}} summary.multigroupfit

###{{{ summary.multigroup

##' @S3method summary multigroup
summary.multigroup <- function(object,...) {
  for (m in object$lvm)
    print(m,...)
  print(object)
  invisible(object)
}

###}}}
