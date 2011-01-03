###{{{ summary.lvm

`summary.lvm` <-
function(object,...) {
  k <- length(vars(object))
  cat("Latent Variable Model \n\twith: ", k, " variables.\n", sep="");
  if (k==0)
    return()
  cat("Npar=", index(object)$npar, "+", index(object)$npar.mean, "\n", sep="")
  cat("\n")
  print(regression(object))
  print(covariance(object))
  print(intercept(object))
  ## printmany(object$cov, printmany(object$covpar, object$covfix, name1="Labels:", name2="Fixed:", print=FALSE), name1="covariance:")
  ## cat("\n")
  ## A <- as(Graph(object), Class="matrix")
  ## printmany(A, printmany(object$par, object$fix, name1="Labels:", name2="Fixed", print=FALSE),
  ##           name1="Adjancency:")
  ## cat("\n")
  ## mu <- matrix(unlist(object$mean),nrow=1); rownames(mu) <- "Mean"; colnames(mu) <- vars(object)
  ## print(mu, print.gap=2, quote=FALSE);
  cat("\n")
}

###}}} summary.lvm

###{{{ summary.lvmfit

`summary.lvmfit` <-
function(object,std="xy", level=9, labels=2, ...) {
  cc <- CoefMat(object,labels=labels,std=std,level=level,...)
  mycoef <- coef(object,level=9)
  nlincon <- attributes(mycoef)$nlincon
  nonexo <- setdiff(vars(object),index(Model(object))$exogenous)
  attributes(mycoef) <- attributes(mycoef)[1:2]
  if (class(object)[1]=="lvm.missing") {
    nn <- unlist(lapply(object$multigroup$data, nrow))
    nc <- nn[object$cc]
    if (length(nc)==0) nc <- 0
    ngroup <- object$multigroup$ngroup
    res <- list(object=object, coef=mycoef, coefmat=cc, nlincon=nlincon, gof=gof(object), n=sum(nn), nc=nc, ngroup=ngroup, varmat=modelVar(object)$P[nonexo,nonexo])
  } else {
    n <- nrow(model.frame(object))
    res <- list(object=object, coef=mycoef, coefmat=cc, nlincon=nlincon, gof=gof(object), n=n, nc=n)##, varmat=modelVar(object)$P[nonexo,nonexo])
  }
  class(res) <- "summary.lvmfit"
  res
}

print.summary.lvmfit <- function(x,varmat=TRUE,...) {
  l2D <- sum(x$object$opt$grad^2)
  rnkV <- qr(vcov(x$object))$rank
  if (l2D>1e-2) warning("Possible problems with convergence!")
  cat("||score||^2=",l2D,"\n",sep="")
  np <- nrow(vcov(x$object))
  if (rnkV<np) warning("Possible problems with identification (rank(informaion)=",rnkV,"<",np,"!")
  cat("Latent variables:", latent(x$object), "\n")
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
  cat("Estimator:",x$object$estimator,"\n")
  cat(rep("-", 50), "\n", sep="");
  if (!is.null(x$gof)) {
    print(x$gof,optim=FALSE)
    cat(rep("-", 50), "\n", sep="");
  }
  invisible(x)
}

coef.summary.lvmfit <- function(object,...) object$coef

###}}} summary.lvmfit

###{{{ summary.multigroupfit

summary.multigroupfit <- function(object,...) {
  cc <- CoefMat.multigroupfit(object,...) 
  res <- list(coef=coef(object), object=object, coefmat=cc, gof=gof(object), object=object)
  class(res) <- "summary.multigroupfit"
  res
}

print.summary.multigroupfit <- function(x,...) {
  l2D <- sum(x$object$opt$grad^2)
  if (l2D>1e-2) warning("Possible problems with convergence!")
  cat("||score||^2=",l2D,"\n")
  cat("Latent variables:", latent(x$object), "\n")
  print(x$object,...)
##  cat(rep("-", 50), "\n\n", sep="");
  ##  print(x$coefmat,quote=FALSE,right=TRUE)
  cat(rep("-", 50), "\n", sep="");
  if (!is.null(attributes(x$coefmat)$nlincon)) {
    cat("Non-linear constraints:\n")
    print(attributes(x$coefmat)$nlincon)
    cat(rep("-", 50), "\n", sep="");
  }
  cat("Estimator:",x$object$estimator,"\n")  
  cat(rep("-", 50), "\n", sep="");
  if (!is.null(x$gof)) {
    print(x$gof)
    cat(rep("-", 50), "\n", sep="");
  }
  invisible(x)  
}

###}}} summary.multigroupfit

###{{{ summary.multigroup

summary.multigroup <- function(object,...) {
  for (m in object$lvm)
    print(m,...)
  print(object)
  invisible(object)
}

###}}}
