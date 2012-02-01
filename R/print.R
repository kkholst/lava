###{{{ print.lvm

`print.lvm` <-
function(x, ...) {
  k <- length(vars(x))
  cat("Latent Variable Model \n\twith: ", k, " variables.\n", sep="");
  if (k==0)
    return()
  cat("Npar=", index(x)$npar, "+", index(x)$npar.mean, "\n", sep="")
  cat("\n")
  ff <- formula(x,TRUE)
  for (f in ff) {
    oneline <- as.character(f); 
    cat(as.character(oneline),"\n")
  }
  invisible(x)
}

###}}} print.lvmz

###{{{ print.lvmfit
`print.lvmfit` <-
function(x,level=2,labels=FALSE,...) {
##    print(signif(coef(x), digits=digits), print.gap=2, quote=FALSE, ...)
##  if (index(x)$npar.reg==0 & level==0) level <- 1
##  coefs <- coef(x, level=level, ...);
##  coefs[coefs[,4]==0,4] <- 1e-16
##  printCoefmat(coefs, na.print="NA", signif.stars=signif.stars,...)
##  cat("---\n")
#  cat("\n", rep("-", 50), "\n\n", sep="");
  print(CoefMat(x,labels=labels,level=level,...),quote=FALSE,right=TRUE)
#  cat("\n", rep("-", 50), "\n\n", sep="");
  invisible(x)
##  invisible(coefs)
}
###}}} print.lvmfit

###{{{ print.lvmfit.randomslope

print.lvmfit.randomslope <- function(x,labels=FALSE,level=2,...) {
  print(CoefMat(x,labels=labels,level=level,...),quote=FALSE,right=TRUE)
  invisible(x)
}

###}}}

###{{{ print.multigroupfit

print.multigroupfit <- function(x,groups=NULL,...)  {
  if (is.null(groups)) {
    if (x$model$missing) {
      groups <- x$model$complete
      if (!is.null(e$model$mnames))
        x$model$names <- x$model$mnames
    } else {
      groups <- seq_len(length(x$model$lvm))
    }  
  }
  res <- coef(x,groups=groups,...)
  counter <- 0
  dots <- list(...)
  dots$groups <- groups
  level <- if (is.null(dots$level)) {
    dots$level <- 2
##    dots$level <- ifelse("lvmfit.randomslope"%in%class(x),2,9)
  }
  myargs <- c(list(x=x), dots)
  CC <- do.call("CoefMat.multigroupfit",myargs)
  for (cc in res) {
    counter <- counter+1
    cat(rep("-",50),"\n",sep="")
    cat("Group ", counter, sep="")
    myname <- x$model$names[counter]
    if (!is.null(myname) && !is.na(myname))
      cat(": ",myname,sep="")
    cat(" (n=",nrow(Model(x)$data[[groups[counter]]]), ")\n", sep="")
    print(CC[[counter]],quote=FALSE,right=TRUE)
  }
  cat("\n")
  invisible(x)
}

###}}} print.multigroupfit

###{{{ print.multigroup
print.multigroup <- function(x,...) {
  cat("\n")
  cat("Number of groups:", x$ngroup, "\n")
  cat("Number of free parameters (not counting mean-parameters):", x$npar,"\n")  
##  cat("Parameter-vector:", unlist(x$parlist), "\n\n")
  cat("Number of free mean parameters:", length(grep("m",x$mean)),"\n")  
##  cat("Mean-vector:", x$mean, "\n\n")
  invisible(x)
}
###}}} print.multigroup

###{{{ printmany

printmany <- function(A,B,nspace=1,name1=NULL,name2=NULL,digits=3,rownames=NULL,emptystr=" ",bothrows=TRUE,print=TRUE,...) {
  cA <- colnames(A); cB <- colnames(B)
  A <- format(A, digits=digits)
  B <- format(B, digits=digits)
  nA <- nrow(A); nB <- nrow(B)
  if (nrow(A)<nrow(B)) {
    rA <- rownames(A)
    A <- rbind(A, matrix("", nrow=nB-nA, ncol=ncol(A)))
  }
  if (nrow(B)<nrow(A)) {
    rB <- rownames(B)
    B <- rbind(B, rep("", nrow=nA-nB, ncol=ncol(B)))
  }
  if (!is.null(rownames) & length(rownames)==nrow(A))
    rownames(A) <- rownames(B) <- rownames
  res <- cbind(A, matrix("", nrow=nrow(A), ncol=nspace)); colnames(res) <- c(colnames(A), rep(emptystr,nspace))
  if (!is.null(name1)) {
    oldname <- colnames(res)
    res <- cbind(rep("",nrow(res)), rownames(res), res); colnames(res) <- c(name1,"",oldname)
    rownames(res) <- rep("",nrow(res))
  }
  if (!is.null(name2)) {
    oldname <- colnames(res)
    res <- cbind(res,rep("",nrow(res))); colnames(res) <- c(oldname,name2)
  }
  if (!identical(rownames(A),rownames(B)) | bothrows)
    res <- cbind(res, rownames(B))
  res <- cbind(res, B)
  if (print) print(res, quote=FALSE,...)
  invisible(res)
}

###}}} printmany
