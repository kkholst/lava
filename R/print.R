###{{{ print.lvm

##' @S3method print lvm
`print.lvm` <-
function(x, ...) {
  res <- NULL
  myhooks <- gethook("print.hooks")
  for (f in myhooks) {
    res <- do.call(f, list(x=x,...))
  }
  if (is.null(res)) {
    k <- length(vars(x))
    L <- rep(FALSE,k); names(L) <- vars(x); L[latent(x)] <- TRUE
    cat("Latent Variable Model\n") ##;" \n\twith: ", k, " variables.\n", sep="");
    if (k==0)
      return()
    ff <- formula(x,TRUE)
    R <- c()
    for (f in ff) {
      oneline <- as.character(f);
      y <- gsub(" ","",strsplit(f,"~")[[1]][1])
##      if (!(y %in% endogenous(m)))
      {
        col1 <- as.character(oneline)
        D <- attributes(distribution(x)[[y]])$family
        col2 <- x$attributes$type[y]
        if (is.null(col2) || is.na(col2)) col2 <- "Normal"
        if (!is.null(D$family)) col2 <- paste(D$family,sep="")
        if (!is.null(D$link)) col2 <- paste(col2,"(",D$link,")",sep="")
        if (!is.null(D$par)) col2 <- paste(col2,"(",paste(D$par,collapse=","),")",sep="")
        if (L[y]) col2 <- paste(col2,", Latent",sep="")  
        R <- rbind(R,c(col1,"  ",col2))
      }
    }
    if (length(R)>0) {
      rownames(R) <- rep("",nrow(R)); colnames(R) <- rep("",ncol(R))
      print(R,quote=FALSE,...)
    }
    cat("\n")
    cat("Number of free parameters: ", with(index(x),npar+npar.mean+npar.ex),"\n", sep="")

##      oneline <- as.character(f); 
##      cat(as.character(oneline),"\n")

  }
  invisible(x)
}

###}}} print.lvm

###{{{ print.lvmfit

##' @S3method print lvmfit
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

##' @S3method print lvmfit.randomslope
print.lvmfit.randomslope <- function(x,labels=FALSE,level=2,...) {
  print(CoefMat(x,labels=labels,level=level,...),quote=FALSE,right=TRUE)
  invisible(x)
}

###}}}

###{{{ print.multigroupfit

##' @S3method print multigroupfit
print.multigroupfit <- function(x,groups=NULL,...)  {
  if (is.null(groups)) {
    if (x$model$missing) {
      modelclass <- attributes(x$model0)$modelclass
      nmis <- attributes(x$model0)$nmis
      orggroup <- unique(modelclass)
      groupn <- unlist(lapply(orggroup,function(i) sum(modelclass==i)))
      cumsumgroup <- cumsum(c(0,groupn))
      groups <- unlist(lapply(orggroup,function(i)
                              which.min(nmis[which(modelclass==i)])+cumsumgroup[i])) ##  groups with max. number of variables
                       ##      suppressMessages(browser())
                       ## if (length(groups)==0)
                       ##   groupedDataps <- seq_len(x$model0$ngroup)      
      for (i in seq_len(length(groups))) {
        if (nmis[groups[i]]>0) warning("No complete cases in group ",i,". Showing results of group with max number of variables. All coefficients can be extracted with 'coef'. All missing pattern groups belonging to this sub-model can be extracted by calling: coef(..., groups=c(",paste(which(modelclass==i),collapse=","),"))")
      }
      if (!is.null(x$model$mnameses))
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
  ##  suppressMessages(browser())
  myargs <- c(list(x=x), dots)
  myargs$groups <- groups
  CC <- do.call("CoefMat.multigroupfit",myargs)  
  for (cc in res) {
    counter <- counter+1
    cat(rep("-",50),"\n",sep="")
    cat("Group ", counter, sep="")
    myname <- x$model$names[counter]
    if (!is.null(myname) && !is.na(myname))
      cat(": ",myname,sep="")
    if (!x$model$missing) cat(" (n=",nrow(Model(x)$data[[groups[counter]]]), ")", sep="")
    cat("\n")
    print(CC[[counter]],quote=FALSE,right=TRUE)
  }
  cat("\n")
  invisible(x)
}

###}}} print.multigroupfit

###{{{ print.multigroup

##' @S3method print multigroup
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
