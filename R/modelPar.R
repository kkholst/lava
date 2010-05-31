
`modelPar` <-
  function(x,p,...) UseMethod("modelPar")

###{{{ modelPar.lvmfit
modelPar.lvmfit <- function(x, p=pars(x), ...) modelPar(Model(x),p=p,...)

###}}} modelPar.lvmfit

###{{{ modelPar

modelPar.lvm <- function(x,p, debug=FALSE, ...) {
  npar <- index(x)$npar
  npar.mean <- index(x)$npar.mean
  Debug(list("npar=",npar),debug)
  Debug(list("npar.mean=",npar.mean),debug)
  Debug(list("p=",p),debug)

  if (length(p)!=npar & length(p)!=(npar+npar.mean)) stop("Wrong dimension of parameter vector!")  
  if (length(p)!=npar) { ## if meanstructure
      meanpar <- p[1:npar.mean]
      p. <- p[-c(1:npar.mean)]
    } else {
      meanpar <- NULL
      p. <- p
    }
  return(list(p=p.,meanpar=meanpar))
}

###}}} modelpar.lvm

###{{{ modelPar.multigroupfit
modelPar.multigroupfit <- function(x,p=pars(x),...) {
  modelPar(Model(x),p,...)
}
###}}}

###{{{ modelPar.multigroup
modelPar.multigroup <- function(x,p, debug=FALSE, ...) {
  npar <- x$npar
  npar.mean <- x$npar.mean
  k <- x$ngroup
  if (length(p)!=npar & length(p)!=(npar+npar.mean)) stop("Wrong dimension of parameter vector!")  
  if (length(p)!=npar) { ## if meanstructure
      meanpar <- p[1:npar.mean]
      p. <- p[-c(1:npar.mean)]
    } else {
      meanpar <- NULL
      p. <- p
    }


  parlist <- list(); for (i in 1:k) parlist[[i]] <- numeric(length(x$parlist[[i]]))
  if (!is.null(meanpar)) {
    meanlist <- list(); for (i in 1:k) meanlist[[i]] <- numeric(length(x$meanlist[[i]]))
  }

  if (length(p.)>0)
  for (i in 1:length(p.)) {
    for (j in 1:k) {
      idx <- match(paste("p",i,sep=""), x$parlist[[j]])
      if (!is.na(idx))
        parlist[[j]][idx] <- p.[i]
      if (!is.null(meanpar)) {
        midx <- match(paste("p",i,sep=""), x$meanlist[[j]])
        if (!is.na(midx))
          meanlist[[j]][midx] <- p.[i]
      }
    }
  }
  
  if (!is.null(meanpar)) {
    for (i in 1:length(meanpar)) {
      for (j in 1:k) {
        idx <- match(paste("m",i,sep=""), x$meanlist[[j]])
        if (!is.na(idx))
          meanlist[[j]][idx] <- meanpar[i]
      }
    }
  } else {
    meanlist <- NULL
  }
  p0 <- parlist
  for (i in 1:length(p0))
    p0[[i]] <- c(meanlist[[i]],parlist[[i]])  
  return(list(p=p0, par=parlist, mean=meanlist))
}

###}}}
