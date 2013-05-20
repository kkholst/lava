##' @export
`score` <-
function(x,...) UseMethod("score")

###{{{ score.lvm

##' @S3method score lvm
score.lvm <- function(x, data, p, model="gaussian", S, n, mu=NULL, weight=NULL, weight2=NULL, debug=FALSE, reindex=FALSE, mean=TRUE, constrain=TRUE, indiv=TRUE,...) {

  cl <- match.call()
  lname <- paste(model,"_score.lvm",sep="")
  if (!exists(lname)) {
    lname <- paste(model,"_gradient.lvm",sep="")
    mygrad <- get(lname)    
    scoreFun <- function(...) -mygrad(...)
    if (is.null(mygrad)) {
      stop("Missing gradient")
    }      
  } else {
    scoreFun <- get(lname)
  }
  
  if (missing(data) || is.null(data)) {
    cl[[1]] <- scoreFun
    score <- eval.parent(cl)
    return(rbind(score))
  }    
  
  if (is.null(index(x)$dA) | reindex)
    x <- updatelvm(x,zeroones=TRUE,deriv=TRUE)
  
  xfix <- colnames(data)[(colnames(data)%in%parlabels(x,exo=TRUE))]
  xconstrain <- intersect(unlist(lapply(constrain(x),function(z) attributes(z)$args)),index(x)$manifest)

  Debug(xfix,debug)
  if (missing(n)) {
    n <- nrow(data)
  }

  if (length(xfix)>0 | length(xconstrain)>0) { ##### Random slopes!
    x0 <- x
    if (length(xfix)>0) {
      Debug("random slopes...",debug)
      nrow <- length(vars(x))
      xpos <- lapply(xfix,function(y) which(regfix(x)$labels==y))
      colpos <- lapply(xpos, function(y) ceiling(y/nrow))
      rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
      myfix <- list(var=xfix, col=colpos, row=rowpos)
      for (i in 1:length(myfix$var))
        for (j in 1:length(myfix$col[[i]])) {
          regfix(x0, from=vars(x0)[myfix$row[[i]][j]],to=vars(x0)[myfix$col[[i]][j]]) <-
            data[1,myfix$var[[i]]]
        }
      index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)    
    }
    pp <- modelPar(x0,p)
    p0 <- with(pp, c(meanpar,p))
    k <- length(index(x0)$manifest)

    myfun <- function(ii) {
      if (length(xfix)>0)
        for (i in 1:length(myfix$var)) {
          index(x0)$A[cbind(myfix$row[[i]],myfix$col[[i]])] <- data[ii,myfix$var[[i]]]
        }
      return(scoreFun(x0,data=data[ii,], p=with(pp,c(meanpar,p)),weight=weight[ii,,drop=FALSE],weight2=weight2[ii,,drop=FALSE],model=model,debug=debug,indiv=indiv,...))
    }
    score <- t(sapply(1:nrow(data),myfun))
    if (!indiv) {
      score <- colSums(rbind(score))
    }
    if (length(score)<length(p))  score <- c(score,rep(0,length(p)-length(score)))
    return(score)
  }
  cl$constrain <- FALSE
  cl[[1]] <- scoreFun
  score <- eval.parent(cl)
  if (is.null(dim(score))) score <- rbind(score)
  if (NCOL(score)<length(p))  score <- cbind(rbind(score),rep(0,length(p)-NCOL(score)))

#  score <- eval(cl,parent.frame())
  return(score)
}

###}}} score.lvm
  
###{{{ score.lvm.missing

##' @S3method score lvm.missing
score.lvm.missing <- function(x,
                          p=pars(x), estimator=x$estimator,
                              weight=Weight(x$estimate),
                              indiv=FALSE,
                              list=FALSE,
                              ...) {
  S <- score(x$estimate$model0, p=p, model=estimator, weight=weight,
               indiv=indiv,...)
  if (indiv & !list) {
    S0 <- matrix(ncol=length(p),nrow=length(x$order))
    rownames(S0) <- 1:nrow(S0)
    myorder <- x$orderlist
    if (length(x$allmis)>0)
      myorder[[x$allmis]] <- NULL
    for (i in 1:length(S))
      S0[myorder[[i]],] <- S[[i]]
    if (length(x$allmis)>0) {
      S0 <- S0[-x$orderlist[[x$allmis]],]
    }
    S0[is.na(S0)] <- 0
    return(S0)
  }
  return(S)
}

###}}} score.lvm.missing

###{{{ score.multigroupfit

##' @S3method score multigroupfit
score.multigroupfit <- function(x,p=pars(x), weight=Weight(x), estimator=x$estimator, ...) {
  score(x$model0, p=p, weight=weight, model=estimator,...)
}

###}}} score.multigroupfit

###{{{ score.multigroup

##' @S3method score multigroup
score.multigroup <- function(x,data=x$data,weight=NULL,p,indiv=FALSE,...) {
  rm <- procrandomslope(x)
  pp <- with(rm, modelPar(model,p)$p)
  parord <- modelPar(rm$model,1:with(rm$model,npar+npar.mean))$p
  S <- list()
  for (i in 1:x$ngroup) {
    S0 <- rbind(score(x$lvm[[i]],p=pp[[i]],data=data[[i]],weight=weight[[i]],indiv=indiv,...))
    S1 <- matrix(ncol=length(p),nrow=nrow(S0))
    S1[,parord[[i]]] <- S0
    S <- c(S, list(S1))
  }
  if (indiv) return(S)
  res <- matrix(0,nrow=1,ncol=length(p))
  for (i in 1:x$ngroup)
    res[,parord[[i]]] <- res[,parord[[i]]]  + S[[i]][,parord[[i]]]
  return(as.vector(res))
}

###}}} score.multigroup

###{{{ score.lvmfit

##' @S3method score lvmfit
score.lvmfit <- function(x, data=model.frame(x), p=pars(x), model=x$estimator, weight=Weight(x), weight2=x$data$weight2, ...) {
  score(x$model0,data=data,p=p,model=model,weight=weight,weight2=weight2,...)
}

###}}} score.lvmfit
