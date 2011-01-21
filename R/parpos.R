## parpos <- function(x,mean=TRUE,p,...) {
##   if (class(x)[1]%in%c("multigroup","multigroupfit")) {
##     if (class(x)[1]%in%"multigroupfit")
##       return(parpos.multigroup(Model(x),mean=mean,p=p,...))
##     return(parpos.multigroup(x,mean=mean,p=p,...))
##   }
##   if (!missing(p)) {
##     if (!is.character(p)) p <- names(p)
##     idx1 <- match(p,coef(Model(x),mean=mean,fix=FALSE))
##     idx11 <- match(p,coef(Model(x),mean=mean,fix=FALSE,labels=TRUE))
##     res <- (union(idx1,idx11)); names(res) <- p
##     res <- idx1; res[!is.na(idx11)] <- idx11[!is.na(idx11)]
##     names(res) <- p
##     ord <- order(res)
##     res <- sort(res)
##     attributes(res)$ord <- ord    
##     return(res)
##   }
##   if (mean)
##     nn <- with(index(x),matrices(x,1:npar+npar.mean,meanpar=1:npar.mean)) ## Position of parameters
##   else nn <- with(index(x),matrices(x,1:npar,NULL))
##   nn$A[index(x)$M0!=1] <- 0
##   nn$P[index(x)$P0!=1] <- 0
##   nn$v[index(x)$v0!=1] <- 0
##   nn
## }

parpos <- function(x,mean=TRUE,p,...) {
  if (class(x)[1]%in%c("multigroup","multigroupfit")) {
    if (class(x)[1]%in%"multigroupfit")
      return(parpos.multigroup(x$model0,mean=mean,p=p,...))
    return(parpos.multigroup(x,mean=mean,p=p,...))
  }
  if (!missing(p)) {    
    if (!is.character(p)) p <- names(p)
    cc1 <- coef(Model(x),mean=mean,fix=FALSE)
    cc2 <- coef(Model(x),mean=mean,fix=FALSE,labels=TRUE)
    idx1 <- na.omit(match(p,cc1))
    idx11 <- na.omit(match(p,cc2))
    res <- (union(idx1,idx11));
    if (length(res)!=length(p)) {
      names(res) <- cc1[res]
    } else {      
      names(res) <- p
    }
    ##    res <- idx1; res[!is.na(idx11)] <- idx11[!is.na(idx11)]
    ##    names(res) <- p
    ord <- order(res)
    res <- sort(res)
    attributes(res)$ord <- ord    
    return(res)
  }
  if (mean)
    nn <- with(index(x),matrices(x,1:npar+npar.mean,meanpar=1:npar.mean)) ## Position of parameters
  else nn <- with(index(x),matrices(x,1:npar,NULL))
  nn$A[index(x)$M0!=1] <- 0
  nn$P[index(x)$P0!=1] <- 0
  nn$v[index(x)$v0!=1] <- 0
  nn
}

parpos.multigroup <- function(x,mean=TRUE,p,...) {
  if (missing(p)) {
    p <- unique(unlist(lapply(x$lvm, function(z) setdiff(parlabels(z),names(constrain(z))) )))
  }
  if (!is.character(p)) p <- names(p)
  p0 <- rep(NA,with(x,npar+npar.mean));
  names(p0) <- c(x$mean,x$par)
  for (i in 1:length(x$lvm)) {
    cur <- parpos(x$lvm[[i]],p=p)
    if (length(cur)>0) {
      p0[c(x$meanpos[[i]],x$parpos[[i]])[cur]] <- names(cur)
      p <- p[-match(names(cur),p)]
    }    
    if (length(p)==0) break;
  }
  return(p0)    
}

   
