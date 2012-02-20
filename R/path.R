
path <- function(object,to,from,...) UseMethod("path")
path.lvmfit <- function(object,to=NULL,from,...) {
  mypath <- path(Model(object),to,from,...)
  cc <- coef(object,level=9) ## All parameters (fixed and variable)

  #cc0 <- coef(object,level=1) ## Estimated parameters
  cc0 <- coef(object,level=2) ## Estimated parameters
  i1 <- na.omit(match(rownames(cc),rownames(cc0)))
  idx.cc0 <-  which(rownames(cc)%in%rownames(cc0)); ## Position of estimated parameters among all parameters
  S <- matrix(0,nrow(cc),nrow(cc)); rownames(S) <- colnames(S) <- rownames(cc)
  V <- object$vcov
  npar.mean <- index(object)$npar.mean
#  if (object$control$meanstructure & npar.mean>0)
#    V <- V[-c(1:npar.mean),-c(1:npar.mean)]
  S[idx.cc0,idx.cc0] <- V[i1,i1]  ## "Covariance matrix" of all parameters
  

  
  idx <- list()
  coefs <- list()
  V <- list()
  for (i in 1:length(mypath)) {
    xx <- mypath[[i]]
    ii <- c()
    for (j in 1:(length(xx)-1)) {
      st <- paste(xx[j+1], "<-", xx[j],sep="")
      ii <- c(ii, match(st,rownames(cc)))
    }
    idx <- c(idx, list(ii)) 
    V <- c(V, list(S[ii,ii]))
    coefs <- c(coefs, list(cc[ii]))
  }

  edges <- list()
  for (i in 1:length(mypath)) {
    p0 <- mypath[[i]]
    ee <- c()
    for (i in 1:(length(p0)-1)) {
      ee <- c(ee, paste(p0[i],p0[i+1],sep="~"))
    }
    edges <- c(edges, list(ee))
  }
  res <- list(idx=idx,V=V,coef=coefs, path=mypath, edges=edges)
  return(res)
}

path.lvm <- function(object,to=NULL,from,...) path(Graph(object),to=to,from=from,...)
path.graphNEL <- function(object,to,from,...) {
  if (class(to)[1]=="formula") {
    fvar <- extractvar(to)
    if (length(fvar$x)==1 & length(fvar$y)==1)
      return(path(object,to=fvar$y,from=fvar$x))
    res <- list()
    for (y in fvar$y) {
      for (x in fvar$x) {
        cat("x=",x, " y=",y, "\n")
        res <- c(res, list(path(object,to=y,from=x)))
      }
    }
    return(res)
  }
  ff <- function(g,from=1,to=NULL,res=list()) {    
    M <- edgeMatrix(g)
    i1 <- which(M[1,]==from)
    for (i in i1) {
      e <- M[,i]; newto <- e[2];
      if (is.null(to) || M[2,i]==to) {
        res <- c(res, list(M[,i]))
      }
      newpath <- ff(g,from=newto,to=to,list())
      if (length(newpath)>0)
      for (j in 1:length(newpath)) {
        if (is.null(to) || (tail(newpath[[j]],1)==to))
          res <- c(res, list(c(M[,i],newpath[[j]][-1])))
      }
    }
    return(res)
  }
  idxfrom <- ifelse(is.numeric(from),from,which(from==nodes(object)))
  reachable <- acc(object,nodes(object)[idxfrom])[[1]]

  if (is.null(to)) {
    idxto <- reachable
  } else {
    idxto <- ifelse(is.numeric(to),to,which(to==nodes(object)))
  }

  if (!(nodes(object)[idxto] %in% names(reachable)))
##    return(structure(NULL,to=to[1],from=from[1]))
    return(NULL)
  ##    stop("No directional relationship between variables")
  
  mypaths <- ff(object,idxfrom,idxto)
  res <- list()
  for (i in 1:length(mypaths)) {
    res <- c(res, list(nodes(object)[mypaths[[i]]]))
  }
  return(res)
}
