
'%s%' <-
function(st1,st2) {
  paste(st1, st2, sep="")
}

###{{{ parlabels

parlabels <- function(x) {
  c(unlist(intfix(x)[unlist(lapply(intfix(x), function(y) !is.na(y) & !is.numeric(y)))]),
    regfix(x)$labels[!is.na(regfix(x)$labels)],
    covfix(x)$labels[!is.na(covfix(x)$labels)])
}

###}}} parlabels

###{{{ describecoef

describecoef <- function(x,par,from,to,mean=TRUE) {  
  p <- coef(x, mean=mean)  
  if (!missing(from)) {
    st1 <- paste(to,"<-",from,sep="")
    st2 <- paste(to,"<->",from,sep="")
    st3 <- paste(from,"<->",to,sep="")
    pos <- na.omit(match(unique(c(st1,st2,st3)),p))
    attributes(pos) <- NULL
    return(pos)
  }
  res <- strsplit(p,"<->")
  var.idx <- which(unlist(lapply(res,length))>1) ## Variance parameters
  rest.idx <- setdiff(1:length(p),var.idx)
  res[rest.idx] <- strsplit(p[rest.idx],"<-")
  mean.idx <- which(unlist(lapply(res,length))==1) ## Mean parameters
  reg.idx <- setdiff(rest.idx,mean.idx)
  names(res)[mean.idx] <- paste("m",1:length(mean.idx),sep="")
  for (i in var.idx)
    attr(res[[i]],"type") <- "cov"
  for (i in mean.idx)
    attr(res[[i]],"type") <- "mean"
  for (i in reg.idx)
    attr(res[[i]],"type") <- "reg"
  if (missing(par))
    return(res)
  return(res[par])
}

###}}}

###{{{ substArg

substArg <- function(x,env,...) {
  if (!missing(env)) {
    a <- with(env,substitute(x))
#    a <- substitute(x,environment(env))
  } else {
    a <- substitute(x)
  }
  myclass <- tryCatch(class(eval(a)),error=function(e) NULL)
  if (is.null(myclass) || myclass=="name") {
#  if (is.null(myclass)) {
    res <- unlist(sapply(as.character(a),
                         function(z) {
                           trimmed <- gsub(" ","",z,fixed=TRUE)
                           val <- strsplit(trimmed,"+",fixed=TRUE)
                           if (val[1]=="") val <- NULL
                           val
                         })); attributes(res)$names <- NULL
    return(res)
  }
  return(eval(a))
}

## g <- function(zz,...) {
##   env=new.env(); assign("x",substitute(zz),env)
##   substArg(zz,env=env)
## }
## h <- function(x,...) {
##   env=new.env(); assign("x",substitute(x),env)  
##   substArg(x,env=TRUE)
## }

###}}}

###{{{ procrandomslope

procrandomslope <- function(object,data=object$data,...) {
  Xfix <- FALSE
  xfix <- myfix <- list()
  xx <- object
  for (i in 1:object$ngroup) {
    x0 <- object$lvm[[i]]
    data0 <- data[[i]]
    xfix0 <- colnames(data0)[(colnames(data0)%in%parlabels(x0))]
    xfix <- c(xfix, list(xfix0))
    if (length(xfix0)>0) { ## Yes, random slopes
      Xfix<-TRUE
    }
    xx$lvm[[i]] <- x0    
  }   
  if (Xfix) {
    for (k in 1:object$ngroup) {
      x1 <- x0 <- object$lvm[[k]]
      data0 <- data[[k]]        
      nrow <- length(vars(x0))
      xpos <- lapply(xfix[[k]],function(y) which(regfix(x0)$labels==y))
      colpos <- lapply(xpos, function(y) ceiling(y/nrow))
      rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
      myfix0 <- list(var=xfix[[k]], col=colpos, row=rowpos)
      myfix <- c(myfix, list(myfix0))        
      for (i in 1:length(myfix0$var))
        for (j in 1:length(myfix0$col[[i]])) 
          regfix(x0,
                 from=vars(x0)[myfix0$row[[i]][j]],to=vars(x0)[myfix0$col[[i]][j]]) <-
                   colMeans(data0[,myfix0$var[[i]],drop=FALSE],na.rm=TRUE)
      index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
      object$lvm[[k]] <- x0
      yvars <- endogenous(x0)
      #parkeep <- c(parkeep, parord[[k]][coef(x1,mean=TRUE)%in%coef(x0,mean=TRUE)])
    }
#    parkeep <- sort(unique(parkeep))
    object <- multigroup(object$lvm,data,fix=FALSE,exo.fix=FALSE)
  }
  return(list(model=object,fix=myfix))
}

###}}}

###{{{ fixsome function

fixsome <- function(x, exo.fix=TRUE, measurement.fix=TRUE, S, mu, n, data, x0=FALSE, debug=TRUE, na.method="complete.obs") {

  var.missing <- c()
  if (!missing(data) | !missing(S)) {
        
    if (!missing(data)) {
      dd <- procdata.lvm(x,data=data,na.method=na.method)
    } else {
      dd <- procdata.lvm(x, list(S=S,mu=mu,n=n))
    }
    S <- dd$S; mu <- dd$mu; n <- dd$n    
    var.missing <- setdiff(manifest(x),colnames(S))
  } else { S <- NULL; mu <- NULL }

  if (measurement.fix) {    
    if (length(var.missing)>0) {## Convert to latent:
      new.lat <- setdiff(var.missing,latent(x))
      if (length(new.lat)>0)
      x <- latent(x, new.lat)
    }
    etas <- latent(x)
    ys <- endogenous(x)
    M <- as(Graph(x), Class="matrix")
    for (e in etas) { ## Makes sure that at least one arrow from latent variable is fixed (identification)
      ys. <- names(which(M[e,ys]==1))
      if (length(ys.)>0 & is.na(diag(covfix(x)$values)[e]) )
        if (all(is.na(x$fix[e, ys.]==1))) {
          if (debug)
            cat("Fixing '", ys.[1], "' on '", e, "'.\n", sep="")
          regfix(x,from=e,to=ys.[1]) <- 1
##          x$fix[e,ys.[1]] <- 1
        }
    }
    
    latintNA <- unlist(lapply(intfix(x)[latent(x)],is.na))
    if (length(latintNA)>0) {
      if (any(latintNA))
      intfix(x, latent(x)[which(latintNA)]) <- 0 ## For identifiality we fix mean of latent variables to zero unless already fixed
    }
  }

  if (is.null(S)) x0 <- TRUE
  if (exo.fix) {
    if (x0) {
      S0 <- diag(length(manifest(x)))
      mu0 <- rep(0,nrow(S0))      
    }
    else {
      S0 <- S
      mu0 <- mu
    }
    exo <- exogenous(x);
    exo.idx <- match(exo,manifest(x)); exo_all.idx <- match(exo, vars(x))
    if (length(exo.idx)>0)
      for (i in 1:length(exo.idx))
        for (j in 1:length(exo.idx)) {
          i. <- exo_all.idx[i]; j. <- exo_all.idx[j]
          myval <- S0[exo.idx[i],exo.idx[j]];
          if (i.==j. & myval==0) myval <- 1
          else if (is.na(myval) || is.nan(myval)) myval <- 0
          x$covfix[i.,j.] <- x$covfix[j.,i.] <- myval
        }
      x$mean[exo_all.idx] <- mu0[exo.idx]
  }
  
  index(x) <- reindex(x)
  return(x)
}

###}}}

###{{{ commutation

### Finds the unique commutation matrix:
### K%*%as.vector(A) = as.vector(t(A))
commutation <- function(m, n=m,sparse=FALSE) { 
  H <- function(i,j) { ## mxn-matrix with 1 at (i,j)
    Hij <- matrix(0, nrow=m, ncol=n)
    Hij[i,j] <- 1
    Hij
  }
  K <- matrix(0,m*n,m*n)  
  for (i in 1:m)
    for (j in 1:n)
      K <- K + H(i,j)%x%t(H(i,j))
  if (sparse)
    return(as(K, "sparseMatrix"))
  K  
}

###}}}

###{{{ izero

izero <- function(i,n) { ## n-1 zeros and 1 at ith entry
  x <- rep(0,n); x[i] <- 1
  x
}

###}}}

###{{{ Debug

`Debug` <-
  function(msg, cond) {
    if (cond)
      print(paste(msg, collapse=" "))
  }

###}}}

###{{{ categorical2dummy

categorical2dummy <- function(x,data,silent=TRUE,...) {
  x0 <- x
  X <- exogenous(x)
  catX <- c()  
  for (i in X) {
    if (!is.numeric(data[,i])) catX <- c(catX,i)
  }
  if (length(catX)==0) return(list(x=x,data=data))
  f <- as.formula(paste("~ 1+", paste(catX,collapse="+")))
  opt <- options(na.action="na.pass")
  M <- model.matrix(f,data)
  
  options(opt)
  Mnames <- colnames(M)
  Mpos <- attributes(M)$assign
  A <- index(x)$A
  F <- regfix(x)
  count <- 0
  for (i in catX) {
    count <- count+1
    Y <- colnames(A)[A[i,]==1]
    kill(x0) <- i  
    x0 <- regression(x0,to=Y,from=Mnames[Mpos==count],silent=silent)
  }
  index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
  return(list(x=x0,data=cbind(data,M)))
}

###}}}

###{{{ procdata.lvm

`procdata.lvm` <-
  function(x,data,debug=FALSE,categorical=FALSE,na.method="complete.obs") {    
    if (is.numeric(data) & !is.list(data)) {
      data <- rbind(data)
    }
    if (is.data.frame(data) | is.matrix(data)) {
      nn <- colnames(data)
      data <- as.data.frame(data); colnames(data) <- nn; rownames(data) <- NULL
      obs <- setdiff(intersect(vars(x), colnames(data)),latent(x))
      Debug(obs, debug)
      mydata <- subset(data, select=obs)
      for (i in 1:ncol(mydata)) {
        if (is.Surv(mydata[,i]))
          mydata[,i] <- mydata[,i][,1]
        if (is.character(mydata[,i]) | is.factor(mydata[,i]))
          mydata[,i] <- as.numeric(as.factor(mydata[,i]))-1
      }
      
##      mydata <- data[,obs]
      if (any(is.na(mydata))) {
##        warning("Discovered missing data. Going for a complete-case analysis. For data missing at random see 'missingMLE'.\n", immediate.=TRUE)
        mydata <- na.omit(mydata)
      }
      S <- NULL
      n <- nrow(mydata)
      if (n==1) {
        S <- diag(ncol(mydata)); colnames(S) <- rownames(S) <- obs
      }
      if (na.method=="pairwise.complete.obs") {
        mu <- colMeans(mydata,na.rm=TRUE)
        if (is.null(S))
          S <- cov(na.omit(mydata),use=na.method)        
      }
      if (na.method=="complete.obs") {
        mu <- colMeans(na.omit(mydata))
        if (is.null(S))
          S <- (n-1)/n*cov(na.omit(mydata)) ## MLE variance matrix of observed variables
      }      
    }
    else
      if (is.list(data)) {
        S <- reorderdata.lvm(x,data$S)
        mu <- reorderdata.lvm(x,data$mu)
        n <- data$n
        ##      if (is.null(n)) stop("n was not specified");
      }
      else
        stop("Unexpected type of data!");
    if (nrow(S)!=ncol(S)) stop("Wrong type of data!");
    return(list(S=S,mu=mu,n=n))
  }

###}}}    

###{{{ reorderdata.lvm
`reorderdata.lvm` <-
  function(x, data) {
    if (is.vector(data)) {
      nn <- names(data)
      ii <- na.omit(match(manifest(x), nn))
      data[ii,drop=FALSE]
    } else {
      nn <- colnames(data)
      ii <- na.omit(match(manifest(x), nn))
      data[ii,ii,drop=FALSE]    
    }
  }
###}}}

###{{{ tr

## Trace function. Returns the sum of the diagonal of symmetric matrix
`tr` <-
  function(A) {
    
    ##    if(!is.matrix(A) )
    ##      stop("argument of 'tr' should be a matrix.")    
##     n <- nrow(A)
##         if (!n) 
##       stop("0 x 0 matrix")
##     if (n != ncol(A)) 
##       stop("non-square matrix")
##     if (any(!is.finite(A))) 
##       stop("infinite or missing values")
##     if (any(!is.numeric(A)))
##       stop("numeric or complex values required")
    
    return(sum(diag(A)))
  }

###}}}

###{{{ symmetrize
`symmetrize` <-
function(M, upper=TRUE) {
  if (length(M)==1) return(M)
  if (!is.matrix(M) | ncol(M)!=nrow(M)) stop("Only implemented for square matrices.")
  if (upper) {
    for (i in 1:(ncol(M)-1))
      for (j in (i+1):nrow(M))
        M[i,j] <- M[j,i]
    return(M)
  } else {
    for (i in 1:ncol(M))
      for (j in 1:nrow(M))
        if (M[i,j]==0)          
          M[i,j] <- M[j,i]
        else
          M[j,i] <- M[i,j]
    return(M)
  }
}
###}}}

###{{{ toformula
toformula <- function (y = ".", x = ".") 
{
    xst <- x[1]
    xn <- length(x)
    if (xn > 1) 
        for (i in 2:length(x)) {
            xst <- paste(xst, "+", x[i])
        }
    yst <- y[1]
    yn <- length(y)
    if (yn > 1) {
        yst <- paste("c(", yst, sep = "")
        for (i in 2:length(y)) {
            yst <- paste(yst, ", ", y[i], sep = "")
        }
        yst <- paste(yst, ")", sep = "")
    }
    ff <- paste(yst, "~", xst)
    return(as.formula(ff))
}
###}}} toformula

###{{{ getvars

## getvars <- function(x,env=parent.frame()) {
##   vars <- substitute(x,env=env)
## ##  vars <- eval(substitute(x),env=parent.frame())
##   if (class(vars)[1]=="formula") {
##     vars <- all.vars(vars)
##   }
##   return(as.character(vars))
## }

###}}} getvars

###{{{ frobnorm

# Frobenius norm af matrice x. Hvis x er vektor er dette lig 2-normen
frobnorm <- function(x) {
  frob2 <- tr(t(x)%*%x)
  return(sqrt(frob2))
}
mdist <- function(x,y) { frobnorm(x-y) }
meq <- function(A,tol=1e-9) { frobnorm(A)<tol }

###}}} frobnorm

###{{{ printR

printR <- function(x,eol="\n",...) {
  if (is.vector(x)) {
    row <- paste("c(",paste(x,collapse=","),")",sep="")
    cat(row,eol,sep="")
  }
  if (is.matrix(x)) {
    row <- "rbind("
    cat(row,eol)
    for (i in 1:nrow(x)) {
      printR(x[i,],"")
      cat(ifelse (i<nrow(x),",",""),eol,sep="")      
    }
    row <- ")"
    cat(row,eol)
  }
  invisible(x)
}

###}}}

extractvar <- function(f) {
    yy <- getoutcome(f)
    xx <- attributes(terms(f))$term.labels
    myvars <- all.vars(f)
    return(list(y=yy,x=xx,all=myvars))
}
getoutcome <- function(formula) {
  aa <- attributes(terms(formula))
  if (aa$response==0)
    return(NULL)
  paste(deparse(formula[[2]]),collapse="")
}
decomp.specials <- function(x,pattern="[()]") {
  st <- gsub(" ","",x)
  vars <- rev(unlist(strsplit(st,pattern)))[1]
  unlist(strsplit(vars,","))
}
Decomp.specials <- function(x,pattern="[()]") {
  st <- gsub(" ","",x)
  st <- gsub("\n","",st)
  mysplit <- rev(unlist(strsplit(st,pattern)))
  type <- mysplit[2] 
  vars <- mysplit[1]
  res <- unlist(strsplit(vars,","))
  if (type=="s" | type=="seq") {
    return(paste(res[1],seq(as.numeric(res[2])),sep=""))
  } 
  unlist(strsplit(vars,","))
}

as.numeric.list <- function(x,...) {
  res <- list()
  asnum <- as.numeric(x)
  lapply(x,function(y) ifelse(is.na(as.numeric(y)),y,as.numeric(y)))
}

edge2pair <- function(e) {
  sapply(e,function(x) strsplit(x,"~"))
}
numberdup <- function(xx) { ## Convert to numbered list
  dup.xx <- duplicated(xx)
  dups <- xx[dup.xx]
  xx.new <- numeric(length(xx))
  count <- 0
  for (i in 1:length(xx)) {
    if (!dup.xx[i]) {
      count <- count+1
      xx.new[i] <- count
    } else {
      xx.new[i] <- xx.new[match(xx[i],xx)[1]]
    }
  }
  return(xx.new)
}

logit <- function(p) log(p/(1-p))
tigol <- function(z) 1/(1+exp(-z))

###{{{ Inverse/pseudo

Inverse <- function(X,tol=0,det=TRUE) {
  n <- nrow(X)
  if (nrow(X)==1) {
    res <- 1/X
    if (det) attributes(res)$det <- X
    return(res)
  } 
  svdX <- svd(X)
  id0 <- numeric(n)
  id0[svdX$d>tol] <- 1/svdX$d[svdX$d>tol]
  res <- with(svdX, v%*%diag(id0)%*%t(u))
  if (det)
    attributes(res)$det <- prod(svdX$d)
  return(res)
}

###}}}

###{{{ naiveGrad

naiveGrad <- function(f, x, h=1e-9) {
  nabla <- numeric(length(x))
  for (i in 1:length(x)) {
    xh <- x; xh[i] <- x[i]+h
    nabla[i] <- (f(xh)-f(x))/h
  }
  return(nabla)
}

###}}}

###{{{ whichentry
## X k-dim. array of dimension (d1,d2,...,dk)
## Element x at entry (x1,...,xk).
## position (via which) := x1 + d1*(x2-1) + d1*d2*(x3-1) + ... + prod(d1,...,d[k-1])*(xk-1)
whichentry <- function(x) {
  idx <- which(x)
  D <- dim(x)
  if (length(D)<2)
    return(idx)
  K <- rev(cumprod(D))[-1]
  t(sapply(idx,
           function(ii) {
             entry <- c()
             pn <- ii
             for (i in 1:(length(K))) {
               xn <- floor((pn-1)/K[i])+1
               pn <- pn-(xn-1)*K[i]
               entry <- c(entry,xn)
             }; entry <- rev(c(entry,pn))
             return(entry)
           }))
}

###}}} whichentry

###{{{ blockdiag
`%+%` <-
function(x,...) {
  blockdiag(x,...)
}

blockdiag <- function(x,...,pad=0) {
  xx <- list(x,...)
  rows <- unlist(lapply(xx,nrow))
  crows <- c(0,cumsum(rows))
  cols <- unlist(lapply(xx,ncol))
  ccols <- c(0,cumsum(cols))
  res <- matrix(pad,nrow=sum(rows),ncol=sum(cols))
  for (i in 1:length(xx)) {
    idx1 <- 1:rows[i]+crows[i]; idx2 <- 1:cols[i]+ccols[i]
    res[idx1,idx2] <- xx[[i]]
  }
  colnames(res) <- unlist(lapply(xx,colnames)); rownames(res) <- unlist(lapply(xx,rownames))
  return(res)
}

###}}} blockdiag

###{{{ CondMom

# conditional on Compl(idx)
CondMom <- function(mu,S,idx,X) {
  idxY <- idx
  
  idxX <- setdiff(1:ncol(S),idxY)
  SXX <- S[idxX,idxX,drop=FALSE];
  SYY <- S[idxY,idxY,drop=FALSE]
  SYX <- S[idxY,idxX,drop=FALSE]
  muY <- mu[idxY,drop=FALSE]
  muX <- mu[idxX,drop=FALSE]
  iSXX <- solve(SXX)  
  Z <- apply(X,1,function(xx) xx-muX)
  SZ <- t(SYX%*%iSXX%*%Z)
  condmean <- matrix(apply(SZ,1,function(x) muY+x),ncol=ncol(SZ),nrow=nrow(SZ))
  condvar <- SYY-SYX%*%iSXX%*%t(SYX)
  return(list(mean=condmean,var=condvar))
}

###}}} CondMom

npar.lvm <- function(x) {
  return(index(x)$npar+ index(x)$npar.mean)
}
