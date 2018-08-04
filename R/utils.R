char2num <- function(x,...) {
    idx <- grep("^[-]*[0-9\\.]+",x,perl=TRUE,invert=TRUE)
    if (length(idx)>0) x[idx] <- NA
    as.numeric(x)
}

###{{{ substArg

substArg <- function(x,env,...) {
  if (!missing(env)) {
    a <- with(env,substitute(x))
  } else {
    a <- substitute(x)
  }
  myclass <- tryCatch(class(eval(a)),error=function(e) NULL)
  if (is.null(myclass) || myclass=="name") {
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
  for (i in seq_len(object$ngroup)) {
    x0 <- object$lvm[[i]]
    data0 <- data[[i]]
    xfix0 <- colnames(data0)[(colnames(data0)%in%parlabels(x0,exo=TRUE))]
    xfix <- c(xfix, list(xfix0))
    if (length(xfix0)>0) { ## Yes, random slopes
      Xfix<-TRUE
    }
    xx$lvm[[i]] <- x0
  }
  if (Xfix) {
    for (k in seq_len(object$ngroup)) {
      x0 <- object$lvm[[k]]
      data0 <- data[[k]]
      nrow <- length(vars(x0))
      xpos <- lapply(xfix[[k]],function(y) which(regfix(x0)$labels==y))
      colpos <- lapply(xpos, function(y) ceiling(y/nrow))
      rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
      myfix0 <- list(var=xfix[[k]], col=colpos, row=rowpos)
      myfix <- c(myfix, list(myfix0))
      for (i in seq_along(myfix0$var))
        for (j in seq_along(myfix0$col[[i]]))
          regfix(x0,
                 from=vars(x0)[myfix0$row[[i]][j]],to=vars(x0)[myfix0$col[[i]][j]]) <-
                   colMeans(data0[,myfix0$var[[i]],drop=FALSE],na.rm=TRUE)
      index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
      object$lvm[[k]] <- x0
    }
    object <- multigroup(object$lvm,data,fix=FALSE,exo.fix=FALSE)
  }
  return(list(model=object,fix=myfix))
}

###}}} procrandomslope

###{{{ kronprod

## ' Calculate matrix product with kronecker product
## '
## ' \deqn{(A\crossprod B) Y}
## ' @title Calculate matrix product with kronecker product
## ' @param A
## ' @param B
## ' @param Y
## ' @author Klaus K. Holst
kronprod <- function(A,B,Y) {
    if (missing(Y)) {
        ## Assume 'B'=Identity, (A otimes B)Y
        k <- nrow(B)/ncol(A)
        res <- rbind(apply(B,2,function(x) matrix(x,nrow=k)%*%t(A)))
        return(res)
    }
    rbind(apply(Y,2,function(x) B%*%matrix(x,nrow=ncol(B))%*%t(A)))
}

###}}} kronprod

###{{{ izero

izero <- function(i,n) { ## n-1 zeros and 1 at ith entry
  x <- rep(0,n); x[i] <- 1
  x
}

###}}}

###{{{ Debug

`Debug` <-
  function(msg, cond=lava.options()$debug) {
    if (cond)
      print(paste(msg, collapse=" "))
  }

###}}}

###{{{ categorical2dummy

categorical2dummy <- function(x,data,messages=0,...) {
  x0 <- x
  X <- intersect(index(x)$exogenous,colnames(data))
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
    mnames <- Mnames[Mpos==count]
    kill(x0) <- i
    Y <- colnames(A)[A[i,]==1]
    if (length(mnames)==1) {
      fix <- as.list(F$labels[i,])
      fixval <- F$values[i,]
      fix[which(!is.na(fixval))] <- fixval[na.omit(fixval)]
      regression(x0,to=Y,from=mnames,messages=messages) <- fix[Y]
    } else {
      x0 <- regression(x0,to=Y,from=mnames,messages=messages)
    }
  }
  index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
  return(list(x=x0,data=cbind(data,M)))
}

###}}}

###{{{ procdata.lvm

`procdata.lvm` <-
  function(x,data,categorical=FALSE,
    na.method=ifelse(any(is.na(data[,intersect(colnames(data),manifest(x))])),"complete.obs","pairwise.complete.obs"),
    missing=FALSE
    ) {
    if (is.numeric(data) & !is.list(data)) {
      data <- rbind(data)
    }
     if (is.data.frame(data) | is.matrix(data)) {
      nn <- colnames(data)
      data <- as.data.frame(data); colnames(data) <- nn; rownames(data) <- NULL
      obs <- setdiff(intersect(vars(x), colnames(data)),latent(x))
      Debug(obs)
      mydata <- subset(data, select=obs)
      if (NROW(mydata)==0) stop("No observations")
      for (i in seq_len(ncol(mydata))) {
        if (inherits(mydata[,i],"Surv"))
          mydata[,i] <- mydata[,i][,1]
        if (is.character(mydata[,i]) | is.factor(mydata[,i]))
          mydata[,i] <- as.numeric(as.factor(mydata[,i]))-1
      }
      
      S <- NULL
      n <- nrow(mydata)
      if (n==1) {
        S <- diag(nrow=ncol(mydata)); colnames(S) <- rownames(S) <- obs
      }
      if (na.method=="complete.obs" && !missing) {
        mydata0 <- na.omit(mydata)
        n <- nrow(mydata0)
        mu <- colMeans(mydata0)
        if (is.null(S) && n>2) 
            S <- (n-1)/n*cov(mydata0) ## MLE variance matrix of observed variables
        rm(mydata0)
      }
      nS <- is.null(S) || any(is.na(S))
      if (na.method=="pairwise.complete.obs" || nS) {
          mu <- colMeans(mydata,na.rm=TRUE)
          if (nS) {
              n <- nrow(mydata)
              S <- (n-1)/n*cov(mydata,use="pairwise.complete.obs")
              S[is.na(S)] <- 1e-3
          }
      }
    }
    else
      if (is.list(data)) {
        if ("cov"%in%names(data)) data$S <- data$cov
        if ("var"%in%names(data)) data$S <- data$var
        if ("mean"%in%names(data)) data$mu <- data$mean
        n <- data$n
        S <- reorderdata.lvm(x,data$S)
        mu <- reorderdata.lvm(x,data$mu)
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
      ii <- na.omit(match(index(x)$manifest, nn))
      data[ii,drop=FALSE]
    } else {
      nn <- colnames(data)
      ii <- na.omit(match(index(x)$manifest, nn))
      data[ii,ii,drop=FALSE]
    }
  }

###}}}

###{{{ symmetrize

`symmetrize` <-
function(M, upper=TRUE) {
  if (length(M)==1) return(M)
  if (!is.matrix(M) | ncol(M)!=nrow(M)) stop("Only implemented for square matrices.")
  if (upper) {
    for (i in seq_len(ncol(M)-1))
      for (j in seq(i+1,nrow(M)))
        M[i,j] <- M[j,i]
    return(M)
  } else {
    for (i in seq_len(ncol(M)))
      for (j in seq_len(nrow(M)))
        if (M[i,j]==0)
          M[i,j] <- M[j,i]
        else
          M[j,i] <- M[i,j]
    return(M)
  }
}

###}}}

###{{{ naiveGrad

naiveGrad <- function(f, x, h=1e-9) {
  nabla <- numeric(length(x))
  for (i in seq_along(x)) {
    xh <- x; xh[i] <- x[i]+h
    nabla[i] <- (f(xh)-f(x))/h
  }
  return(nabla)
}

###}}}

###{{{ CondMom

# conditional on Compl(idx)
CondMom <- function(mu,S,idx,X) {
  idxY <- idx

  idxX <- setdiff(seq_len(ncol(S)),idxY)
  SXX <- S[idxX,idxX,drop=FALSE];
  SYY <- S[idxY,idxY,drop=FALSE]
  SYX <- S[idxY,idxX,drop=FALSE]
  iSXX <- solve(SXX)
  condvar <- SYY-SYX%*%iSXX%*%t(SYX)
  if (missing(mu)) return(condvar)

  muY <- mu[,idxY,drop=FALSE]
  muX <- mu[,idxX,drop=FALSE]
  if (is.matrix(mu))
    Z <- t(X-muX)
  else
    Z <- apply(X,1,function(xx) xx-muX)
  SZ  <- t(SYX%*%iSXX%*%Z)
##  condmean <- matrix(
  if (is.matrix(mu))
    condmean <- SZ+muY
  else
    condmean <- t(apply(SZ,1,function(x) muY+x))
##  ,ncol=ncol(SZ),nrow=nrow(SZ))
  return(list(mean=condmean,var=condvar))
}

###}}} CondMom

###{{{ Depth-First/acc (accessible)

DFS <- function(M,v,explored=c()) {
  explored <- union(explored,v)
  incident <- M[v,]
  for (v1 in setdiff(which(incident==1),explored)) {
    explored <- DFS(M,v1,explored)
  }
  return(explored)
}

acc <- function(M,v) {
  if (is.character(v)) v <- which(colnames(M)==v)
  colnames(M)[setdiff(DFS(M,v),v)]
}

###}}} Depth-First/acc (accessible)


npar.lvm <- function(x) {
  return(index(x)$npar+ index(x)$npar.mean+index(x)$npar.ex)

}

as.numeric.list <- function(x,...) {
  lapply(x,function(y) ifelse(is.na(as.numeric(y)),y,as.numeric(y)))
}

edge2pair <- function(e) {
  sapply(e,function(x) strsplit(x,"~"))
}
numberdup <- function(xx) { ## Convert to numbered list
  dup.xx <- duplicated(xx)
  ## dups <- xx[dup.xx]
  xx.new <- numeric(length(xx))
  count <- 0
  for (i in seq_along(xx)) {
    if (!dup.xx[i]) {
      count <- count+1
      xx.new[i] <- count
    } else {
      xx.new[i] <- xx.new[match(xx[i],xx)[1]]
    }
  }
  return(xx.new)
}

extractvar <- function(f) {
    yy <- getoutcome(f)
    xx <- attributes(terms(f))$term.labels
    myvars <- all.vars(f)
    return(list(y=yy,x=xx,all=myvars))
}

##' @export
getoutcome <- function(formula,sep,...) {
  aa <- attributes(terms(formula,...))
  if (aa$response==0) {
    res <- NULL
  } else {
    res <- paste(deparse(formula[[2]]),collapse="")
  }
  if (!missing(sep) && length(aa$term.labels)>0) {
      attributes(res)$x <- lapply(strsplit(aa$term.labels,"\\|")[[1]],
                                  function(x) as.formula(paste0("~",x)))
  } else {
      attributes(res)$x <- aa$term.labels
  }
  return(res)
}


##' @export
Specials <- function(f,spec,split2="+",...) {
  tt <- terms(f,spec)
  pos <- attributes(tt)$specials[[spec]]
  if (is.null(pos)) return(NULL)
  x <- rownames(attributes(tt)$factors)[pos]
  st <- gsub(" ","",x)
  res <- unlist(strsplit(st,"[()]"))[2]
  if (is.null(split2)) return(res)
  unlist(strsplit(res,"+",fixed=TRUE))
}


##' @export
decomp.specials <- function(x,pattern="[()]",pattern2=NULL, pattern.ignore=NULL, sep="[,\\+]",perl=TRUE,reverse=FALSE,...) {
  st <- gsub(" |^\\(|)$","",x) # Remove white space and leading/trailing parantheses
  if (!is.null(pattern.ignore)) {
      if (grepl(pattern.ignore,st,perl=perl,...)) return(st)
  }
  if (!is.null(pattern)) {
    st <- rev(unlist(strsplit(st,pattern,perl=perl,...)))[1]
  }
  if (!is.null(pattern2)) {
    st <- (unlist(strsplit(st,pattern2,perl=perl,...)))
    if (reverse) st <- rev(st)
  }
  unlist(strsplit(st,sep,perl=perl,...))
}

Decomp.specials <- function(x,pattern="[()]") {
  st <- gsub(" ","",x)
  st <- gsub("\n","",st)
  mysplit <- rev(unlist(strsplit(st,pattern)))
  type <- mysplit[2]
  vars <- mysplit[1]
  res <- unlist(strsplit(vars,","))
  if (type=="s" | type=="seq") {
    return(paste0(res[1],seq(char2num(res[2]))))
  }
  unlist(strsplit(vars,","))

}

printline <- function(n=70) {
    cat(rep("_", n), "\n", sep="");

}
