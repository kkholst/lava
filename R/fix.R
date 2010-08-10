###{{{ print.fix

print.fix <- function(x,exo=FALSE,...) {
  switch(attributes(x)$type,
        reg = cat("Regression parameters:\n"),
        cov = cat("Regression parameters:\n"),
        mean = cat("Regression parameters:\n"))
  M <- linconstrain(x,print=TRUE)
  ## idx <- 1:attributes(x)$nvar
  ## if (!exo & attributes(x)$type!="reg") idx <- setdiff(idx,attributes(x)$exo.idx)
  ## if (attributes(x)$type=="mean") {
  ##   for (i in idx) {
  ##       cat(names(x)[i],"\t")
  ##   }
  ##   cat("\n")
  ##   for (i in idx) {
  ##       cat(x[[i]],"\t")
  ##   }
  ##   cat("\n")
  ## } else {
  ##   with(x, print(rel[idx,idx,drop=FALSE]))
  ##   with(x, printmany(labels[idx,idx,drop=FALSE], values[idx,idx,drop=FALSE], name1="labels=", name2="values="))
  ## }
  invisible(x)
}


linconstrain <- function(x,print=TRUE,indent="  ",exo=FALSE,...) {  
  idx <- 1:attributes(x)$nvar
  if (!exo & attributes(x)$type!="reg")
    idx <- setdiff(idx,attributes(x)$exo.idx)
  if (attributes(x)$type=="mean") {
    M <- rbind(unlist(x[idx]))
    rownames(M) <- ""
    M[is.na(M)] <- "*"    
  } else {  
    M <- x$rel[idx,idx,drop=FALSE]
    M[M==0] <- NA
    M[M==1] <- "*"
    M[which(!is.na(x$labels[idx,idx]))] <- x$labels[idx,idx][which(!is.na(x$labels[idx,idx]))]
    M[which(!is.na(x$values[idx,idx]))] <- x$values[idx,idx][which(!is.na(x$values[idx,idx]))]
  }
  if (print) {
    M0 <- M
    rownames(M0) <- paste(indent,rownames(M))
    print(M0,quote=FALSE,na.print="",...)
  }
  invisible(M)
}

###}}} print.fix

###{{{ intfix

"intfix" <- function(object,...) UseMethod("intfix")
"intercept" <- function(object,...) UseMethod("intercept")

intercept.lvm <- intfix.lvm <- function(object,...) {
  res <- object$mean; attr(res,"type") <- "mean"
  attr(res,"exo.idx") <- index(object)$exo.idx
  attr(res,"nvar") <- length(res)  
  class(res) <- "fix"
  return(res)
}
"intercept<-" <- "intfix<-" <- function(object,...,value) UseMethod("intfix<-")
"intercept<-.lvm" <- "intfix<-.lvm" <- function(object, vars,...,value) {
  if (class(value)[1]=="formula") {
    lhs <- getoutcome(value)
    yy <- decomp.specials(lhs)
    if ((class(value[[3]])=="logical" && is.na(value[[3]]))) {
      intfix(object,yy) <- NA
      return(object)
    }        
    tt <- terms(value)
    xf <- attributes(terms(tt))$term.labels
    res <- lapply(xf,decomp.specials)[[1]]
    myvalue <- suppressWarnings(as.numeric.list(as.list(res)))
    myvalue <- lapply(myvalue, function(x) ifelse(x=="NA",NA,x))
    intfix(object,yy) <- myvalue
    object$parpos <- NULL
    return(object)
  }
  if (class(vars)[1]=="formula") {
    vars <- all.vars(vars)
  }    
  object$mean[vars] <- value
  newindex <- reindex(object)
  object$parpos <- NULL
  index(object)[names(newindex)] <- newindex
  return(object)
}

###}}} intfix

###{{{ covfix

"covfix" <- function(object,...) UseMethod("covfix")
covfix.lvm <- function(object,...) {
  res <- list(rel=object$cov, labels=object$covpar, values=object$covfix); attr(res,"type") <- "cov"
  attr(res,"exo.idx") <- index(object)$exo.idx
  attr(res,"nvar") <- NROW(res$rel)
  class(res) <- "fix"
  return(res)
}
"covfix<-" <- function(object,...,value) UseMethod("covfix<-")
"covfix<-.lvm" <- function(object, var1, var2=var1, diag=(length(var1)==1),...,value) {

  if (class(var1)[1]=="formula") {
    var1 <- all.vars(var1)
  }    
  if (class(var2)[1]=="formula")
    var2 <- all.vars(var2)  
  
  object <- addvar(object,c(var1,var2),...)
  
  if (is.null(var2)) {
    if (length(value)==1)
      value <- rep(value,length(var1))
    if (length(value)!=length(var1)) stop("Wrong number of parameters")
    for (i in 1:length(var1)) {
      vali <- suppressWarnings(as.numeric(value[[i]]))
      if (is.na(value[[i]]) | value[[i]]=="NA") {
        object$covfix[var1[i],var1[i]] <- object$covpar[var1[i],var1[i]] <- NA
      }
      else {
        if (is.numeric(value[[i]]) | !is.na(vali)) {
          object$covfix[var1[i],var1[i]] <- vali
          object$covpar[var1[i],var1[i]] <- NA
        } else {
          object$covfix[var1[i],var1[i]] <- NA
          object$covpar[var1[i],var1[i]] <- value[[i]]
        }
      }
    }
    newindex <- reindex(object)
    object$parpos <- NULL
    index(object)[names(newindex)] <- newindex    
    return(object)    
  }

  ##object$cov[var1,var2] <-  object$cov[var2,var1] <- 1  
  K <- length(var1)*length(var2)
  if (length(value)==1)
    value <- rep(value,K)
  if (length(value)!=K) stop("Wrong number of parameters")

  for (i in 1:length(var1)) {
    for (j in 1:length(var2)) {
      if (diag | var1[i]!=var2[j]) {
        p <- (i-1)*length(var2) + j
        valp <- suppressWarnings(as.numeric(value[[p]]))
        if (is.na(value[[p]]) | value[[p]]=="NA") {
          object$covfix[var1[i],var2[j]] <- object$covpar[var1[i],var2[j]] <- NA
          object$covfix[var2[j],var1[i]] <- object$covpar[var2[j],var1[i]] <- NA
        }
        else {
          object$cov[var1[i],var2[j]] <-  object$cov[var2[j],var1[i]] <- 1  
          ##        cancel(object) <- c(var1[i],var2[j]) ## Remove old associations
          if (is.numeric(value[[p]]) | !is.na(valp)) {
            object$covfix[var1[i],var2[j]] <- object$covfix[var2[j],var1[i]] <- valp
            object$covpar[var1[i],var2[j]] <- object$covpar[var2[j],var1[i]] <- NA
          } else {
            object$covpar[var1[i],var2[j]] <- object$covpar[var2[j],var1[i]] <- value[[p]]
            object$covfix[var1[i],var2[j]] <- object$covfix[var2[j],var1[i]] <- NA    
          }
        }
      }
    }
  }
  newindex <- reindex(object)
  object$parpos <- NULL
  index(object)[names(newindex)] <- newindex
  return(object)
}

###}}} covfix

###{{{ regfix

"regfix" <- function(object,...) UseMethod("regfix")
regfix.lvm <- function(object,...) {
  res <- list(rel=index(object)$M, labels=object$par, values=object$fix); attr(res,"type") <- "reg"
  attr(res,"exo.idx") <- index(object)$exo.idx
  attr(res,"nvar") <- NROW(res$rel)
  class(res) <- "fix"
  return(res)
}
"regfix<-" <- function(object,...,value) UseMethod("regfix<-")
"regfix<-.lvm" <- function(object, to, from, ..., value) {
  if (is.null(to)) stop("variable list needed")
  if (class(to)[1]=="formula") {
    lhs <- decomp.specials(getoutcome(to))
    if (is.null(lhs)) {
      to <- all.vars(to)
      if (is.null(from)) stop("predictor list needed")
      if (class(from)[1]=="formula")
        from <- all.vars(from)      
    } else {
      from <- all.vars(to)
      to <- lhs ##decomp.specials(lhs)      
      ##      from <- all.vars(extractvar(to)$x
      from <- setdiff(from,to)
    }
  }
  
  object <- addvar(object,c(to,from),...)
  if (length(from)==length(to) & length(from)==length(value)) {
##    if (length(value)!=length(from)) stop("Wrong number of parameters")
    for (i in 1:length(from)) {
      if (!isAdjacent(Graph(object), from[i], to[i])) {
        covfix(object,to[i],from[i]) <- NA## Remove any old correlation specification
        object <- regression(object, to=to, from=from[i])
      }
      vali <- suppressWarnings(as.numeric(value[[i]]))
      if (is.na(value[[i]]) | value[[i]]=="NA") {
        object$fix[from[i],to[i]] <- object$par[from[i],to[i]] <- NA
      }
      else {
        if (is.numeric(value[[i]]) | !is.na(vali)) {
          object$fix[from[i],to[i]] <- vali
          object$par[from[i],to[i]] <- NA
        } else {
          object$par[from[i],to[i]] <- value[[i]]
          object$fix[from[i],to[i]] <- NA
        }      
      }
    }
    newindex <- reindex(object)
    object$parpos <- NULL
    index(object)[names(newindex)] <- newindex
#    index(object) <- reindex(object)
    return(object)
  }
  
  for (i in from) {
    for (j in to) {
      if (!isAdjacent(Graph(object), i, j)) {
##        cancel(object) <- c(i,j) ## Remove old associations
        object <- regression(object,to=j,from=i)
      }
    }
  }
  
  K <- length(from)*length(to)
  if (length(value)==1)
    value <- rep(value,K)  
  if (length(value)!=K) stop("Wrong number of parameters")

  for (j in 1:length(to)) {
    for (i in 1:length(from)) {
      ##      p <- (i-1)*length(to) + j
      p <- (j-1)*length(from) + i
      valp <- suppressWarnings(as.numeric(value[[p]]))
      if (is.na(value[[p]]) | value[[p]]=="NA")
        object$fix[from[i],to[j]] <- object$par[from[i],to[j]] <- NA
      else {
        if (is.numeric(value[[p]]) | !is.na(valp)) {
        object$fix[from[i],to[j]] <- valp
        object$par[from[i],to[j]] <- NA
      } else {
        object$par[from[i],to[j]] <- value[[p]]
        object$fix[from[i],to[j]] <- NA
      }
      }
    }
  }
  newindex <- reindex(object)
  object$parpos <- NULL
  index(object)[names(newindex)] <- newindex
  #index(object) <- reindex(object)
  return(object)
}

###}}} regfix

###{{{ parfix

"parfix<-" <- function(x,...,value) UseMethod("parfix<-")
"parfix<-.lvm" <- function(x,idx,...,value) {
  parfix(x,idx,value,...)
}

"parfix" <- function(x,...) UseMethod("parfix")
parfix.lvm <- function(x,idx,value,fix=TRUE,...) {
  object <- Model(x)
  if (fix)
    object <- fixsome(object)
  if (length(idx)!=length(value))
    value <- rep(value,length.out=length(idx))
  I <- index(object)
  V <- with(I, matrices(Model(object),npar.mean+1:npar,1:npar.mean))
  V$A[I$M0!=1] <- 0; V$P[I$P0!=1] <- 0
  v.fix <- which(V$v%in%idx)
  vval <- V$v[v.fix]
  v.ord <- match(vval,idx)
  Pval <- V$P[V$P%in%idx]
  P.fix <- whichentry(matrix(V$P%in%idx,nrow=nrow(V$P)))
  P.ord <- match(Pval,idx)
  Aval <- V$A[which(V$A%in%idx)]
  A.fix <- whichentry(matrix(V$A%in%idx,nrow=nrow(V$A)))
  A.ord <- match(Aval,idx)
  count <- 0
  if (length(v.fix)) {
    for (i in 1:length(v.fix)) {
      count <- count+1
      object$mean[[v.fix[i]]] <- value[[v.ord[i]]]
    }
  }
  if (length(A.fix)>0) {
    for (i in 1:nrow(A.fix)) {
      count <- count+1
##      if (is.numeric(value[[count]])) {
      if (is.numeric(value[[ A.ord[i] ]])){
        object$fix[A.fix[i,1],A.fix[i,2]] <- value[[A.ord[i]]]
        object$par[A.fix[i,1],A.fix[i,2]] <- NA
      } else {
        object$par[A.fix[i,1],A.fix[i,2]] <- value[[A.ord[i]]]
        object$fix[A.fix[i,1],A.fix[i,2]] <- NA
      }
    }
  }
  if (length(P.fix)>0) {
    for (i in 1:nrow(P.fix)) {
      count <- count+1
      ##      if (is.numeric(value[[count]])) {
      if (is.numeric(value[[ P.ord[i] ]])) {
        object$covfix[P.fix[i,1],P.fix[i,2]] <- value[[P.ord[i]]]
        object$covpar[P.fix[i,1],P.fix[i,2]] <- NA
      } else {
        object$covpar[P.fix[i,1],P.fix[i,2]] <- value[[P.ord[i]]]
        object$covfix[P.fix[i,1],P.fix[i,2]] <- NA
      }
    }
  }
  newindex <- reindex(object)
  object$parpos <- NULL
  index(object)[names(newindex)] <- newindex
  attributes(object)$fixed <- list(v=v.fix,A=A.fix,P=P.fix)
  return(object)
}

###}}} parfix
