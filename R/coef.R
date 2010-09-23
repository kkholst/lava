###{{{ coef.lvm

`coef.lvm` <-
function(object, mean=TRUE, fix=TRUE, symbol=c("<-","<->"," on "," with "), silent=TRUE, p, data, level=9, labels=FALSE, debug=FALSE,...) {
 if (fix) ## 12/7-2010
   object <- fixsome(object,measurement.fix=FALSE)
  if (!missing(p)) {
    coefs <- matrix(NA,nrow=length(p),ncol=4); coefs[,1] <- p
    rownames(coefs) <- c(coef(object,mean=TRUE)[c(1:index(object)$npar.mean)],paste("p",1:index(object)$npar,sep=""))
    if (!is.null(data)) {
      I <- information(object,p=p,data=data,type="E")
      myvcov <- solve(I)
    } else {
      myvcov <- matrix(NA,length(p),length(p))
    }
    object$vcov <- myvcov
    coefs[,2] <- sqrt(diag(myvcov))
    coefs[,3] <- coefs[,1]/coefs[,2]
    coefs[,4] <-  2*(1-pnorm(abs(coefs[,3])))    
    colnames(coefs) <- c("Estimate","Std. Error", "Z value", "Pr(>|z|)")
    object$coefficients <- coefs;
    return(coef.lvmfit(object,level=level,...))
  }  

  AP <- matrices(object, paste("p",1:index(object)$npar, sep=""))
  
  A <- AP$A; A[index(object)$M1==0] <- "0" ## Only free parameters
  P <- AP$P; P[index(object)$P1==0] <- "0"
  ##  P[lower.tri(P, diag=FALSE)] <- "0"
  nn <- vars(object)
  
  counter <- 0
  res <- c()
  resname <- c()
  for (i in 1:ncol(A))
    for (j in 1:nrow(A)) {
      val <- A[j,i]
      if (val!="0") {
        if (labels & !is.na(regfix(Model(object))$labels[j,i]))
          res <- c(res, regfix(Model(object))$labels[j,i])
        else
          res <- c(res, paste(nn[i],symbol[1],nn[j],sep=""))
        counter <- counter+1
        resname <- c(resname, val)      
      }
    }

  for (i in 1:ncol(P))
    for (j in i:nrow(P))
    {
      val <- P[j,i]
      if (val!="0") {
        counter <- counter+1
        if (labels & !is.na(covfix(Model(object))$labels[j,i]))
          res <- c(res, covfix(Model(object))$labels[j,i])
        else
          res <- c(res, paste(nn[i],symbol[2],nn[j],sep=""))
        resname <- c(resname, val)
      }
    }
  
  names(res) <- resname
  resnum <- sapply(resname, function(s) as.numeric(substr(s,2,nchar(s))))
  res <- res[order(resnum)]
 if (mean) {
   nmean <- sum(index(object)$v1==1)
   if (nmean>0) {

     if (!labels)
       res <- c(vars(object)[index(object)$v1==1], res)
     else {
       mres <- c()
       for (i in 1:length(index(object)$v1)) {
       val <- index(object)$v1[i]
       if (val==1) {
         if (!is.na(intfix(Model(object))[[i]])) {
           mres <- c(mres, intfix(Model(object))[[i]])
         }
         else
           mres <- c(mres, vars(object)[i])
       }
     }
       res <- c(mres,res)     
     }
     names(res)[1:nmean] <- paste("m",1:nmean,sep="")
   }
 }
 
  if (!silent) {
    cat(paste(res, collapse="\n")); cat("\n")
  }
  res
}

###}}}

###{{{ coef.lvmfit

`coef.lvmfit` <-
function(object, level=ifelse(missing(type),-1,2), symbol=c("<-","<->"," on "," with "), data,
         std=NULL, labels=TRUE, vcov, 
         debug=FALSE, type, reliability=FALSE, second=FALSE, ...) {

  ## object0 <- object
  ## xfix <- colnames(model.frame(object))[(colnames(model.frame(object))%in%parlabels(Model(object)))]
  ## if (missing(fix)) {
  ##   fix <- ifelse(length(xfix)>0,FALSE,TRUE)
  ## }
  if (is.null(object$control$meanstructure)) meanstructure <- TRUE
  else meanstructure <- object$control$meanstructure
  npar <- index(object)$npar; npar.mean <- index(object)$npar.mean*meanstructure

  if (class(object)[1]=="lvm.missing") {
    ##npar <- object$estimate$model$npar; npar.mean <- object$estimate$model$npar.mean    ##    myorder <- modelPar(object$estimate$model,1:(npar+npar.mean))$p[[object$cc]]
    ##    myorder.reg <- modelPar(object$estimate$model,1:(npar))$p[[object$cc]]
    if (length(object$cc)==0) {## No complete cases
      coefs <- coef(object$estimate)
      p <- pars(object)
      pn <- names(p)
      pp <- modelPar(object$multigroup,pn)$p
      coefnames <- c()
      for (i in 1:length(pars(object))) {
        for (j in 1:length(pp)) {
          idx <- which(pn[i]==pp[[j]])
          if (length(idx)>0) {
            coefnames <- c(coefnames, rownames(coefs[[j]])[idx[1]])
            break;
          }
        }
      }
      c1 <- coef(Model(object),mean=TRUE)
      c1. <- coef(Model(object),mean=FALSE)      
      ##      varpar1 <- which(sapply(c1,function(x) length(grep("<->",x)))==1)
      ##      varpar2 <- which(sapply(coefnames,function(x) length(grep("<->",x)))==1)
      myorder <- match(coefs,c1)
      myorder.reg <- na.omit(match(coefs,c1.))
##      myorder <- match(c2,c1)
##      myorder.reg <- na.omit(match(c2,c1.))
      
##      myorder <- match(c1,c2)
##      myorder.reg <- match(c1.,c2)-(length(c1)-length(c1.))
    } else {
      myorder <- modelPar(object$multigroup,1:(npar+npar.mean))$p[[object$cc]]
      myorder.reg <- modelPar(object$multigroup,1:(npar))$p[[object$cc]]
    }
  } else {
    myorder <- 1:(npar+npar.mean)
    myorder.reg <- 1:npar
  }
  
  if (level<0) {
    res <- pars.default(object)
    names(res) <- coef(Model(object), mean=meanstructure)[myorder]
    return(res)
  }
  latent.var <- latent(object)
  latent.idx <- which(vars(object)%in%latent.var)
  Type <- Var <- From <- VarType <- FromType <- c()

  

  Astd <- Pstd <- vstd <- mytype <- NULL
  if (!is.null(std)) {
    stdCoef <- stdcoef(object)
    {
      switch(tolower(std),
             latent = {Astd=stdCoef$Astar; Pstd=stdCoef$Pstar; vstd=stdCoef$vstar},
             y = {Astd=stdCoef$AstarY; Pstd=stdCoef$PstarY; vstd=stdCoef$vstarY},
             xy = {Astd=stdCoef$AstarXY; Pstd=stdCoef$PstarXY; vstd=stdCoef$vstarXY},
             yx = {Astd=stdCoef$AstarXY; Pstd=stdCoef$PstarXY; vstd=stdCoef$vstarXY}
             )
    }
  }
  
  myparnames <- paste("p",1:npar,sep="")[myorder.reg]
  myparnames <- paste("p",1:npar,sep="")[order(myorder.reg)]
  
  p <- matrices(Model(object), myparnames)
  A <- p$A
  P <- p$P
  ## coefs <- object$coef[myorder,]
  mycoef <- object$coef
  if (!missing(type) | !missing(vcov)) {
    if (!missing(vcov)) {
      mycoef[,2] <- sqrt(diag(vcov))
    } else {
      if (!missing(data)) 
        I <- informaton(object,type=type,data=data)
      else
        I <- information(object,type=type)    
      myvcov <- solve(I)
      mycoef[,2] <- sqrt(diag(myvcov))
    }
    mycoef[,3] <- mycoef[,1]/mycoef[,2]
    mycoef[,4] <-  2*(1-pnorm(abs(mycoef[,3])))
  }
  coefs <- mycoef[order(myorder),,drop=FALSE]
  nn <- colnames(A)
  
  free <- A!="0" 
  free[index(object)$M1!=1] <- FALSE
  nlincon <- matrix(Model(object)$par%in%names(constrain(Model(object))),nrow(A))
  if (missing(data)) {
    data <- matrix(0,ncol=length(manifest(object))); colnames(data) <- manifest(object)
  }
  nlincon.estimates.full<- constraints(object,second=second,data=data)  
  nlincon.estimates <- nlincon.estimates.full[,-(5:6),drop=FALSE]
  matched <- c()
  res <- c()
  
  for (i in 1:ncol(A))
    for (j in 1:nrow(A)) {
      val <- A[j,i]
      if (val!="0") {        
        matching <- match(val,rownames(coefs))
        matched <- c(matched,matching)
        if (!is.na(matching)) {
          if (free[j,i])
            newrow <- matrix(coefs[matching,],nrow=1)
          else {
            newrow <- matrix(c(coefs[matching,1],NA,NA,NA), nrow=1)
          }
        } else {
          Debug(list("(i,j)", i, ",", j),debug)
          if (nlincon[j,i]) {
            newrow <- matrix(nlincon.estimates[Model(object)$par[j,i],],nrow=1)
          } else {
            newrow <- matrix(c(Model(object)$fix[j,i], NA, NA, NA), nrow=1)
          }
        }
        if (!is.null(std)) {
          newrow <- cbind(newrow,Astd[j,i])
        }
        if (labels & !is.na(regfix(Model(object))$labels[j,i])) {
          rownames(newrow) <- regfix(Model(object))$labels[j,i]
          if (labels>1) {
            rownames(newrow) <- paste(rownames(newrow),paste(nn[i],symbol[1],nn[j],sep=""),sep=":")
          }
        } else {       
          rownames(newrow) <- paste(nn[i],symbol[1],nn[j],sep="")
        }
##        if (is.na(Model(object)$fix[j,i]) | level>0) {
        if (free[j,i] | level>2) {
          res <- rbind(res, newrow)
          Type <- c(Type,"regression")
          Var <- c(Var, nn[i])
          From <- c(From, nn[j])
        }
      }
    }
  free.var <- P!="0" 
  free.var[index(object)$P1!=1] <- FALSE
  nlincon.var <- matrix(Model(object)$covpar%in%names(constrain(Model(object))),nrow(P))  


    if (level>0)   
      ## Variance estimates:
      for (i in 1:ncol(p$P))
        for (j in i:nrow(p$P)) {
          val <- p$P[j,i]
          if (val!="0" & !any(vars(object)[c(i,j)]%in%exogenous(object)))
            if (level>1 | !all(vars(object)[c(i,j)]%in%manifest(object)))
            {
            matching <- match(val,rownames(coefs))
            matched <- c(matched,matching)
            
            if (!is.na(matching)) {
              if (free.var[j,i])
                newrow <- matrix(coefs[matching,],nrow=1)
              else
                newrow <- matrix(c(coefs[matching,1],NA,NA,NA), nrow=1)
                ## We don't want to report p-values of tests on the boundary of the parameter space
              if (i==j)
                newrow[,4] <- NA
            } else {
              Debug(list("(i,j)", i, ",", j),debug)
              if (nlincon.var[j,i]) {
                newrow <- matrix(nlincon.estimates[Model(object)$covpar[j,i],],nrow=1)
              } else {
                newrow <- matrix(c(Model(object)$covfix[j,i], NA, NA, NA), nrow=1)
              }
            }               
###            if (is.na(Model(object)$covfix[j,i]) | level>0)
            if (!missing(std)) {
              newrow <- cbind(newrow,Pstd[i,j])
            }
            if (labels & !is.na(covfix(Model(object))$labels[j,i])) {
              rownames(newrow) <- covfix(Model(object))$labels[j,i]
              if (labels>1) {
                rownames(newrow) <- paste(rownames(newrow),paste(nn[i],symbol[2],nn[j],sep=""),sep=":")
              }
            } else {       
              rownames(newrow) <- paste(nn[i],symbol[2],nn[j],sep="")
            }
            if ((free.var[j,i]) | level>2) {
              res <- rbind(res, newrow)
              Type <- c(Type,"variance")
              Var <- c(Var, nn[i])
              From <- c(From, nn[j])
            }
          }
        }
  res0 <- res

  ##  res <- res0
  ## Mean parameter:
  nlincon.mean <- lapply(Model(object)$mean, function(x) x%in%names(constrain(Model(object))) )

  if (level>0 & npar.mean>0) {
    mu.estimated <- coefs[1:index(object)$npar.mean,,drop=FALSE] ##coefs[setdiff(1:nrow(coefs),matched),]
    munames <- rownames(coefs)[1:index(object)$npar.mean]
    meanpar <- matrices(Model(object), myparnames, munames)$v
    for (i in 1:length(meanpar)) {
      if (!vars(object)[i]%in%exogenous(object)) {
        val <- meanpar[i]
        matching <- match(val,rownames(coefs))
        if (!is.na(matching)) {
          if (index(object)$v1[i]==1)  ## if free-parameter
            newrow <- matrix(coefs[matching,],nrow=1)
          else
            newrow <- matrix(c(coefs[matching,1],NA,NA,NA), nrow=1)
        } else {
          if (nlincon.mean[[i]]) {
            newrow <- matrix(nlincon.estimates[Model(object)$mean[[i]],],nrow=1)
          } else {
            newrow <- matrix(c(as.numeric(meanpar[i]), NA, NA, NA), nrow=1)          }
        }
        if (!missing(std)) {
          newrow <- cbind(newrow,vstd[i])
        }
        if (labels & !(is.na(intfix(Model(object))[[i]]) | is.numeric(intfix(Model(object))[[i]]))) {
          rownames(newrow) <- intfix(Model(object))[[i]]
          if (labels>1) {
            rownames(newrow) <- paste(rownames(newrow),vars(Model(object))[i],sep=":")
          }          
        } else {       
          rownames(newrow) <- vars(Model(object))[i]
        }        
        if ((index(object)$v1[i]==1) | level>2) {
          res <- rbind(res, newrow)
          Type <- c(Type,"intercept")
          Var <- c(Var, vars(Model(object))[i])
          From <- c(From, NA)
        }
      }
    }
  }

  mycolnames <- colnames(coefs)
  if (!is.null(std)) mycolnames <- c(mycolnames, paste("std",std,sep="."))
  colnames(res) <- mycolnames
  attributes(res)$type <- Type
  attributes(res)$var <- Var
  attributes(res)$from <- From
  attributes(res)$latent <- latent.var
  attributes(res)$nlincon <- nlincon.estimates.full

  return(res)
}

###}}} coef.lvmfit

###{{{ coef.multigroup

coef.multigroup <- function(object,...) {
  return(object$parpos)
}

###}}} coef.multigroup

###{{{ coef.multigroupfit

coef.multigroupfit <-
  function(object, level=1,vcov,
           debug=FALSE,labels=FALSE,...) {

    if (level==0) {
      theta <- pars(object)
      if (missing(vcov))
        theta.sd <- sqrt(diag(object$vcov))
      else
        theta.sd <- sqrt(diag(vcov))
      res <- cbind(theta,theta.sd,(Z <- theta/theta.sd),2*(1-pnorm(abs(Z))))
      colnames(res) <- c("Estimate","Std. Error", "Z value", "Pr(>|z|)")
      return(res)
    }
    
##     cc <- coef(object, level=0)
##     object0 <- object
##     object0$model$lvm <- object$model$models.orig
##     object0$model$ngroup <- length(object0$model$lvm)
##     parpos <- modelPar(Model(object0), 1:nrow(cc))$p
    
##     cc <- coef(object0, level=0)

##     model <- Model(object0)
##     parpos <- modelPar(model, 1:nrow(cc))$p
##     res <- list()
##     for (i in 1:model$ngroup) {
##       newcoef <- cc[parpos[[i]],,drop=FALSE]
##       rownames(newcoef) <- coef(Model(model)[[i]], mean=object$meanstructure, silent=TRUE)
##       ## Position of variance parameters:
##       varpos <- variances(Model(model)[[i]])
##       ## Number of parameters resp mean-parameters
##       p <- nrow(newcoef); p0 <- length(coef(Model(model)[[i]], mean=FALSE, silent=TRUE))
##       newcoef[(p-p0) + varpos,4] <- NA
##       res <- c(res, list(newcoef))
##     }    

    
    cc <- coef(object, level=0, ...)
    model <- Model(object)
    parpos <- modelPar(model, 1:nrow(cc))$p
    npar.mean <- object$model$npar.mean
    npar <- object$model$npar
    mynames <- c()
    if (npar.mean>0) {
      mynames <- unlist(object$model$meanlist)
      mynames <- names(mynames)[!duplicated(mynames)]
    }
    if (npar>0) {
      mynames <- c(mynames,object$model$par)
    }

    res <- list()
    for (i in 1:model$ngroup) {
      newcoef <- cc[parpos[[i]],,drop=FALSE]
  ##    rownames(newcoef) <- mynames[parpos[[i]]]      
      rownames(newcoef) <- orignames <- coef(Model(model)[[i]], mean=object$meanstructure, silent=TRUE)
##      attributes(newcoef)$pnames <- orignames
##      attributes(newcoef)$names <- coef(Model(model)[[i]], mean=object$meanstructure, silent=TRUE)
##      rownames(newcoef) <- coef(Model(model)[[i]], mean=object$meanstructure, silent=TRUE)
      ## Position of variance parameters:
      varpos <- variances(Model(model)[[i]],mean=FALSE)
      ## Number of parameters resp mean-parameters
      p <- nrow(newcoef); p0 <- length(coef(Model(model)[[i]], mean=FALSE, silent=TRUE))
      newcoef[(p-p0) + varpos,4] <- NA
      res <- c(res, list(newcoef))
    }
    return(res)
}

###}}}

###{{{ CoefMat

CoefMat.multigroupfit <- function(x,level=9,labels=FALSE,...) {
  cc <- coef(x,level=level)
  parpos <- modelPar(x, 1:length(pars(x)))$p
  res <- c()
  nlincon.estimates <- c()
  nlincon.names <- c()
  for (i in 1:Model(x)$ngroup) {
    m0 <- Model(Model(x))[[i]]
    mycoef <- cc[[i]]
    npar <- index(m0)$npar
    npar.mean <- index(m0)$npar.mean
    if (npar>0)
      rownames(mycoef)[(1:npar)+npar.mean] <- paste("p",1:npar,sep="")    
    m0$coefficients <- mycoef
    Vcov <- vcov(x)[parpos[[i]],parpos[[i]],drop=FALSE]; colnames(Vcov) <- rownames(Vcov) <- rownames(mycoef)
    m0$vcov <- Vcov
    cc0 <- coef.lvmfit(m0,level=level,labels=labels)
    res <- c(res, list(CoefMat(cc0)))
    newnlin <- attributes(cc0)$nlincon
    if (length(newnlin)>0)
    if (i==1) {
      nlincon.estimates <- newnlin
      nlincon.names <- rownames(newnlin)
    } else {
      for (j in 1:NROW(newnlin)) {
        if (!(rownames(newnlin)[j]%in%nlincon.names)) {
          nlincon.estimates <- rbind(nlincon.estimates,newnlin[j,])
          nlincon.names <- c(nlincon.names,rownames(newnlin))
        }
      }
    }
  }
  rownames(nlincon.estimates) <- nlincon.names
  attributes(res)$nlincon <- nlincon.estimates
  return(res)
}

CoefMat <- function(x,digits=5,scientific=0,level=9,...) {
  cc <- x
  if (!is.matrix(x)) {
    cc <- coef(x,level=level,...)
  }
  res <- c()
  mycoef <- format(cc,digits=digits,scientific=scientific,...)
##  mycoef <- formatC(cc,digits=digits,scientific=scientific,...)
  mycoef[is.na(cc)] <- ""
  mycoef[cc[,4]<1e-16,4] <- "  <2e-16"
  

  M <- ncol(cc)
  N <- nrow(cc)
  Nreg <- sum(attributes(cc)$type=="regression")
  Nvar <- sum(attributes(cc)$type=="variance")
  Nint <- sum(attributes(cc)$type=="intercept")
  latent.var <- attributes(cc)$latent
  
  if (Nreg>0) {
    reg.idx <- which(attributes(cc)$type=="regression")
    latent.from <- which(attributes(cc)$from[reg.idx]%in%latent.var)
    latent.from <- latent.from[which(is.na(match(attributes(cc)$var[latent.from],latent.var)))]
    
    reg.idx <- setdiff(reg.idx,latent.from)
    Nmeas <- length(latent.from)
    if (Nmeas>0) {
      first.entry <- c()
      for (i in latent.var) {
        pos <- match(i,attributes(cc)$from[latent.from])
        if (!is.na(pos))
          first.entry <- c(first.entry, pos)
      }
      res <- rbind(res, c("Measurements:",rep("",M)))
      count <- 0
      Delta <- FALSE
      for (i in latent.from) {
        count <- count+1
        newrow <- mycoef[i,]        
        newname <- rownames(cc)[i]
        if (count%in%first.entry) Delta <- !Delta
        prefix <- ifelse(Delta,"  ","   ")
        res <- rbind(res,c(paste(prefix,newname),newrow))
        ##    res <- rbind(res,c(newname,newrow))
      }      
    }    
    if ((Nreg-Nmeas)>0) {
      responses <- unique(attributes(cc)$var[reg.idx])
      first.entry <- c()
      for (i in responses) {
        pos <- match(i,attributes(cc)$var[reg.idx])
        first.entry <- c(first.entry, pos)
      }
      res <- rbind(res, c("Regressions:",rep("",M)))
      count <- 0
      Delta <- FALSE
      for (i in reg.idx) {
        count <- count+1
        newrow <- mycoef[i,]
        newname <- rownames(cc)[i]
        if (count%in%first.entry) Delta <- !Delta
        prefix <- ifelse(Delta,"  ","   ")
        res <- rbind(res,c(paste(prefix,newname),newrow))
        ##    res <- rbind(res,c(newname,newrow))
      }
    }
  }
  if (Nint>0) {
    int.idx <- which(attributes(cc)$type=="intercept")
    res <- rbind(res, c("Intercepts:",rep("",M)))
    for (i in int.idx) {
      newrow <- mycoef[i,]
      newname <- rownames(cc)[i]
      res <- rbind(res,c(paste("  ",newname),newrow))
      ##    res <- rbind(res,c(newname,newrow))
    }
  }
  if (Nvar>0) {
    var.idx <- which(attributes(cc)$type=="variance")
    res <- rbind(res, c("Residual Variances:",rep("",M)))
    for (i in var.idx) {
      newrow <- mycoef[i,]
      newname <- rownames(cc)[i]
      res <- rbind(res,c(paste("  ",newname),newrow))
      ##    res <- rbind(res,c(newname,newrow))
    }
  }
   res0 <- res[,-1]
  rownames(res0) <- format(res[,1],justify="left")
  res0
}

###}}} CoefMat

###{{{ standardized coefficients

stdcoef <- function(x,p=coef(x),...) {
  ##  M <- lisrel(Model(e),pars(e))
  M0 <- moments(Model(x),p)
  A <- t(M0$A)
  P <- M0$P
  v <- M0$v
  C <- M0$Cfull
  N <- diag(sqrt(diag(C)),ncol=nrow(C)); colnames(N) <- rownames(N) <- vars(x)
  iN <- N; diag(iN)[diag(N)>0] <- 1/diag(iN)[diag(N)>0]
  diag(iN)[diag(N)==0] <- NA 
  Nn <- N; Nn[] <- 0; diag(Nn) <- 1
  Nn[latent(x),latent(x)] <- N[latent(x),latent(x)]
  iNn <- Nn; diag(iNn) <- 1/diag(Nn)
  Ny <- Nn;
  Ny[endogenous(x),endogenous(x)] <- N[endogenous(x),endogenous(x)]
  iNy <- Ny; diag(iNy) <- 1/diag(Ny)
  ## Standardized w.r.t. latent,y and x:
  AstarXY <- t(iN%*%A%*%N)
  PstarXY <- iN%*%P%*%iN
  if (!is.null(v))
    vstarXY <- iN%*%v
  else
    vstarXY <- NULL
  pstdXY <- pars(Model(x),A=AstarXY,P=PstarXY,v=vstarXY)
  ## Standardized w.r.t. latent, y:
  AstarY <- t(iNy%*%A%*%Ny)
  PstarY <- iNy%*%P%*%iNy
  if (!is.null(v))
    vstarY <- iNy%*%v
  else
    vstarY <- NULL
  pstdY <- pars(Model(x),A=AstarY,P=PstarY,v=vstarY)
  ## Standardized w.r.t. latent only:
  Astar <- t(iNn%*%A%*%Nn)
  Pstar <- iNn%*%P%*%iNn
  if (!is.null(v))
    vstar <- iNn%*%v
  else
    vstar <- NULL  
  pstd <- pars(Model(x),A=Astar,Pstar,v=vstar)
  
  res <- list(par=cbind(p,pstd,pstdXY),
              AstarXY=AstarXY, PstarXY=PstarXY, vstarXY=vstarXY,
              AstarY=AstarY, PstarY=PstarY, vstarY=vstarY,
              Astar=Astar, Pstar=Pstar, vstar=vstar)
  return(res)
}

###}}} standardized coefficients
