"constrain<-" <- function(x,...,value) UseMethod("constrain<-")
"constrain" <- function(x,...) UseMethod("constrain")

constrain.default <- function(x,estimate=FALSE,...) {
  if (estimate) {
    return(constraints(x,...))
  }
  if (class(Model(x))[1]=="multigroup" ) {
    res <- list()
    for (m in Model(x)$lvm) {
      if (length(constrain(m))>0)
        res <- c(res, constrain(m))
    }
    return(res)
  }  
  return(Model(x)$constrain)
}

"constrain<-.multigroupfit" <-
  "constrain<-.multigroup" <- function(x,par,k=1,...,value) {
    constrain(Model(x)$lvm[[k]],par=par,...) <- value
    return(x)
}

"constrain<-.default" <- function(x,par,args,...,value) {
  if (class(par)[1]=="formula") {
    lhs <- getoutcome(par)
    xf <- attributes(terms(par))$term.labels
    par <- lhs
    args <- xf
  }
  Model(x)$constrain[[par]] <- value
  attributes(Model(x)$constrain[[par]])$args <- args
  index(Model(x)) <- reindex(Model(x))
  return(x)    
}

constraints <- function(object,vcov=object$vcov,level=0.95,data=model.frame(object),...) {
  if (!require("numDeriv")) stop("package Rgraphviz not available")
##  if (class(object)[1]=="lvm") {
##    return(constrain(object,estimate=FALSE))
##    }
  if (class(object)[1]=="multigroupfit") return(attributes(CoefMat.multigroupfit(object))$nlincon)   
  if (length(index(object)$constrain.par)<1) return(NULL)
  parpos <- Model(object)$parpos
  if (is.null(parpos)) {
    parpos <- with(index(object),matrices(Model(object),1:npar+npar.mean,meanpar=1:npar.mean))
    parpos$A[index(object)$M0==0] <- 0
    parpos$P[index(object)$P0==0] <- 0
    parpos$v[index(object)$v1==0] <- 0
  }
  myidx <- unlist(lapply(parpos$parval, function(x) {
    if (!is.null(attributes(x)$reg.idx)) {
      return(parpos$A[attributes(x)$reg.idx[1]])
    } else if (!is.null(attributes(x)$cov.idx)) {
      return(parpos$P[attributes(x)$cov.idx[1]])
    }
    else if (!is.null(attributes(x)$m.idx)) {
      return(parpos$v[attributes(x)$m.idx[1]])
    } else NA
  }))
##  myidx <- na.myidx)
  names(myidx) <- names(parpos$parval)    
  mynames <- c()
  N <- length(index(object)$constrain.par)
  res <- matrix(nrow=N,ncol=6)
  count <- 0
  for (pp in index(object)$constrain.par) {
    count <- count+1
    myc <- constrain(Model(object))[[pp]]
    mycoef <- numeric(6)
    val.idx <- myidx[attributes(myc)$args]
    val.idx0 <- na.omit(val.idx)
    mydata <- rbind(numeric(length(manifest(object))))
    colnames(mydata) <- manifest(object)
    data <- rbind(data)
    iname <- intersect(colnames(mydata),colnames(data))
    mydata[1,iname] <- unlist(data[1,iname])
    M <- modelVar(Model(object),p=pars.default(object),data=as.data.frame(mydata))
    vals <- M$parval[attributes(myc)$args]
    fval <- mycoef[1] <- myc(unlist(vals))
    myc0 <- function(theta) {
      theta0 <- unlist(vals);
##    theta0[val.idx0] <- theta[val.idx0];
      theta0[!is.na(val.idx)] <- theta
      return(myc(theta0))
    }
    vals0 <- unlist(vals)[!is.na(val.idx)]    
##  vals0 <- unlist(vals)[na.omit(val.idx)]
    if (length(vals0)==0)
      mycoef[2] <- NA
    else {
      if (!is.null(attributes(fval)$grad)) {
        Gr <- cbind(attributes(fval)$grad(unlist(vals0)))
      } else {
        Gr <- cbind(as.numeric(jacobian(myc0, unlist(vals0))))
      }
      V <- vcov[val.idx0,val.idx0]
      mycoef[2] <- (t(Gr)%*%V%*%Gr)^0.5
    }
    ## if (second) {
    ##   if (!is.null(attributes(fval)$hessian)) {
    ##     H <- attributes(fval)$hessian(unlist(vals))
    ##   } else {
    ##     H <- hessian(myc, unlist(vals))
    ##   }
    ##   HV <- H%*%vcov[val.idx,val.idx]
    ##   mycoef[1] <- mycoef[1] + 0.5*sum(diag(HV))
    ##   mycoef[2] <- mycoef[2] + 0.5*sum(diag(HV%*%HV))
    ## }
    mycoef[3] <- mycoef[1]/mycoef[2]
    mycoef[4] <- 2*(1-pnorm(abs(mycoef[3])))
    mycoef[5:6] <- mycoef[1] + c(1,-1)*qnorm((1-level)/2)*mycoef[2]
    res[count,] <- mycoef
    mynames <- c(mynames,pp)
  }
  rownames(res) <- mynames
  colnames(res) <- c("Estimate","Std. Error", "Z value", "Pr(>|z|)", "2.5%", "97.5%")
  return(res)
}
