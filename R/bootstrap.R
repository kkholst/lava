bootstrap <- function(x,...) UseMethod("bootstrap")

bootstrap.lvmfit <- function(x,R=100,data=model.frame(x),
                             control=list(start=coef(x)),
                             p=coef(x), parametric=FALSE,
                             estimator=x$estimator,weight=Weight(x),...)
  bootstrap.lvm(Model(x),R=R,data=data,control=control,estimator=estimator,weight=weight,parametric=parametric,p=p,...)

bootstrap.lvm <- function(x,R=100,data,fun=NULL,control=list(),
                          p, parametric=FALSE,
                          constraints=TRUE,sd=FALSE,silent=FALSE,...) {
  
  coefs <- sds <- c()
  on.exit(list(coef=coefs[-1,], sd=sds[-1,], coef0=coefs[1,], sd0=sds[1,], model=x))
  
  pmis <- missing(p)
  bootfun <- function(i) {
    if (i==0) {
      d0 <- data
    } else {
      if (!parametric | pmis) {
        d0 <- data[sample(1:nrow(data),replace=TRUE),]
      } else {
        d0 <- sim(x,p=p,n=nrow(data))
      }
    }
    e0 <- estimate(x,data=d0,control=control,silent=TRUE,index=FALSE
                   )##,...)    
    if (!silent) cat(".")
    if (!is.null(fun)) {
      coefs <- fun(e0)
      newsd <- NULL
    } else {    
      coefs <- coef(e0,symbol=c("<-","<->"))
      newsd <- c()
      if (sd) {
        newsd <- e0$coef[,2]
      }
      if (constraints & length(constrain(x))>0) {
        cc <- constraints(e0,...)
        coefs <- c(coefs,cc[,1])
        names(coefs)[seq(length(coefs)-length(cc[,1])+1,length(coefs))] <- rownames(cc)
        if (sd) {
          newsd <- c(newsd,cc[,2])
        }
      }
      ##      coefs <- rbind(coefs,newcoef)
      ##      if (sd)
      ##        sds <- rbind(sds, newsd)
    }
    return(list(coefs=coefs,sds=newsd))
  }

  i <- 0
  if (require(foreach) & lava.options()$parallel) {
    res <- foreach (i=0:R) %dopar% bootfun(i)
  } else {
    res <- lapply(0:R,bootfun)
  }
  if (!silent) cat("\n")
  
  coefs <- matrix(unlist(lapply(res, function(x) x$coefs)),nrow=R+1,byrow=TRUE)
  sds <- NULL
  if (sd)
    sds <- matrix(unlist(lapply(res, function(x) x$sds)),nrow=R+1,byrow=TRUE)

  if (!is.null(fun)) {
    rownames(coefs) <- c()
    res <- list(coef=coefs[-1,,drop=FALSE],coef0=coefs[1,],model=x) 
  } else {    
    if (constraints& length(constrain(x))>0) colnames(coefs)[-(1:length(coef(fitted)))] <- names(res[[1]]$coefs)
    colnames(coefs) <- names(res[[1]]$coefs)
    rownames(coefs) <- c(); if (sd) colnames(sds) <- colnames(coefs)
    res <- list(coef=coefs[-1,,drop=FALSE], sd=sds[-1,,drop=FALSE], coef0=coefs[1,], sd0=sds[1,], model=x)
  }
  class(res) <- "bootstrap.lvm"
  return(res)
}


"print.bootstrap.lvm" <- function(x,idx,level=0.95,...) {
  cat("Non-parametric bootstrap statistics (R=",nrow(x$coef),"):\n\n",sep="")
  uplow <-(c(0,1) + c(1,-1)*(1-level)/2)
  nn <- paste(uplow*100,"%")
  c1 <- t(apply(x$coef,2,function(x) c(mean(x), sd(x), quantile(x,uplow))))
  c1 <- cbind(c1[,1],c1[,1]-x$coef0,c1[,-1,drop=FALSE])
  colnames(c1) <- c("Estimate","Bias","Std.Err",nn)
  if (missing(idx)) {
    print(c1)
  } else {
    print(c1[idx,,drop=FALSE])
  }
  if (length(x$sd)>0) {
      c2 <- t(apply(x$sd,2,function(x) c(mean(x), sd(x), quantile(x,c(0.025,0.975)))))
      c2 <- cbind(c2[,1],c2[,1]-x$sd0,c2[,-1])
      colnames(c2) <- c("Estimate","Bias","Std.Err","2.5%","97.5%")
      cat("\nStandard errors:\n")
      if (missing(idx)) {        
        print(c2)
      } else {
        print(c2[idx,,drop=FALSE])
      }
    }
  cat("\n")
  invisible(x)
}
