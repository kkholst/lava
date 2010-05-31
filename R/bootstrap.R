bootstrap <- function(x,...) UseMethod("bootstrap")

bootstrap.lvmfit <- function(x,R=100,data=model.frame(x),control=list(start=coef(x)),...) bootstrap.lvm(Model(x),R=R,data=data,control=control,...)

bootstrap.lvm <- function(x,R=100,data,control=list(),constraints=TRUE,sd=FALSE,...) {
  m <- Model(x)
  coefs <- sds <- c()
  on.exit(list(coef=coefs[-1,], sd=sds[-1,], coef0=coefs[1,], sd0=sds[1,], model=x))
  for (i in 0:R) {
    if (i==0) {
      d0 <- data
    } else {
      d0 <- data[sample(1:nrow(data),replace=TRUE),]
    }
    e0 <- estimate(m,data=d0,control=control,silent=TRUE,index=FALSE,...)    
    cat(".")
    newcoef <- coef(e0)
    newsd <- c()
    if (sd) {
      newsd <- e0$coef[,2]
    }
    if (constraints & length(constrain(x))>0) {
      cc <- constraints(e0,...)
      newcoef <- c(newcoef,cc[,1])
      if (sd) {
        newsd <- c(newsd,cc[,2])
      }
    }
    coefs <- rbind(coefs,newcoef)
    if (sd)
      sds <- rbind(sds, newsd)
  }; cat("\n")

  if (constraints& length(constrain(x))>0) colnames(coefs)[-(1:length(coef(e0)))] <- rownames(cc)
  rownames(coefs) <- c(); if (sd) colnames(sds) <- colnames(coefs)
  res <- list(coef=coefs[-1,], sd=sds[-1,], coef0=coefs[1,], sd0=sds[1,], model=x)
  class(res) <- "bootstrap.lvm"
  return(res)
}


"print.bootstrap.lvm" <- function(x,idx,...) {
  cat("Non-parametric bootstrap statistics (R=",nrow(x$coef),"):\n\n",sep="")
  c1 <- t(apply(x$coef,2,function(x) c(mean(x), sd(x), quantile(x,c(0.025,0.975)))))
  c1 <- cbind(c1[,1],c1[,1]-x$coef0,c1[,-1])
  colnames(c1) <- c("Estimate","Bias","Std.Err","2.5%","97.5%")
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
