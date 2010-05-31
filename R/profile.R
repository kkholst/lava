profile.lvmfit <- function(fitted,idx,tau,...) {
  mm <- parfix(Model(fitted),idx,tau)
  index(mm) <- reindex(mm,zeroones=TRUE,deriv=TRUE)
  fixed <- attributes(mm)$fixed
  plogl <- function(tau0) {
    for (i in fixed$v) {
      mm$mean[[i]] <- tau0
    }
    if (length(fixed$A)>0) {
      for (i in 1:nrow(fixed$A)) {
        index(mm)$A[fixed$A[i,1],fixed$A[i,2]] <-
          mm$fix[fixed$A[i,1],fixed$A[i,2]] <- tau0
      }
    }
    if (length(fixed$P)>0) {
      for (i in 1:nrow(fixed$P))
        index(mm)$P[fixed$P[i,1],fixed$P[i,2]] <-
          mm$covfix[fixed$P[i,1],fixed$P[i,2]] <- tau0
    }
    ee <- estimate(mm,model.frame(fitted),fix=FALSE,silent=TRUE,index=FALSE,control=list(start=coef(fitted),method="NR",trace=0))
    return(logLik(ee))
  }
  val <- sapply(tau,plogl)
  attributes(val) <- NULL
  val
}

profci.lvmfit <- function(x,parm,level=0.95,...) {
  ll <- logLik(x)-qchisq(level,1)/2
  pp <- function(tau) (profile.lvmfit(x,parm,tau) - ll)
  tau0 <- coef(x)[parm]
  tau0.sd <- x$vcov[parm,parm]^0.5
  tau.range <- tau0 + 6*c(-1,1)*tau0.sd
  if (parm%in%(variances(x)+index(x)$npar.mean))
    tau.range[1] <- max(1e-6,tau.range[1])
  lower <- uniroot(pp,interval=c(tau.range[1],tau0))
  upper <- uniroot(pp,interval=c(tau0,tau.range[2]))
  return(c(lower$root,upper$root))
}
