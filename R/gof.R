satmodel <- function(object,logLik=TRUE,data=model.frame(object),
##                     control=list(start=coef(object),trace=1),
                     control=list(trace=1),
                     weight=Weight(object),estimator=object$estimator,
                     missing="lvm.missing"%in%class(object),
                     regr=FALSE,
                     ...) {
  if (object$estimator=="gaussian" & logLik & !missing)
    return(logLik(object, type="sat"))
  covar <- exogenous(object)
  y <- endogenous(object)
  m0 <- Model(object)
  if (length(covar)>0)
    suppressWarnings(m0 <- regression(m0,y,covar))
  if (length(latent(m0))>0)
    kill(m0) <- latent(m0)
  cancel(m0) <- y
  if (!regr)
    suppressWarnings(covariance(m0) <- y)
  else {
    if (length(y)>1) {
      for (i in 1:(length(y)-1))
       for (j in (i+1):length(y)) {
         m0 <- regression(m0,y[i],y[j])
       }
    }
    exogenous(m0) <- covar
  }
  if (is.null(control$start)) {
    mystart <- rep(0,with(index(m0), npar.mean+npar))
    mystart[variances(m0,mean=TRUE)] <- 1
    control$start <- mystart
  }  
  cat("Calculating MLE of saturated model:\n")
  e0 <- estimate(m0,data=data,weight=weight,estimator=estimator,silent=TRUE,control=control,missing=missing,...)
  if (logLik)
    return(logLik(e0))
  return(e0)
}

condition <- function(A) {
  suppressWarnings(with(eigen(A),tail(values,1)/head(values,1)))
}

`gof` <-
  function(object,...) UseMethod("gof")

## gof.multigroupfit <- function(object,...) {
##   L0 <- logLik(object); df0 <- attributes(L0)$df
##   L1 <- logLik(object,type="sat"); df1 <- attributes(L1)$df

##   df <- df1-df0; names(df) <- "df"
##   Q <- -2*(L0-L1); attributes(Q) <- NULL; names(Q) <- "chisq";
##   pQ <- 1-pchisq(Q,df)
##   values <- c(L0,L1); names(values) <- c("log likelihood (model)", "log likelihood (saturated model)")
##   res <- list(statistic = Q, parameter = df,
##               p.value=pQ, method = "Likelihood ratio test",
##               estimate = values)
##   class(res) <- "htest"
##   return(res)    
## }


gof.lvmfit <- function(object,chisq=FALSE,...) {
  n <- object$data$n
  loglik <- logLik(object,...)
  
  df <- attributes(loglik)$df
  nobs <- attributes(loglik)$nall*length(endogenous(object))
  myAIC <- -2*(loglik - df); attributes(myAIC) <- NULL
  myBIC <- -2*loglik + df*log(nobs); attributes(myBIC) <- NULL
  
  if (class(object)[1]=="lvmfit" & (object$estimator=="gaussian" | chisq)   )
    res <- list(fit=compare(object), n=n, logLik=loglik, BIC=myBIC, AIC=myAIC, model=object)
  else
    res <- list(n=n, logLik=loglik, BIC=myBIC, AIC=myAIC, model=object)

  l2D <- sum(object$opt$grad^2)
  rnkV <- qr(vcov(object))$rank
  res <- c(res, L2score=l2D, rankV=rnkV, cond=condition(vcov(object)), k=nrow(vcov(object)))
  class(res) <- "gof.lvmfit"
  return(res)       
}

print.gof.lvmfit <- function(x,optim=TRUE,...) {
  if (!is.null(x$n))
    with(x,       
         cat("Number of observations =", n, "\n"))
  with(x,
       cat(" Log-Likelihood =", logLik, "\n",
           "BIC =", BIC, "\n",
           "AIC =", AIC, "\n"))
  if (!is.null(x$fit))
  with(x,
       cat(" log-Likelihood of model =", fit$estimate[1], "\n",
           "log-Likelihood of saturated model =", fit$estimate[2], "\n",
           "Chi-squared statistic: Q =", fit$statistic, 
           ", df =", fit$parameter, 
           ", P(Q>q) =", fit$p.value, "\n"))
  if (optim) {
    cat("rank(Information) = ",x$rankV," (p=", x$k,")\n",sep="")
    cat("condition(Information) = ",x$cond,"\n",sep="")
    cat("||score||^2 =",x$L2score,"\n")
  }

  invisible(x)
}


