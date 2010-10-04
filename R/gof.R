satmodel <- function(object,logLik=TRUE,data=model.frame(object),
                     control=list(start=coef(object),trace=1),
                     weight=Weight(object),estimator=object$estimator,
                     ...) {
  if (object$estimator=="gaussian" & logLik)
    return(logLik(object, type="sat"))
  covar <- exogenous(object)
  y <- endogenous(object)
  m0 <- Model(object)
  if (length(covar)>0)
    suppressWarnings(m0 <- regression(m0,y,covar))
  if (length(latent(m0))>0)
    kill(m0) <- latent(m0)
  cancel(m0) <- y
  suppressWarnings(covariance(m0) <- y)
  ## if (length(y)>1) {
  ##   for (i in 1:(length(y)-1))
  ##     for (j in (i+1):length(y)) {
  ##       m0 <- regression(m0,y[i],y[j])
  ##     }
  ## }
  cat("Calculating MLE of saturated model:\n")
  missing <- "lvm.missing"%in%class(object)  
  e0 <- estimate(m0,data=data,weight=weight,estimator=estimator,silent=TRUE,missing=missing,control=control,...)
  if (logLik)
    return(logLik(e0))
  return(e0)
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

  class(res) <- "gof.lvmfit"
  return(res)       
}

print.gof.lvmfit <- function(x,...) {
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
  invisible(x)
}
