satmodel <- function(object,logLik=TRUE,data=model.frame(object),
##                     control=list(start=coef(object),trace=1),
                     control=list(trace=1),
                     weight=Weight(object),estimator=object$estimator,
                     missing="lvm.missing"%in%class(object),
                     regr=FALSE,
                     ...) {
  if (object$estimator=="gaussian" & logLik & !missing) {
    if (class(object)[1]%in%c("multigroupfit","multigroup")) {

      ll <- structure(0,nall=0,nobs=0,df=0,class="logLik")
      for (i in seq_len(Model(object)$ngroup)) {
        l0 <- logLik(Model(Model(object))[[i]],data=model.frame(object)[[i]],type="sat")

        ll <- ll+l0
        for (atr in c("nall","nobs","df"))
          attributes(ll)[[atr]] <- attributes(ll)[[atr]]+attributes(l0)[[atr]]
      }
      
    } 
    return(logLik(object, type="sat"))
  }
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



##' Extract model summaries and GOF statistics for model object
##' 
##' Calculates various GOF statistics for model object including global
##' chi-squared test statistic and AIC. Extract model-specific mean and variance
##' structure, residuals and various predicitions.
##' 
##' 
##' @aliases gof gof.lvmfit moments moments.lvm information information.lvmfit
##' score score.lvmfit logLik.lvmfit
##' @param object Model object
##' @param x Model object
##' @param p Parameter vector used to calculate statistics
##' @param data Data.frame to use
##' @param weight Optional weight matrix
##' @param n Number of observations
##' @param conditional If TRUE the conditional moments given the covariates are
##' calculated. Otherwise the joint moments are calculated
##' @param model String defining estimator, e.g. "gaussian" (see
##' \code{estimate})
##' @param debug Debugging only
##' @param chisq Boolean indicating whether to calculate chi-squared
##' goodness-of-fit (always TRUE for estimator='gaussian')
##' @param level Level of confidence limits for RMSEA
##' @param \dots Additional arguments to be passed to the low level functions
##' @usage
##' 
##' gof(object, ...)
##'
##' \method{gof}{lvmfit}(object, chisq=FALSE, level=0.90, ...)
##' 
##' moments(x,...)
##' 
##' \method{moments}{lvm}(x, p, debug=FALSE, conditional=FALSE, data=NULL, ...)
##' 
##' \method{logLik}{lvmfit}(object, p=coef(object),
##'                       data=model.frame(object),
##'                       model=object$estimator,
##'                       weight=Weight(object),
##'                           ...)
##' 
##' \method{score}{lvmfit}(x, data=model.frame(x), p=pars(x), model=x$estimator, weight=Weight(x), ...)
##' 
##' \method{information}{lvmfit}(x,p=pars(x),n=x$data$n,data=model.frame(x),model=x$estimator,weight=Weight(x),...)
##' 
##' @return A \code{htest}-object.
##' @author Klaus K. Holst
##' @keywords methods models
##' @export
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

##' @S3method gof lvmfit
gof.lvmfit <- function(object,chisq=FALSE,level=0.90,...) {
  n <- object$data$n
  loglik <- logLik(object,...)
  
  df <- attributes(loglik)$df
  nobs <- attributes(loglik)$nall*length(endogenous(object))
  myAIC <- -2*(loglik - df); attributes(myAIC) <- NULL
  myBIC <- -2*loglik + df*log(nobs); attributes(myBIC) <- NULL

  xconstrain <- intersect(unlist(lapply(constrain(object),function(z) attributes(z)$args)),manifest(object))
  
  if (class(object)[1]=="lvmfit" & (object$estimator=="gaussian" | chisq) & length(xconstrain)==0 ) {
    res <- list(fit=compare(object), n=n, logLik=loglik, BIC=myBIC, AIC=myAIC)
    q <- res$fit$statistic
    qdf <- res$fit$parameter
    epsilon <- function(lambda) sapply(lambda,function(x)
                                       ifelse(x>0 & qdf>0,sqrt(x/(qdf*(n-1))),0))
                                       ##sqrt(max(0,x/(qdf*(n-1)))))
    opf <- function(l,p) (p-pchisq(q,df=qdf,ncp=l))^2
    alpha <- (1-level)/2
    hi <- list(par=0)
    RMSEA <- epsilon(q-qdf)
    start <- RMSEA
    if (RMSEA>0) {
      hi <- optimize(function(x) opf(x,p=1-alpha),c(0,q-qdf)); hi$par <- hi$minimum
      hi <- tryCatch(nlminb(hi$par^0.5,function(x) opf(x^2,p=1-alpha)),error=function(...) list(par=NA)); hi$par <- hi$par^2
    }
    lo <- optimize(function(x) opf(x,p=alpha),c(q-qdf,n)); lo$par <- lo$minimum
    lo <- tryCatch(nlminb(lo$par^0.5,function(x) opf(x^2,p=alpha)),error=function(...) list(par=NA)); lo$par <- lo$par^2
    ci <- c(epsilon(c(hi$par,lo$par)))    
    RMSEA <- c(RMSEA=RMSEA,ci);
    names(RMSEA) <- c("RMSEA",paste(100*c(alpha,(1-alpha)),"%",sep=""))
    res <- c(res,list(RMSEA=RMSEA, level=level))
  } else {
    res <- list(n=n, logLik=loglik, BIC=myBIC, AIC=myAIC)
  }

  l2D <- sum(object$opt$grad^2)
  rnkV <- tryCatch(qr(vcov(object))$rank,error=function(...) NULL)
  condnum <- tryCatch(condition(vcov(object)),error=function(...) NULL)
  res <- c(res, L2score=l2D, rankV=rnkV, cond=condnum, k=nrow(vcov(object)))
  class(res) <- "gof.lvmfit"
  return(res)       
}

##' @S3method print gof.lvmfit
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
           "Chi-squared statistic: q =", fit$statistic, 
           ", df =", fit$parameter, 
           ", P(Q>q) =", fit$p.value, "\n"))
  if (optim) {
    if (!is.null(x$RMSEA)) {
      rr <- round(x$RMSEA*10000)/10000
      rmsea <- paste(rr[1]," (",rr[2],";",rr[3],")",sep="")
      cat(" RMSEA (",x$level*100,"% CI): ", rmsea,"\n",sep="")
    }

    cat("rank(Information) = ",x$rankV," (p=", x$k,")\n",sep="")
    cat("condition(Information) = ",x$cond,"\n",sep="")
    cat("||score||^2 =",x$L2score,"\n")
  }

  invisible(x)
}


