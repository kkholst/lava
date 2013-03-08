##' Dose response calculation for binomial regression models
##'
##' @title Dose response calculation for binomial regression models
##' @param model Model object
##' @param intercept Index of intercept parameters
##' @param slope Index of intercept parameters
##' @param prob Index of mixture parameters (only relevant for
##' \code{curereg} models)
##' @param x Optional weights
##' length(x)=length(intercept)+length(slope)+length(prob)
##' @param level Probability at which level to calculate dose
##' @param ci.level Level of confidence limits
##' @param EB Optional ratio of treatment effect and adverse effects
##' used to find optimal dose (regret-function argument)
##' @author Klaus K. Holst
##' @export
PD <- function(model,intercept=1,slope=2,prob=NULL,x,level=0.5,ci.level=0.95,
               EB=NULL) { 
  if (length(intercept)<length(coef(model))) {
    B.intercept <- rep(0,length(coef(model)));    
    B.intercept[intercept] <- 1
  } else {
    B.intercept <- intercept
  }
  if (length(slope)<length(coef(model))) {
    B.slope <- rep(0,length(coef(model)));    
    B.slope[slope] <- 1
  } else {
    B.slope <- slope
  }  
  if (!is.null(prob)) {
    if (length(prob)<length(coef(model))) {
      B.prob <- rep(0,length(coef(model)));    
      B.prob[prob] <- 1
    } else {
      B.prob <- prob
    }
  }
  if (is.null(prob)) B.prob <- NULL
  B <- rbind(B.intercept,B.slope,B.prob)
  S <- B%*%vcov(model)%*%t(B)
  b <- as.vector(B%*%coef(model))

  f <- function(b) {
    mylevel <- level
    if (!is.null(EB)) {      
      if (is.null(prob)) stop("Index of mixture-probability parameters needed")
      pi0 <- family(model)$linkinv(b[3])
      mylevel <- 1-(1-pi0)/pi0*(EB)/(1-EB)
    }
    return(structure((family(model)$linkfun(mylevel)-b[1])/b[2],level=mylevel))
  }

  xx <- f(b)
  Dxx <- -1/b[2]*rbind(1,xx)
  if (!is.null(EB))
    Dxx <- grad(f,b)
  se <- diag(t(Dxx)%*%S%*%Dxx)^0.5
  res <- cbind(Estimate=xx,"Std.Err"=se)
  alpha <- 1-ci.level
  alpha.str <- paste(c(alpha/2,1-alpha/2)*100,"",sep="%")
  res <- cbind(res,res[,1]-qnorm(1-alpha/2)*res[,2],res[,1]+qnorm(1-alpha/2)*res[,2])
  colnames(res)[3:4] <- alpha.str  
  rownames(res) <- paste(round(1000*attributes(xx)$level)/10,"%",sep="")
  structure(res,b=b)
}

##' Regression model for binomial data with unkown group of immortals
##' 
##' @title Regression model for binomial data with unkown group of immortals
##' @param formula Formula specifying 
##' @param cureformula Formula for model of disease prevalence
##' @param data data frame
##' @param fam Distribution family (see the help page \code{family})
##' @param offset Optional offset
##' @param start Optional starting values
##' @param var Type of variance (robust, expected, hessian, outer)
##' @param ... Additional arguments to lower level functions
##' @author Klaus K. Holst
##' @export
curereg <- function(formula,cureformula=~1,data,fam=binomial(),offset=NULL,start,var="robust",...) {
  md <- model.frame(formula,data)
  y <- md[,1]
  X <- model.matrix(formula,data)
  Z <- model.matrix(cureformula,data)
  beta.idx <- seq(ncol(X)); gamma.idx <- seq(ncol(Z))+ncol(X)
  if (missing(start)) start <- rep(0,ncol(X)+ncol(Z))
  op <- nlminb(start,function(x)
               -curereg_logL(x[beta.idx],x[gamma.idx],y,X,Z),
               grad=function(x)
               -curereg_score(x[beta.idx],x[gamma.idx],y,X,Z))
  beta <- op$par[beta.idx]; gamma <- op$par[gamma.idx]
  cc <- c(beta,gamma)
  names(cc) <- c(colnames(X),paste("pr:",colnames(Z),sep=""))
  I <- curereg_information(beta,gamma,y,X,Z,offset,type=var,...)
  V <- lava:::Inverse(I); colnames(V) <- rownames(V) <- names(cc)
  res <- list(coef=cc,opt=op,beta=beta,gamma=gamma,
              beta.idx=beta.idx,gamma.idx=gamma.idx,
              I=I,formula=formula,cureformula=cureformula, y=y, X=X, Z=Z, offset=offset, fam=fam, vcov=V, model.frame=md)
  class(res) <- "curereg"
  res$fitted.values <- predict(res)
  return(res)
}

##' @S3method vcov curereg
vcov.curereg <- function(object,...) object$vcov
##' @S3method coef curereg
coef.curereg <- function(object,...) object$coef
##' @S3method family curereg
family.curereg <- function(object,...) object$fam
##' @S3method predict curereg
predict.curereg <- function(object,beta=object$beta,gamma=object$gamma,newdata,link=TRUE,subdist=FALSE,...) {
  newf <- as.formula(paste("~",as.character(object$formula)[3]))
  if (missing(newdata)) {
    X <- object$X; Z <- object$Z
  } else {
    X <- model.matrix(newf,newdata)
    Z <- model.matrix(object$cureformula,newdata)
  }
  if (length(beta)==length(object$beta)+length(object$gamma)) {
    gamma <- beta[object$gamma.idx]
    beta <- beta[object$beta.idx]
  }    
  g <- object$fam$linkfun
  ginv <- object$fam$linkinv
  dginv <- object$fam$mu.eta ## D[linkinv]  
  Xbeta <- as.vector(X%*%beta)
  Zgamma <- as.vector(Z%*%gamma)
  Pred <- ginv(Xbeta)
  if (subdist) return(Pred)
  p0 <- ginv(Zgamma)
  Pred <- p0*Pred
  A1 <- p0*dginv(Xbeta)
  A2 <- ginv(Xbeta)*dginv(Zgamma)  
  dgamma <- apply(Z,2,function(z) A2*z)
  dbeta <- apply(X,2,function(x) A1*x)  
  attributes(Pred)$grad <- cbind(dbeta,dgamma)
  return(Pred)  
}

##' @S3method residuals curereg
residuals.curereg <- function(object,newdata,...) {
  if (missing(newdata)) {
    y <- object$y
  } else {
    y <- model.frame(object$formula,newdata)[,1]
  }
  y-predict(object,newdata=newdata,...)
}

##' @S3method summary curereg
summary.curereg <- function(object,level=0.95,...) {
  alpha <- 1-level
  alpha.str <- paste(c(alpha/2,1-alpha/2)*100,"",sep="%")
  cc <- cbind(coef(object),diag(vcov(object))^0.5)
  pval <- 2*(1-pnorm(abs(cc[,1]/cc[,2])))
  cc <- cbind(cc[,1],cc[,1]-qnorm(1-alpha/2)*cc[,2],cc[,1]+qnorm(1-alpha/2)*cc[,2],pval)
  colnames(cc) <- c("Estimate",alpha.str,"P-value")

  
  return(structure(list(coef=cc),class="summary.curereg"))
}

##' @S3method print summary.curereg
print.summary.curereg <- function(x,...) {
  printCoefmat(x$coef,...)
}

##' @S3method print curereg
print.curereg <- function(x,...) {
 print(summary(x,...))
}

##' @S3method logLik curereg
logLik.curereg <- function(object,beta=object$beta,gamma=object$gamma,data,offset=object$offset,indiv=FALSE,...) {
  if (!missing(data)) {
    y <- model.frame(object$formula,data)[,1]
    X <- model.matrix(object$formula,data)
    Z <- model.matrix(object$cureformula,data)
    return(curereg_logL(beta,gamma,y,X,Z,offset,object$fam,indiv=indiv,...))
  }    
  curereg_logL(beta,gamma,object$y,object$X,object$Z,offset,object$fam,indiv=indiv,...)
}
curereg_logL <- function(beta,gamma,y,X,Z,offset=NULL,fam=binomial(),indiv=FALSE,...) {
  g <- fam$linkfun
  ginv <- fam$linkinv
  dginv <- fam$mu.eta ## D[linkinv]
  n <- nrow(X)  
  Xbeta <- as.vector(X%*%beta)
  Zgamma <- as.vector(Z%*%gamma)
  p0 <- ginv(Zgamma)
  if (!is.null(offset)) Xbeta <- Xbeta+offset
  Pr <- p0*ginv(Xbeta)
  loglik <- y*log(Pr)+(1-y)*log(1-Pr)
  if (indiv) return(loglik)
  loglik <- sum(loglik)
  structure(loglik,nobs=n,df=length(beta)+1,class="logLik")
}

##' @S3method score curereg
score.curereg <- function(x,beta=x$beta,gamma=x$gamma,data,offset=x$offset,indiv=FALSE,...) {
  if (!missing(data)) {
    y <- model.frame(object$formula,data)[,1]
    X <- model.matrix(object$formula,data)
    Z <- model.matrix(object$cureformula,data)
    s <- curereg_score(beta,gamma,y,X,Z,offset,x$fam,indiv=indiv,...)
  } else {    
    s <- curereg_score(beta,gamma,x$y,x$X,x$Z,offset,x$fam,indiv=indiv,...)
  }
  if (indiv) colnames(s) <- names(x$coef) else names(s) <- names(x$coef)
  return(s)
}

curereg_score <- function(beta,gamma,y,X,Z,offset=NULL,fam=binomial(),indiv=FALSE,...) {
  g <- fam$linkfun
  ginv <- fam$linkinv
  dginv <- fam$mu.eta ## D[linkinv]
  n <- nrow(X)  
  Xbeta <- as.vector(X%*%beta)
  Zgamma <- as.vector(Z%*%gamma)
  p0 <- ginv(Zgamma)
  if (!is.null(offset)) Xbeta <- Xbeta+offset
  Pr <- p0*ginv(Xbeta)
  A0 <- (y/Pr  - (1-y)/(1-Pr))
  A1 <- A0*p0*dginv(Xbeta)
  A2 <- A0*ginv(Xbeta)*dginv(Zgamma)  
  dbeta <- apply(X,2,function(x) A1*x)
  dgamma <- apply(Z,2,function(z) A2*z)
  ss <- cbind(dbeta,dgamma)
  if (indiv) return(ss)
  colSums(ss)
}

##' @S3method information curereg
information.curereg <- function(x,beta=x$beta,gamma=x$gamma,data,offset=x$offset,type=c("robust","outer","obs"),...) {
  if (!missing(data)) {
    y <- model.frame(object$formula,data)[,1]
    X <- model.matrix(object$formula,data)
    Z <- model.matrix(object$cureformula,data)
    I <- curereg_information(beta,gamma,y,X,Z,offset,x$fam,type=type,...)
  } else {
    I <- curereg_information(beta,gamma,x$y,x$X,x$Z,offset,x$fam,type=type,...)
  }
  colnames(I) <- rownames(I) <- names(x$coef)
  return(I)
}

curereg_information <- function(beta,gamma,y,X,Z,offset=NULL,fam=binomial(),type=c("outer","obs","robust"),...) {
  if (tolower(type[1])%in%c("obs")) {
    beta.idx <- seq(ncol(X)); gamma.idx <- seq(ncol(Z))+ncol(X)
    I <- -jacobian(function(x)
                   curereg_score(x[beta.idx],x[gamma.idx],y,X,Z,offset,fam,...),c(beta,gamma))
    return(I)
  }
  if (tolower(type[1])%in%c("robust","sandwich")) {
    I <- curereg_information(beta,gamma,y,X,Z,offset,fam,type="obs")
    J <- curereg_information(beta,gamma,y,X,Z,offset,fam,type="outer")
    return(J%*%solve(I)%*%J)
  }
  S <- curereg_score(beta,gamma,y,X,Z,offset,fam,indiv=TRUE,...)
  crossprod(S)
}

