glm.estimate.hook <- function(x,estimator,...) {
  yy <- c()
  if (estimator=="glm") {
    for (y in endogenous(x)) {
      fam <- attributes(distribution(x)[[y]])$family
      if (is.null(fam)) fam <- gaussian()
      if (!(tolower(fam$family)%in%
            c("gaussian","gamma","inverse.gaussian"))) {
        yy <- c(yy,y)
      }
    }
    if (length(yy)>0) covariance(x,yy) <- 1
  }
  return(c(list(x=x,estimator=estimator,...))) 
}

GLMest <- function(m,data,control=list(),...) {
  v <- vars(m)
  yvar <- endogenous(m)  
  res <- c()
  count <- 0
  V <- NULL
  mymsg <- c()
  for (y in yvar) {
    count <- count+1
    xx <- parents(m,y)
    fam <- attributes(distribution(m)[[y]])$family
    if (is.null(fam)) fam <- gaussian()
    mymsg <- c(mymsg, with(fam, paste(family,"(",link,")",sep="")))
    if (length(xx)==0) xx <- 1
    g <- glm(toformula(y,xx),family=fam,data=data)
    p <- coef(g)
    V0 <- vcov(g)
    names(p)[1] <- y
    if (length(p)>1) {
      names(p)[-1] <- paste(y,xx,sep="<-")
    }
    colnames(V0) <- rownames(V0) <- names(p)
    if (tolower(fam$family)%in%c("gaussian","gamma","inverse.gaussian")) {
      p <- c(p,summary(g)$dispersion)
      V1 <- matrix(0) ## Not estimated!
      colnames(V1) <- rownames(V1) <- names(p)[length(p)] <- paste(y,y,sep="<->")
      V0 <- V0%+%V1
    }
    if (is.null(V)) {
      V <- V0
    } else {
      V <- V%+%V0
    }
  res <- c(res, list(p));
  }
  coefs <- unlist(res)
  idx <- match(coef(m),names(coefs))
  coefs <- coefs[idx]
  V <- V[idx,idx]
  mymsg <- noquote(cbind(mymsg))
  colnames(mymsg) <- "Family(Link)"; rownames(mymsg) <- paste(yvar,":")
  list(estimate=coefs,vcov=V,summary.message=function(...)  {
    mymsg }, dispname="Dispersion:")
}

GLMscore <- function(x,p,data,indiv=FALSE,...) {
  v <- vars(x)
  yvar <- endogenous(x)  
  S <- res <- c()
  count <- 0
  for (y in yvar) {
    count <- count+1
    xx <- parents(x,y)
    fam <- attributes(distribution(x)[[y]])$family
    if (is.null(fam)) fam <- gaussian()
    g <- glm(toformula(y,xx),family=fam,data=data)
    pdispersion <- NULL
    p0 <- p
    if (tolower(fam$family)%in%c("gaussian","gamma","inverse.gaussian")) {
      pdispersion <- tail(p,1)
      p0 <- p[-length(p)]
    }
    S0 <- rbind(score.glm(g,p=p0,indiv=indiv,...))
    if (!is.null(pdispersion)) S0 <- cbind(S0,0)    
    colnames(S0) <- names(p)
    if (is.null(S)) {
      S <- S0
    } else {
      S <- cbind(S,S0)
    }    
    res <- c(res, list(p));
  }
  coefs <- unlist(res)
  idx <- na.omit(match(coef(x),names(coefs)))
  S <- S[,idx,drop=FALSE]
  if (!indiv) S <- as.vector(S)
  return(S)
  
}


##' @S3method score glm
score.glm <- function(x,p=coef(x),data,indiv=FALSE,
                      y,X,link,dispersion,offset=NULL,...) {


  response <- all.vars(formula(x))[1]
  if (inherits(x,"glm")) {
    link <- family(x)
    a.phi <- 1
    if (tolower(family(x)$family)%in%c("gaussian","gamma","inverse.gaussian")) {
      a.phi <- summary(x)$dispersion
    }
    if (missing(data)) {
      X <- model.matrix(x)
      y <- model.frame(x)[,1]      
    } else {
      X <- model.matrix(formula(x),data=data)
      y <- model.frame(formula(x),data=data)[,1]
    }    
    offset <- x$offset
  } else {
    if (missing(link)) stop("Family needed")
    if (missing(data)) stop("data needed")
    X <- model.matrix(formula(x),data=data)
    y <- model.frame(formula(x),data=data)[,1]
  }
  n <- nrow(X)  
  g <- link$linkfun
  ginv <- link$linkinv
  dginv <- link$mu.eta ## D[linkinv]
  ##dg <- function(x) 1/dginv(g(x)) ## Dh^-1 = 1/(h'(h^-1(x)))
  canonf <- do.call(link$family,list())           
  caninvlink <- canonf$linkinv
  canlink <- canonf$linkfun
  Dcaninvlink <- canonf$mu.eta           
  Dcanlink <- function(x) 1/Dcaninvlink(canlink(x))
  ##gmu <- function(x) g(caninvlink(x))
  ##invgmu <- function(z) canlink(ginv(z))
  h <- function(z) Dcanlink(ginv(z))*dginv(z)                                
  if(any(is.na(p))) stop("Over-parameterized model")
  Xbeta <- X%*%p
  if (!is.null(offset)) Xbeta <- Xbeta+offset
  pi <- ginv(Xbeta)  
  ##res <- as.vector(y/pi*dginv(Xbeta)-(1-y)/(1-pi)*dginv(Xbeta))*X
  ##return(res)
  r <- y-pi
  A <- as.vector(h(Xbeta)*r)/a.phi 
  S <- apply(X,2,function(x) x*A)
  if (!indiv) return(colSums(S))
  return(S)
}

##' @S3method pars glm
pars.glm <- function(x,...) {
  if (tolower(family(x)$family)%in%c("gaussian","gamma","inverse.gaussian")) {
    res <- c(coef(x),summary(x)$dispersion)
    names(res)[length(res)] <- "Dispersion"
    return(res)
  }
  return(coef(x))
}

logL.glm <- function(x,p=pars.glm(x),indiv=FALSE,...) {
  f <- family(x)
  ginv <- f$linkinv
  X <- model.matrix(x)
  n <- nrow(X)  
  disp <- 1; p0 <- p
  if (tolower(family(x)$family)%in%c("gaussian","gamma","inverse.gaussian")) {
    disp <- tail(p,1)
    p0 <- p[-length(p)]
  }
  if(any(is.na(p))) stop("Over-parametrized model")
  Xbeta <- X%*%p0
  if (!is.null(x$offset)) Xbeta <- Xbeta+x$offset
  y <- model.frame(x)[,1]
  mu <- ginv(Xbeta)
  w <- x$prior.weights
  dev <-  f$dev.resids(y,mu,w)
  if (indiv) {
    
  } 
  loglik <- length(p)-(f$aic(y,n,mu,w,sum(dev))/2+x$rank)
  structure(loglik,nobs=n,df=length(p),class="logLik")
}

hessian.glm <- function(x,p=coef(x),...) {  
  jacobian(function(theta) score.glm(x,p=theta,...),p)
}

##' @S3method information glm
information.glm <- function(x,...) hessian.glm(x,...)

robustvar <- function(x,id=NULL,...) {
  U <- score(x,indiv=TRUE)
  II <- unique(id)
  K <- length(II)
  J <- 0
  if (is.null(id)) {
    J <- crossprod(U)
  } else {
    for (ii in II) {
      J <- J+tcrossprod(colSums(U[which(id==ii),,drop=FALSE]))
    }
    J <- K/(K-1)*J
  }
  iI <- vcov(x)
  V <- iI%*%J%*%iI
  return(V)  
}


glm_method.lvm <- NULL
glm_objective.lvm <- function(x,p,data,...) {
  GLMest(x,data,...)
}
glm_gradient.lvm <- function(x,p,data,...) {
  -GLMscore(x,p,data,...)
}

glm_variance.lvm <- function(x,p,data,opt,...) {
  opt$vcov
}

