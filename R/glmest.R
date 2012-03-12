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
      V1 <- matrix(0)
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

glm_method.lvm <- NULL
glm_objective.lvm <- function(x,p,data,...) {
  GLMest(x,data,...)
}
glm_gradient.lvm <- NULL
glm_variance.lvm <- function(x,p,data,opt,...) {
  opt$vcov
}
