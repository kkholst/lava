Isqrt <- function(X) {
    eX <- eigen(X);
    with(eX, vectors %*% diag(1/sqrt(values),nrow=length(values)) %*% t(vectors))
}


##' @export
residuals.multigroupfit <- function(object,data=model.frame(object),p=coef(object), k, ...) {
  pp <- modelPar(object,p,...)
  if (!missing(k)) return(residuals(object$model$lvm[[k]],data=data[[k]],p=pp$p[[k]],...))
  res <- c()
  for (i in seq(length(pp$p))) {
    res <- c(res, list(residuals(object$model$lvm[[i]],data=data[[i]],p=pp$p[[i]],...)))
  }
  return(res)
}


##' @export
residuals.lvmfit <- function(object,data=model.frame(object),p=coef(object),...) {
  residuals(Model(object), data=data, p=p, ...)
}

##' @export
residuals.lvm <- function(object,data=model.frame(object),std=FALSE,p=coef(object),...) {
  Y <- setdiff(manifest(object), X <- exogenous(object))
  Pr <- predict(object,p=p,data=data)
  PrY <- Pr[,Y,drop=FALSE]
  ##  y <- endogenous(object)[match(endogenous(object),manifest(object))]
  r <- as.matrix(data[,Y,drop=FALSE]-(PrY))
  res <- r

  if (std) {
    S <- attributes(Pr)$cond.var;
    if (length(Y)>1) {
      res <- r%*%Isqrt(S)
    } else res <- 1/sqrt(S[1,1])*r
  }
  colnames(res) <- colnames(r)
  res
}

