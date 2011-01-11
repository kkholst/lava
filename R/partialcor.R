
partialcor <- function(formula,data,level=0.95) {
  if (attributes(terms(formula))$response==0) {
    preds <- all.vars(formula)
    yy <- setdiff(names(data),preds)
    if (length(yy)<2)
      return(NULL)
    res <- c()
    for (i in 1:(length(yy)-1))
      for (j in (i+1):length(yy)) {
        f <- as.formula(paste("cbind(",yy[i],",",yy[j],")", paste(as.character(formula),collapse="")))
        res <- rbind(res, partialcor(f,data,level=level))
        rownames(res)[nrow(res)] <- paste(yy[i],yy[j],sep="~")
      }
    return(res)
  }
  l <- lm(formula,data)
  k <- ncol(model.matrix(l))
  n <- nrow(model.matrix(l))
  r <- residuals(l)
  rho <- cor(r)[1,2]
  zrho <- atanh(rho)
  var.z <- 1/(n-k-3)
  ci.z <- zrho + c(-1,1)*qnorm(1-(1-level)/2)*sqrt(var.z)
  ci.rho <- tanh(ci.z)
  z <- 1/sqrt(var.z)*zrho
  p.z <- 2*(pnorm(-abs(z))) # p-value using z-transform for H_0: rho=0.
  return(c(cor=rho,z=z,pval=p.z,lowerCI=ci.rho[1],upperCI=ci.rho[2]))
}
