cluster.post.hook <- function(x,...) {
  if (!is.null(x$cluster)) {
    mylev <- unique(x$cluster)
    K <- length(mylev)
    S <- score(x,indiv=TRUE) #,...)
    I <- information(x,type="hessian") #,...)
    iI <- solve(I)
    S0 <- matrix(0,ncol=ncol(S),nrow=K)
    count <- 0
    for (i in mylev) {
      count <- count+1
      S0[count,] <- colSums(S[which(x$cluster==i),])
    }
    J <- crossprod(S0)
    asVar <- iI%*%J%*%iI
  } else {
    asVar <- x$vcov
  }

  mycoef <- x$opt$estimate
  x$vcov <- asVar
  SD <- sqrt(diag(asVar))
  Z <- mycoef/SD
  pval <- 2*(1-pnorm(abs(Z)))
  newcoef <- cbind(mycoef, SD, Z, pval);

  nparall <- index(x)$npar + ifelse(x$control$meanstructure, index(x)$npar.mean,0)
  mycoef <- matrix(NA,nrow=nparall,ncol=4)
  mycoef[x$pp.idx,] <- newcoef
  colnames(mycoef) <- c("Estimate","Std. Error", "Z value", "Pr(>|z|)")    
  mynames <- c()
  if (x$control$meanstructure) {
    mynames <- vars(x)[index(x)$v1==1]
  }
  if (index(x)$npar>0) {
    mynames <- c(mynames,paste("p",1:index(x)$npar,sep=""))
  }
  rownames(mycoef) <- mynames
  x$coef <- mycoef  
  return(x)
}
