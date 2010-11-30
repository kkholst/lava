cluster.post.hook <- function(x,...) {
  if (class(x)[1]=="multigroupfit") {
    if (is.null(x$cluster)) return(NULL)
    allclusters <- unlist(x$cluster)
    uclust <- unique(allclusters)
    K <- length(uclust)
    G <- x$model$ngroup
    S0 <- lapply(score(x,indiv=TRUE), function(x) { x[which(is.na(x))] <- 0; x })    
    S <- matrix(0,length(pars(x)),nrow=K)
    aS <- c() ##matrix(0,S[[1]]
    for (i in uclust) {
      for (j in 1:G) {
        idx <- which(x$cluster[[j]]==i)
        if (length(idx)>0)
          S[i,] <- S[i,] + colSums(S0[[j]][idx,,drop=FALSE])
      }
    }
    J <- crossprod(S)
    I <- information(x,type="hessian",...)
    iI <- solve(I)
    asVar <- iI%*%J%*%iI
    x$vcov <- asVar
    return(x)
  } 
  
  ## lvmfit:
  if (!is.null(x$cluster)) {
    uclust <- unique(x$cluster)
    K <- length(uclust)
    S <- score(x,indiv=TRUE) #,...)
    I <- information(x,type="hessian") #,...)
    iI <- solve(I)
    S0 <- matrix(0,ncol=ncol(S),nrow=K)
    count <- 0
    for (i in uclust) {
      count <- count+1
      S0[count,] <- colSums(S[which(x$cluster==i),,drop=FALSE])
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
