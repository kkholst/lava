
effects.lvmfit <- function(object,to,from,silent=FALSE,...) {
  if (missing(to)) {
    return(summary(object))
  }

  P <- path(object,to=to,from=from,...)
  from <- P$path[[1]][1]
  to <- tail(P$path[[1]],1)

  cc <- coef(object,level=9) ## All parameters (fixed and variable)
  cc0 <- cbind(coef(object)) ## Estimated parameters
  i1 <- na.omit(match(rownames(cc),rownames(cc0)))
  idx.cc0 <-  which(rownames(cc)%in%rownames(cc0)); ## Position of estimated parameters among all parameters
  S <- matrix(0,nrow(cc),nrow(cc)); rownames(S) <- colnames(S) <- rownames(cc)
  V <- object$vcov
  npar.mean <- index(object)$npar.mean
##  if (object$control$meanstructure & npar.mean>0)
##    V <- V[-c(1:npar.mean),-c(1:npar.mean)]
  S[idx.cc0,idx.cc0] <- V[i1,i1] ## "Covariance matrix" of all parameters
  
  idx.orig <- unique(unlist(P$idx))
  coefs.all <- cc[idx.orig]
  S.all <- S[idx.orig,idx.orig]
  idx.all <- numberdup(unlist(P$idx))
  pos <- 1; idx.list <- P$idx; for (i in 1:length(idx.list)) {
    K <- length(idx.list[[i]])
    idx.list[[i]] <- idx.all[pos:(pos+K-1)]; pos <- pos+K
  }
  totalef <- prodsumdelta(coefs.all, idx.list, S.all,...)
  margef <- list(); for (i in 1:length(idx.list)) {
    margef <- c(margef, list(prodsumdelta(coefs.all, idx.list[i], S.all,...)))
  }
  
  directidx <- which(lapply(P$path,length)==2)
  if (length(directidx)==0)
    directef <- list(est=0, sd=NA)
  else
    directef <- margef[[directidx]]
  val <- list(paths=P$path, totalef=totalef, directef=directef, margef=margef, from=from, to=to)
  class(val) <- "effects"
  ##    res <- c(res, list(val))
  val
}

print.effects <- function(x,...) {
  with(x, {
    cat("\nTotal effect of '", from, "' on '", to, "':\n", sep="")
    cat("\t\t", totalef$est, " (Approx. Std.Err = ", totalef$sd, ")\n", sep="")
    cat("Direct effect of '", from, "' on '", to, "':\n", sep="")
    cat("\t\t", directef$est, " (Approx. Std.Err = ", directef$sd, ")\n", sep="")
  
    cat("Indirect effects:\n");
    for (i in 1:length(margef)) {
      if (length(paths[[i]])>2) {
        cat("\tEffect of '", from, "' via ", paste(paths[[i]],collapse="->"), ":\n", sep="");
        cat("\t\t", margef[[i]]$est, " (Approx. Std.Err = ", margef[[i]]$sd, ")\n", sep="")
      }
    }
  }) 
  cat("\n");
  invisible(x) 
}


coef.effects <- function(object,...) {  
  totalef <- with(object$totalef, cbind(est,sd[1]))
  directef <- with(object$directef, cbind(est,sd[1]))
  rownames(totalef) <- "Total"
  rownames(directef) <- "Direct"
  nn <- indirectef <- c()
  K <- seq_len(length(object$margef))
  for (i in K) {
    if (length(object$paths[[i]])>2) {        
      nn <- c(nn,paste(rev(object$paths[[i]]),collapse="<-"))
      indirectef <- rbind(indirectef, with(object$margef[[i]], c(est,sd)))
      }
  }; rownames(indirectef) <- nn  
  mycoef <- rbind(totalef,directef,indirectef)
  mycoef <- cbind(mycoef,mycoef[,1]/mycoef[,2])
  mycoef <- cbind(mycoef,2*(1-pnorm(abs(mycoef[,3]))))
  colnames(mycoef) <- c("Estimate","Std.Err","z value","Pr(>|z|)")
  mycoef
}

confint.effects <- function(object,level=0.95,...) {
  mycoef <- coef(object)
  p <- 1-(1-level)/2
  res <- mycoef[,1] +  + qnorm(p)*cbind(-1,1)%x%mycoef[,2]
  colnames(res) <- paste(c(1-p,p)*100,"%",sep="")
  rownames(res) <- rownames(mycoef)
  res  
}


prodtrans <- function(betas) {
  k <- length(betas)
  res <- prod(betas)
  ##  if (all(betas>0)) {
  ##    attr(res,"gradient") <- res/betas
  ##    return(res)
  ##  }
  nabla <- numeric(k)
  for (i in 1:k)
    nabla[i] <- prod(betas[-i])    
  
  H <- matrix(0,k,k)
  if (k>1)
    for (i in 1:(k-1))
      for (j in (i+1):k)
        H[j,i] <- H[i,j] <- prod(c(1,betas[-c(i,j)]))    
  attr(res,"gradient") <- nabla
  attr(res,"hessian") <- H  
  return(res)
}
prodsumdelta <- function(betas,prodidx,S,order=1) { ## Delta-method
  k <- length(prodidx)
  p <- length(betas)
  if (p==1) {
    return(list(est=betas, sd=sqrt(S), grad=0, hess=0))
  }
  val <- 0; grad <- numeric(p)
  H <- matrix(0,p,p)
  for (i in 1:k) {
    ii <- prodidx[[i]]
    myterm <- prodtrans(betas[ii]);
    if (order>1) {
      H0 <- attributes(myterm)$hessian
      Sigma <- S[ii,ii]
       ## print(Sigma)
       ## print(H0)
       ## print(Sigma%*%H0)
      print(sum(diag(Sigma%*%H0))/2)
      val <- val + (myterm + sum(diag(Sigma%*%H0))/2)
    } else {
      val <- val + myterm      
    }
    grad[ii] <- grad[ii] + attributes(myterm)$gradient    
  }; grad <- matrix(grad,ncol=1)
  return(list(est=val, sd=sqrt(t(grad)%*%S%*%grad), grad=grad, hess=H))
}

