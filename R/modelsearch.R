modelsearch <- function(x,k=1,dir="forward") {
  if (dir!="forward")
    return(backwardsearch(x,k))
  else
    return(forwardsearch(x,k))
}


backwardsearch <- function(x,k=1) {
  if (class(x)!="lvmfit") stop("Expected an object of class 'lvmfit'.")
  p <- pars(x)
  cur <- Model(x)
  pp <- modelPar(cur,p)
  Y <- endogenous(x)
  X <- exogenous(x)
  V <- vars(x)

  p1 <- pp$p
  Tests <- c(); Vars <- list()
  AP <- with(index(cur), symmetrize(A, upper=FALSE) + P)

  parnotvar<- setdiff(1:length(p1), variances(Model(x))) ## We don't want to perform tests on the boundary of the parameter space
  freecomb <- combn(parnotvar, k)
  
  for (i in 1:ncol(freecomb))
    {
      cc0 <- coef(cur, mean=FALSE,silent=TRUE,symbol=c("->","<->"))
      ii <- freecomb[,i]
      p0 <- p1; p0[ii] <- 0
      R <- diag(length(p0)); R <- matrix(R[ii,],nrow=length(ii))
      I <- information(Model(x), p=p1, n=x$data$n)
      if (!is.null(pp$meanpar)) {
        rmidx <- 1:length(pp$meanpar)
        I <- I[-rmidx,-rmidx]
      }
      iI <- solve(I)
      W <- t(rbind(R)%*%p1)%*%solve(R%*%iI%*%t(R))%*%(cbind(R)%*%p1)
      Tests <- c(Tests, W)
      Vars <- c(Vars, list(cc0[ii]))
    }
  ord <- order(Tests, decreasing=TRUE); 
  Tests <- cbind(Tests, 1-pchisq(Tests,k)); colnames(Tests) <- c("Test Statistic", "P-value")
  res <- list(test=Tests[ord,], var=Vars[ord])
  PM <- c()
  for (i in 1:nrow(Tests)) {
    if (!is.na(res$test[i,1])) {
      newrow <- c(formatC(res$test[i,1]), formatC(res$test[i,2]), paste(res$var[[i]],collapse=", "))
      PM <- rbind(PM, newrow)
    }
  }
  colnames(PM) <- c("Wald: W", "P(W>w)", "Index"); rownames(PM) <- rep("",nrow(PM))

  res <- list(res=PM)
  class(res) <- "modelsearch"
  res
}

forwardsearch <- function(x,k=1) {
  if (class(x)!="lvmfit") stop("Expected an object of class 'lvmfit'.")
  p <- coef(x)
  cur <- Model(x)
  pp <- modelPar(cur,p)
  Y <- endogenous(x)
  X <- exogenous(x)
  V <- vars(x)
  q <- length(Y); qx <- length(X)
  npar.sat <- q+q*(q-1)/2 + q*qx
  npar.cur <- index(cur)$npar
  npar.mean <- index(cur)$npar.mean
  nfree <- npar.sat-npar.cur
  if (nfree<k) {
    cat("Cannot free",k,"variables from model.\n");
    return()
  }  
  
  Tests <- c(); Vars <- list()
  AP <- with(index(cur), symmetrize(A, upper=FALSE) + P)
    restricted <- c()
  for (i in 1:(ncol(AP)-1))
    for (j in (i+1):nrow(AP))
      if ( AP[j,i]==0 ) {
        restricted <- rbind(restricted,  c(i,j))
      }
  restrictedcomb <- combn(1:nrow(restricted), k) # Combinations of k-additions to the model

  n <- nrow(model.frame(x))
  S <- (n-1)/n*var(model.frame(x),na.rm=TRUE)
  mu <- colMeans(model.frame(x),na.rm=TRUE)

  cat(ncol(restrictedcomb), "models:\n")
  count <- 0
  for (i in 1:ncol(restrictedcomb))
    {
      cat(".")
      count <- count+1
      if (count==20) {
        cat("\n")
        count <- 0
      }

      varlist <- c()
      altmodel <- cur ## HA: altmodel, H0: cur
      for (j in 1:k) {
        myvar <- restricted[restrictedcomb[j,i],]
        covariance(altmodel) <- V[myvar]
        varlist <- rbind(varlist, V[myvar])
      }
      altmodel$parpos <- NULL
      altmodel <- updatelvm(altmodel,deriv=TRUE,zeroones=TRUE,mean=FALSE)
      ##          mu <- colMeans(model.frame(x))
      ##      altmodel$parpos <- NULL
      ##      index(altmodel)$dA <- NULL      
      ##    index(altmodel) <- reindex(altmodel) ##,deriv=TRUE,zeroones=FALSE)
##      index(altmodel)$Kkk <- index(x)$Kkk
##      index(altmodel)$Ik <- index(x)$Ik
##      index(altmodel)$Im <- index(x)$Im
##      DD <- deriv(altmodel,p=p1)
#      index(x)[names(DD)] <- DD
      
      cc <- coef(altmodel, mean=FALSE,silent=TRUE,symbol=c("->","<->"))
      cc0 <- coef(cur, mean=FALSE,silent=TRUE,symbol=c("->","<->"))
      ##          pos <- match(paste(V[i], "<->", V[j], sep=""),cc)
      ##          if (is.na(pos)) ## Should not be necessary 
      ##            pos <- match(paste(V[i], "<->", V[j], sep=""),cc)
      p1 <- numeric(length(pp$p)+k) ## New parameter basically p1 = (p0, 0)
      ## Need to be sure we place 0 at the correct position
      for (ic in 1:length(cc)) {
        idx <- match(cc[ic],cc0)
        if (!is.na(idx))
            p1[ic] <- pp$p[idx]
      }
      Sc2 <- score(altmodel,p=p1,S=S,data=NULL,n=n)
##      Sc2 <- -gaussian_gradient.lvm(altmodel,p1,S,mu=NULL,n)
      rmidx <- NULL
      if (!is.null(pp$meanpar)) {
        rmidx <- 1:length(pp$meanpar)
      } 
      ##      myidx <- (length(Sc2)-altmodel$index$npar+1):length(Sc2)
##      Sc2 <- Sc2[-rmidx]
      ##          Sc <- colSums(score(altmodel,data=model.frame(x),p=c(pp$meanpar,p1)))
      I <- information(altmodel,p1,n=x$data$n,data=NULL) ##[-rmidx,-rmidx]
      iI <- try(solve(I), silent=TRUE)
##      browser()
      Q <- ifelse (inherits(iI, "try-error"), NA, ## Score test
                   ## rbind(Sc)%*%iI%*%cbind(Sc)
                   (Sc2)%*%iI%*%t(Sc2)
                   )
      Tests <- c(Tests, Q)
##      print(Q)
      Vars <- c(Vars, list(varlist))
    }

  Tests0 <- Tests
  Vars0 <- Vars
  
  cat("\n\n")
  ord <- order(Tests);
  Tests <- cbind(Tests, 1-pchisq(Tests,k)); colnames(Tests) <- c("Test Statistic", "P-value")
  Tests <- Tests[ord,,drop=FALSE]
  Vars <- Vars[ord]
  PM <- c()
  for (i in 1:nrow(Tests)) {
    if (!is.na(Tests[i,1])) {
      vv <- apply(Vars[[i]],1,function(x) paste(x,collapse="<->"))
      newrow <- c(formatC(Tests[i,1]), formatC(Tests[i,2]), paste(vv,collapse=", "))
      PM <- rbind(PM, newrow)
    }
  }
  if (is.null(PM)) {
    cat("Saturated model\n")
    return(invisible(NULL))
  }
  colnames(PM) <- c("Score: S", "P(S>s)", "Index"); rownames(PM) <- rep("",nrow(PM))
  res <- list(res=PM, test=Tests, var=Vars)
  class(res) <- "modelsearch"
  res
}


print.modelsearch <- function(x,...) {
  print(x$res, quote=FALSE, ...)
  invisible(x)
}



