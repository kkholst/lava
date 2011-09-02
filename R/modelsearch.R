equivalence <- function(x,rel,tol=1e-3,k=1,omitrel=TRUE,...) {
  if (class(rel)[1]=="formula") {
    myvars <- all.vars(rel)
  } else {    
    myvars <- rel
  }
  if (length(myvars)!=2) stop("Two variables only")
  x0 <- Model(x)
  cancel(x0) <- rel  
  e0 <- estimate(x0,data=model.frame(x),weight=Weight(x),estimator=x$estimator,...)
  if (k!=1) {
    p0 <- coef(x)
    p0[] <- 0
    p0[match(names(coef(e0)),names(p0))] <- coef(e0)
    S0 <- score(x,p=p0)[,,drop=TRUE]; 
    I0 <- information(x,p=p0)   
    T0 <- rbind(S0)%*%solve(I0)%*%cbind(S0); names(T0) <- "Q"
  } 
  s <- modelsearch(e0,k=k)
  relname <- c(paste(myvars,collapse="<->"),
               paste(rev(myvars),collapse="<->"))
  relidx <- NULL
  if (k==1) {
    relidx <- na.omit(match(relname,s$res[,"Index"]))
    T0 <- s$test[relidx,1]
  }
##  else {
##    if (covariance(Model(x))$rel[myvars[1],myvars[2]]==0) {    
    ## paridx <- which(names(coef(x))%in%relname)
    ## if (length(paridx)==0) {
    ##   paridx <- which(names(coef(x))%in%paste(myvars,collapse="<-"))
    ## }    
##    p0 <- coef(x)
##    p0[paridx] <- 0
  ##   p0[] <- 0
  ##   p0[match(names(coef(e0)),names(p0))] <- coef(e0)
  ##   S0 <- score(x,p=p0)[,,drop=TRUE]
  ##   I0 <- information(x,p=p0,data=NULL,n=nrow(model.frame(x)))   
  ##   T0 <- rbind(S0)%*%solve(I0)%*%cbind(S0)    
  ## }
  T <- s$test[,1]
  Equiv <- setdiff(which(abs(T-T0)<tol),relidx)
  Improve <- which((T-T0)>tol)
  if (omitrel) { ## Don't save models including 'rel'
    keep <- c()
    if (length(Equiv)>0) {
      for (i in 1:length(Equiv)) {
        newvars <- s$var[[Equiv[i]]]
        if (!any(apply(newvars,1,function(z) all(z%in%myvars)))) keep <- c(keep,Equiv[i])
      }
      Equiv <- keep
    }
    keep <- c()
    if (length(Improve)>0) {
      for (i in 1:length(Improve)) {
        newvars <- s$var[[Improve[i]]]
        if (!any(apply(newvars,1,function(z) all(z%in%myvars)))) keep <- c(keep,Improve[i])
      }
      Improve <- keep
    }
  }
  eqvar <- ivar <- NULL
  models <- list()
  if (length(Equiv)>0){
    for (i in 1:length(Equiv)) {      
      xnew <- x0
      newvars <- s$var[[Equiv[i]]]
      for (j in 1:nrow(newvars)) {
        exo.idx <- which(newvars[j,]%in%index(x0)$exogenous)
        if (length(exo.idx)>0) {
          xnew <- regression(xnew,from=newvars[j,exo.idx],to=newvars[j,setdiff(1:2,exo.idx)])
        } else {
          covariance(xnew) <- newvars
        }
      }
      models <- c(models,list(xnew))
    }
    eqvar <- s$var[Equiv]
  }
  if (length(Improve)>0)   {
      for (i in 1:length(Improve)) {
      xnew <- x0
      newvars <- s$var[[Improve[i]]]
      for (j in 1:nrow(newvars)) {
        exo.idx <- which(newvars[j,]%in%index(x0)$exogenous)
        if (length(exo.idx)>0) {
          xnew <- regression(xnew,from=newvars[j,exo.idx],to=newvars[j,setdiff(1:2,exo.idx)])
        } else {
          covariance(xnew) <- newvars
        }
      }
      models <- c(models,list(xnew))
    }
    ivar <- s$var[Improve]
  }  
  res <- list(equiv=eqvar, improve=ivar, scoretest=s, models=models, I=Improve, E=Equiv, T0=T0, vars=myvars)
  class(res) <- "equivalence"
  return(res)
}

print.equivalence <- function(x,...) {
  cat("  0)\t ",paste(x$vars,collapse="<->"),"  (",formatC(x$T0),")\n",sep="")
  cat("Empirical equivalent models:\n")
  if (length(x$E)==0)
    cat("\t none\n")
  else
    for (i in 1:length(x$E)) {        
      cat("  ",i,")\t ",  x$scoretest$res[x$E[i],"Index"],
          "  (",x$scoretest$res[x$E[i],1],")",
          "\n",sep="")
    }
  cat("Candidates for model improvement:\n")
  if (length(x$I)==0)
    cat("\t none\n")
  else
  for (i in 1:length(x$I)) {
      cat("  ",i,")\t ",  x$scoretest$res[x$I[i],"Index"],
          "  (",x$scoretest$res[x$I[i],1],")",
          "\n",sep="")
  }
  invisible(x)
}

holm <- function(p) {
  k <- length(p)
  w <- 1/k
  ii <- order(p)
  po <- p[ii]
  qs <- min(1,po[1]/w)
  for (i in 2:k) {
      qs <- c(qs, min(1, max(qs[i-1],po[i]*(1-w*(i-1))/w)))
    }
  return(qs)
}

modelsearch <- function(x,k=1,dir="forward",...) {
  if (dir=="forward") {
    res <- forwardsearch(x,k,...)
    return(res)
  }
  if (dir=="backstep") {
    res <- backwardeliminate(x,...)
    return(res)
  }
  res <- backwardsearch(x,k,...)
  return(res)
}

backwardeliminate <- function(x,
                              keep=NULL,
                              pthres=0.05,AIC=FALSE,silent=TRUE,
                              missing=FALSE,intercepts=FALSE,
                              maxsteps=Inf,
                              information="E",
                              messages=TRUE,
                              data,
                              ...) {

  if (class(x)[1]=="lvm") { M <- x } else { M <- Model(x) }
  if(missing(data)) data <- model.frame(x)

  if (intercepts) ii <- NULL
  ff <- function() {
    ii <- grep("m",names(coef(M)))
    vv <- variances(M,mean=TRUE)
    cc <- estimate(M,data,quick=TRUE,silent=silent,missing=missing,...)
    I0 <- information(M,p=cc,data=data,type=information)[-c(ii,vv),-c(ii,vv)]
    cc0 <- cc[-c(ii,vv)]
    res <- (1-pnorm(abs(cc0/sqrt(diag(solve(I0))))))*2
    return(res)            
  }
  
  done <- FALSE; i <- 0;
  while (!done & i<maxsteps) {    
    p <- ff(); ordp <- order(p,decreasing=TRUE)
    curp <- p[ordp[1]]
    if (curp<pthres) break;
    var1 <- unlist(strsplit(names(curp),"<-"))
    if (messages) message("Remove: ",names(curp))
    cancel(M) <- var1
  }
  return(M)
}

backwardsearch <- function(x,k=1,...) {
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
      I <- information(Model(x), p=p1, n=x$data$n, data=model.frame(x))
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

forwardsearch <- function(x,k=1,silent=FALSE,...) {
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

  if (!silent)
    cat("Calculating score test for",ncol(restrictedcomb), "models:\n")
  count <- 0
  for (i in 1:ncol(restrictedcomb))
    {
      if (!silent) {
        cat(".")
        count <- count+1
        if (count==20) {          
          cat("\n")
          count <- 0
        }
      }    
      varlist <- c()
      altmodel <- cur ## HA: altmodel, H0: cur
      for (j in 1:k) {
        myvar <- restricted[restrictedcomb[j,i],]
        if (any(wx <- V[myvar]%in%X)) {
          altmodel <- regression(altmodel,V[myvar][which(!wx)],V[myvar][which(wx)])
##          exogenous(altmodel,xfree=FALSE) <- setdiff(X,V[myvar])
        } else {
          covariance(altmodel) <- V[myvar]
        }
        varlist <- rbind(varlist, V[myvar])
      }
      altmodel$parpos <- NULL
##      altmodel <- updatelvm(altmodel,deriv=TRUE,zeroones=TRUE,mean=FALSE)
      altmodel <- updatelvm(altmodel,deriv=TRUE,zeroones=TRUE,mean=TRUE)
      cc <- coef(altmodel, mean=TRUE,silent=TRUE,symbol=c("->","<->"))
      cc0 <- coef(cur, mean=TRUE,silent=TRUE,symbol=c("->","<->"))
      ##          pos <- match(paste(V[i], "<->", V[j], sep=""),cc)
      ##          if (is.na(pos)) ## Should not be necessary 
      ##            pos <- match(paste(V[i], "<->", V[j], sep=""),cc)
      ###      p1 <- numeric(length(pp$p)+k) ## New parameter basically p1 = (p0, 0)
      p1 <- numeric(length(p)+k)
      ## Need to be sure we place 0 at the correct position
      for (ic in 1:length(cc)) {
        idx <- match(cc[ic],cc0)
        if (!is.na(idx))
##            p1[ic] <- pp$p[idx]
          p1[ic] <- p[idx]
      }
#      Sc2 <- score(altmodel,p=p1,S=S,data=NULL,n=n)
      Sc2 <- score(altmodel,p=p1,data=model.frame(x),model=x$estimator,weight=Weight(x))
##      browser()
##      rmidx <- NULL
##      if (!is.null(pp$meanpar)) {
##        rmidx <- 1:length(pp$meanpar)
##      } 
      ##      myidx <- (length(Sc2)-altmodel$index$npar+1):length(Sc2)
##      Sc2 <- Sc2[-rmidx]
##      I <- information(altmodel,p1,n=x$data$n,data=NULL) ##[-rmidx,-rmidx]
      I <- information(altmodel,p1,n=x$data$n,data=model.frame(x),weight=Weight(x),estimator=x$estimator) ##[-rmidx,-rmidx]
      iI <- try(Inverse(I), silent=TRUE)
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

  if (!silent)
    cat("\n")
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
  return(res)
}


print.modelsearch <- function(x,tail=nrow(x$res),adj=c("holm","BH"),...) {
  N <- nrow(x$res)
  if (!is.null(adj)) {
    ##    adjp <- rev(holm(as.numeric(x$test[,2])))
    adjp <- sapply(adj,function(i) p.adjust(x$test[,2],method=i))
    x$res <- cbind(x$res,formatC(adjp))
  }
  print(x$res[(N-tail+1):N,], quote=FALSE, ...)
  invisible(x)
}



