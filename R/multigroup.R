###{{{ multigroup

multigroup <- function(models, datasets, fix, exo.fix=TRUE, keep=NULL, missing=FALSE, ...) {
  nm <- length(models)
  if (nm!=length(datasets)) stop("Supply dataset for each model")
  if (nm<2) stop("Two or more groups neeeded")
  mynames <- names(models)
  
  ## Check for random slopes
  xfix <- list()
  for (i in 1:nm) {
    x0 <- models[[i]]
    data0 <- datasets[[i]]
    xfix0 <- colnames(data0)[(colnames(data0)%in%parlabels(x0,exo=TRUE))]
    xfix <- c(xfix, list(xfix0))
  }
  if (missing(fix)) {
    fix <- !any(unlist(lapply(xfix, function(x) length(x)>0)))  
  }
  for (i in 1:nm) {
    x0 <- models[[i]]
    data0 <- datasets[[i]]
    if (length(exogenous(x0)>0)) {
      catx <- categorical2dummy(x0,data0)
      models[[i]] <- catx$x; datasets[[i]] <- catx$data
    }
  }

  
  models.orig <- NULL
######################
### MLE with MAR mechanism
######################
  if (missing) {

    parcount <- 0
    reservedpars <- c()
    mynpar <- c()
    for (i in 1:nm) {
      ## Fix some parameters (predictors,latent variables,...)
      d0 <- datasets[[i]][1,,drop=FALSE]; d0[,] <- 1
      models[[i]] <- fixsome(models[[i]], exo.fix=exo.fix, measurement.fix=fix, data=d0)
      ## Find named/labelled parameters
      rpar <- unique(parlabels(models[[i]]))
      reservedpars <- c(reservedpars, rpar)
      mynpar <- c(mynpar, with(index(models[[1]]), npar+npar.mean))      
    }; reservedpars <- unique(reservedpars)
    nonamepar <- sum(mynpar)
    ## Find unique parameter-names for all parameters
    newpars <- c()
    i <- 0
    pos <- 1
    while(pos<=nonamepar) {
      i <- i+1
      newname <- paste("par",i,sep="")
      if (!(newname%in%reservedpars)) {
        newpars <- c(newpars,newname)
        pos <- pos+1
      }
    } 

    pos <- 0
    models0 <- list()
    datasets0 <- list()
    for (i in 1:nm) {
      myvars <- unlist(intersect(colnames(datasets[[i]]),c(vars(models[[i]]),xfix[[i]],keep)))
      mydata <- datasets[[i]][,myvars]
      if (any(is.na(mydata))) {
        if (i>1) pos <- pos+mynpar[i-1]
        models[[i]] <- baptize(models[[i]],newpars[pos+1:mynpar[i]] ,overwrite=FALSE)        
        ##        warning("Missing data encountered")  
        val <- missingModel(models[[i]],mydata,fix=FALSE,keep=keep)
        datasets0 <- c(datasets0, val$datasets)
        models0 <- c(models0, val$models)            
      } else {
        datasets0 <- c(datasets0, list(mydata))
        models0 <- c(models0, list(models[[i]]))
      }
    }

    models.orig <- models
    val <- multigroup(models0,datasets0,fix=FALSE,missing=FALSE,...)
    val$models.orig <- models.orig; val$missing <- TRUE
    return(val)
  }

  
######################
### Usual analysis:
######################
  for (i in 1:nm) {
    if (is.data.frame(datasets[[i]])) {
      myvars <- intersect(names(datasets[[i]]),c(vars(models[[i]]),xfix[[i]],keep))
      if (any(is.na(datasets[[i]][,myvars]))) {
        warning(paste("Missing data encountered (group ",i,")... Going for complete-case analysis", sep=""))
        datasets[[i]] <- na.omit(datasets[[i]][,myvars,drop=FALSE])
      }
    }
  }
    
  ##
  exo <- exogenous(models)
  means <- lvms <- As <- Ps <- ps <- datas <- samplestat <- list()
  for (i in 1:nm) {

    if (!is.null(exogenous(models[[i]]))) {
      if (any(is.na(exogenous(models[[i]])))) {
        exogenous(models[[i]]) <- exo
      }
    }
    ##    if (!is.null(exogenous(models[[i]])) || is.na(exogenous(models[[i]])))
    ##exogenous(models[[i]]) <- exo
    
    ##    mydata <- datasets[[i]][,manifest(models[[i]]),drop=FALSE]
    mydata <- datasets[[i]]
    mymodel <- fixsome(models[[i]], data=mydata, measurement.fix=fix, exo.fix=exo.fix)
    ##index(mymodel) <- reindex(mymodel,zeroones=TRUE,deriv=TRUE)
    mymodel <- updatelvm(mymodel,zeroones=TRUE,deriv=TRUE)
    
    P <- index(mymodel)$P1; P[P==0] <- NA
    P[!is.na(P) & !is.na(mymodel$covpar)] <- mymodel$covpar[!is.na(P) & !is.na(mymodel$covpar)]

    A <- index(mymodel)$M1; A[A==0] <- NA
    A[!is.na(A) & !is.na(mymodel$par)] <- mymodel$par[!is.na(A) & !is.na(mymodel$par)]
    p <- pars(mymodel, A, P)
    p[p=="1"] <- NA
    
    mu <- unlist(mymodel$mean)[which(index(mymodel)$v1==1)]
    means <- c(means, list(mu))
    lvms <- c(lvms, list(mymodel))
    datas <- c(datas, list(mydata))
    samplestat <- c(samplestat, list(procdata.lvm(models[[i]],data=mydata)))
    As <- c(As, list(A))
    Ps <- c(Ps, list(P))
    ps <- c(ps, list(p))
  };

######
  pp <- unlist(ps)
  parname <- unique(pp[!is.na(pp)])
  opt <- options(warn=-1); pidx <- is.na(as.numeric(parname)); options(opt)
  parname <- unique(unlist(pp[!is.na(pp)]));
  nfree <- sum(is.na(pp)) + length(parname)
  
  if (nfree>0) {
    pp0 <- lapply(ps, is.na)
    usedname <- cbind(parname, rep(NA,length(parname)))
    counter <- 1
    pres <- pp0
    for (i in 1:length(pp0)) {
      if (length(pp0[[i]]>0))
      for (j in 1:length(pp0[[i]])) {
        pidx <- match(ps[[i]][j],parname)
        if (pp0[[i]][j]) {
          pres[[i]][j] <- paste("p",counter,sep="")
          counter <- counter+1        
        } else if (!is.na(pidx)) {
          if (!is.na(usedname[pidx,2])) {
            pres[[i]][j] <- usedname[pidx,2]
          } else {
            val <- paste("p",counter,sep="")
            pres[[i]][j] <- val
            usedname[pidx,2] <- val
            counter <- counter+1
          }
        } else {
          pres[[i]][j] <- NA
        }
      }
    }
    mypar <- paste("p",1:nfree,sep="")
    myparpos <- pres
    myparlist <- lapply(pres, function(x) x[!is.na(x)])
  } else {
    mypar <- NULL
    myparpos <- NULL
    myparlist <- NULL
  }

  ### Mean parameter
  
  mm <- unlist(means)
  meanparname <- unique(mm[!is.na(mm)])
  opt <- options(warn=-1); midx <- is.na(as.numeric(meanparname)); options(opt)
  meanparname <- meanparname[midx]
  any.mean <- sum(is.na(mm)) + length(meanparname)
  nfree.mean <- sum(is.na(mm)) + length(setdiff(meanparname,parname))
  
  if (any.mean>0) {
    mm0 <- lapply(means, is.na)
    usedname <- cbind(meanparname, rep(NA,length(meanparname)))
    counter <- 1
    res <- mm0
    for (i in 1:length(mm0)) {
      if (length(mm0[[i]])>0)
      for (j in 1:length(mm0[[i]])) {
        midx <- match(means[[i]][j],meanparname)
        if (mm0[[i]][j]) {
          res[[i]][j] <- paste("m",counter,sep="")
          counter <- counter+1        
        } else if (!is.na(midx)) {
          pidx <- match(meanparname[midx],pp)
          if (!is.na(pidx)) {
            res[[i]][j] <- unlist(myparlist)[pidx]
          } else {
            if (!is.na(usedname[midx,2])) {
              res[[i]][j] <- usedname[midx,2]
            } else {
              val <- paste("m",counter,sep="")
              res[[i]][j] <- val
              usedname[midx,2] <- val
              counter <- counter+1
            }
          }
        } else {
          res[[i]][j] <- NA
        }
      }
    }
    mymeanpos <- res
    mymeanlist <- lapply(res, function(x) x[!is.na(x)])
    mymean <- unique(unlist(mymeanlist))   ##paste("m",1:nfree.mean,sep="")    
  } else {
    mymean <- NULL
    mymeanpos <- NULL
    mymeanlist <- NULL
  }  
  ####

##  mymean <- paste("m",1:nfree.mean,sep=""); names(mymean) <- vars(models)[midx]
  
  res <- list(npar=nfree, npar.mean=nfree.mean, ngroup=length(lvms), names=mynames,
              lvm=lvms, data=datas, samplestat=samplestat,
              A=As, P=Ps,
              meanpar=names(mu),
              par=mypar, parlist=myparlist,  parpos=myparpos,
              mean=mymean, meanlist=mymeanlist, meanpos=mymeanpos,
              models.orig=models.orig, missing=missing
              )
  class(res) <- "multigroup"
  checkmultigroup(res)
  return(res)
}

###}}}

checkmultigroup <- function(x) {
    ## Check validity:
  for (i in 1:x$ngroup) {
    if (nrow(x$data[[i]])<2) {
      warning("With only one observation in the group, all parameters should be inherited from another a group!")
    }
  }
}
