###{{{ missingModel

missingModel <- function(model,data,var=endogenous(model),fix=FALSE,type=2,keep=NULL,weight=NULL,weight2=NULL,cluster=NULL,...) {
  if (class(model)!="lvm") stop("Needs a lvm-object")
  if (type==3) {
    var <- manifest(model)
  }

  data0 <- subset(data, select=manifest(model))
  data.mis <- is.na(data0[,var,drop=FALSE])
  patterns <- unique(data.mis,MARGIN=1)
  
  mis.type <- apply(data.mis,1,
                  function(x) which(apply(patterns,1,function(y) identical(x,y))))
  pattern.allmis <- which(apply(patterns,1,all)) ## Remove entry with all missing
##  if (!is.null(pattern.allmis))
##    patterns <- patterns[-pattern.allmis,]
 
##  if (fix) {
##    data0 <- data[1,,drop=FALSE]; data0[,] <- 1
##    model <- fixsome(model, data=data0, measurement.fix=TRUE, exo.fix=FALSE)
##  }
  
  models <- datasets <- weights <- weights2 <- clusters <- c()
  mymodel <- baptize(model)
  pattern.compl <- 0
  count <- 0
  A <- index(model)$A
  topendo <- endogenous(model,top=TRUE)
  exo <- exogenous(model)
  exclude <- c()

  for (i in setdiff(1:nrow(patterns),pattern.allmis)) {
    includemodel <- TRUE
    count <- count+1
    mypattern <- patterns[i,]
    m0 <- mymodel;
    if (any(mypattern)) {
      ##      browser()
      latent(m0,zero=FALSE) <- colnames(data.mis)[mypattern]
      if (type>1) {
        mytop <- intersect(topendo,colnames(data.mis)[mypattern])
        if (!is.null(mytop))
          kill(m0) <- mytop
        for (xx in exo) {          
          ## If exogenous variable only have effect on missing variables,
          ##  then remove it from the model
          if (all(c(rownames(A)[A[xx,]==1])%in%names(mypattern)[mypattern])) {
            kill(m0) <- xx
            ## and is missing,
##            if (all(c(xx,rownames(A)[A[xx,]==1])%in%names(mypattern)[mypattern])) {

          }
        }
      }
    ##      kill(m0) <- colnames(data.mis)[mypattern]
    } else
    pattern.compl <- count
##    d0 <- data[mis.type==i,manifest(m0),drop=FALSE];
    d0 <- data[mis.type==i,c(manifest(m0),keep),drop=FALSE];

    w0.var <- intersect(manifest(m0),colnames(weight))
    w0 <- weight[mis.type==i,w0.var,drop=FALSE];
    if (!is.list(weight2)) {
      w02.var <- intersect(manifest(m0),colnames(weight2))
      w02 <- weight2[mis.type==i,w02.var,drop=FALSE];
    } else {
      weights2 <- weight2
    }
    clust0 <- cluster[mis.type==i]
    modelexo <- exogenous(model)
    exogenous(m0) <- modelexo
    ##    index(m0) <- reindex(m0,deriv=TRUE,zeroones=TRUE)
    if (is.null(intersect(modelexo,latent(m0)))) {
      print("Missing exogenous variables... Going for complete-case analysis in these cases")
    } else {
      if( sum(unlist(index(m0)[c("npar","npar.mean")]))>0 ) {
        models <- c(models, list(m0))
        datasets <- c(datasets, list(d0))
        weights <- c(weights, list(w0))
        if (!is.list(weight2))
          weights2 <- c(weights2, list(w02))
        clusters <- c(clusters, list(clust0))
      } else {
        exclude <- c(exclude,count)
      }
    }
  }

  
##  print(exclude)
  Patterns <- patterns
  if (length(exclude)>0)
    Patterns <- Patterns[-exclude,]
  pattern.allcomp<- which(apply(Patterns,1,function(x) all(!x))) ## Complete cases 

  res <- list(models=models, datasets=datasets,
              weights=weights,
              weights2=weights2,
              clusters=clusters,
              patterns=Patterns,
              pattern.compl=pattern.compl,
              pattern.allmis=pattern.allmis,
              pattern.allcomp=pattern.allcomp,
              mis.type=mis.type)
  return(res)  
}

###}}}

###{{{ estimate.MAR.lvm

estimate.MAR <- function(x,data,which=endogenous(x),fix,type=2,startcc=FALSE,control=list(),silent=FALSE,weight,weight2,cluster,onlymodel=FALSE,estimator="gaussian",hessian=TRUE,keep=NULL,...) {
  cl <- match.call()

  ##  cl[1] <- call("missingModel")
  ##  val <- eval(cl)
  Debug("estimate.MAR")
  redvar <- intersect(intersect(parlabels(x),latent(x)),colnames(data))
  if (length(redvar)>0 & !silent)
    warning(paste("Remove latent variable colnames from dataset",redvar)) 
  
  xfix <- setdiff(colnames(data)[(colnames(data)%in%parlabels(x,exo=TRUE))],latent(x))
##    names(data)[(names(data)%in%parlabels(x))]
##  if (length(xfix)>0) {
##    warning("Random slopes detected. Estimation with data MAR not supported.")   
##  }
  if (missing(fix))
    fix <- ifelse(length(xfix)>0,FALSE,TRUE)  

  S <- diag(length(manifest(x)));
  mu <- rep(0,nrow(S));  
  ##  S <- mu <- NULL
  K <- length(exogenous(x))
  vnames <- index(x)$manifest
  names(mu) <- rownames(S) <- colnames(S) <- vnames
  if (K>0) {
    ##if (!silent)cat("Calculating 1. and 2. moments of exogenous variables...\n")    
    exo.idx <- c(which(manifest(x)%in%exogenous(x)))
    xx <- subset(Model(x),exogenous(x))
    exogenous(xx) <- NULL
    covfix(xx, vars(xx)) <- NA
    xx <- covariance(xx,exogenous(x),exogenous(x))
    datax <- data[,exogenous(x),drop=FALSE]
    mu0 <- colMeans(datax,na.rm=TRUE)
    cov0 <- cov(datax,use="pairwise.complete.obs")*(nrow(datax)-1)/nrow(datax)
    cov0upper <- cov0[upper.tri(cov0,diag=TRUE)]
    exogenous(xx) <- NULL
    coefpos <- matrices(xx,1:(K*(K-1)/2+K))$P
    ii <- coefpos[upper.tri(coefpos,diag=TRUE)]
    start <- c(mu0, cov0upper[order(ii)])
    ##    ex <- estimate(xx,datax,silent=TRUE)
    S[exo.idx,exo.idx] <- cov0
    mu[exo.idx] <- mu0    
    ##    cat("\n")
  }

  x0 <- x
  x <- fixsome(x, measurement.fix=fix, exo.fix=TRUE, S=S, mu=mu, n=1)
  
  ##  print(covfix(x)$values)##[,exogenous(x),exogenous(x)])
  ##  data0 <- data[1,,drop=FALSE]; data0[,] <- 1
  ##  if (fix) {
  ##    x <- fixsome(x, data=data0, measurement.fix=fix, exo.fix=FALSE)
  ##    x <- fixsome(x, data=data0, S=S, mu=mu, measurement.fix=fix, exo.fix=TRUE, x0=TRUE)
  ##  }    
  ##fixsome(x,data=data0,measurement.fix=FALSE,exo.fix=TRUE)

  if (!silent)
    cat("Identifying missing patterns...")

  val <- missingModel(x,data,var=which,type=type,keep=c(keep,xfix),weight=weight,weight2=weight2,cluster=cluster,...)
  if (!silent)
    cat("\n")

##   res <- c()
  
##   x0 <- val$models[[5]]
##   data0 <- val$data[[5]]

##   x1 <- val$models[[5]]
##   data1 <- val$data[[5]]
##   f1 <- function(p) -logLik(x1,data1,p=p)
##   g1 <- function(p) colSums(-score(x1,data1,p))


##   f <- function(p) -logLik(x0,data0,p=p)
##   g <- function(p) colSums(-score(x0,data0,p))
##   res <- rbind(res,g(p))
  
##   ff <- function(p) f(p)+f1(p)
##   gg <- function(p) g(p)+g1(p)

  
##   lower <- c(-Inf,-Inf,-Inf,-Inf,rep(1e-9,4))
##   nlminb(p,ff,gg,lower=lower,control=list(trace=1))

##   optim(p,f,g,control=list(trace=1))
        
  
##   e0 <- estimate(x0,data0,control=list(trace=1,method="nlminb1"))    

    
  ##S <- cov(data,use="pairwise.complete.obs")
  ##mu <- colMeans(data,na.rm=TRUE) 

  
##   exo <- exogenous(model);
##   obs <- manifest(model)  
##   S <- diag(length(obs))
##   mu <- rep(0,nrow(S))
##   exo.idx <- match(exo,obs);
##   S[exo.idx,exo.idx] <- var(data[,exo,drop=FALSE])
##   colnames(S) <- rownames(S) <- obs
##   mu[exo.idx] <- colMeans(data[,exo,drop=FALSE])
##   names(mu) <- obs
##   model <- fixsome(model, measurement.fix=fix, S=S, mu=mu, n=1)
  if (nrow(val$patterns)==1) {
    res <- estimate(x,data=data,fix=fix,weight=weight,weight2=weight2,estimator=estimator,silent=silent,control=control,...)
    return(res)
  }
  if (startcc & is.null(control$start)) {
    if (!silent)
      cat("Obtaining starting value...")
    start0 <- rep(1,sum(unlist(index(x)[c("npar","npar.mean")])))
    mystart <- tryCatch(
                        (estimate(x,data=na.omit(data),silent=TRUE,
                                     weight=weight,weight2=weight2,estimator=estimator,quick=TRUE,... 
                                      )),
                        error=function(e) rep(1,sum(unlist(index(x)[c("npar","npar.mean")])))
                        )
    control$start <- mystart
    if (!silent)
      cat("\n")
  }
  if (is.null(control$meanstructure))
    control$meanstructure <- TRUE
  if (is.null(control$information))
    control$information <- "obs"

  mg0 <- with(val, suppressWarnings(multigroup(models,datasets,fix=FALSE,exo.fix=FALSE,missing=FALSE)))

  if (onlymodel) return(list(mg=mg0,val=val,weight=val$weights,weight2=val$weights2,cluster=val$clusters))

  
##  e.mis <- estimate(mg0,control=list(start=p,trace=1,method="nlminb1"))
  e.mis <- estimate(mg0,control=control,silent=silent,weight=val$weights,weight2=val$weights2,
                    cluster=val$clusters,estimator=estimator,...)
  ##  return(e.mis)
  ##cc <- coef(e.mis,level=1)[[pattern.compl]]
  cc <- coef(e.mis,level=0)
  mynames <- c()
  if (e.mis$model$npar.mean>0)
  ##   mynames <- c(mynames,coef(x,mean=TRUE)[1:e.mis$model$npar.mean])
    mynames <- c(mynames,paste("m",1:e.mis$model$npar.mean,sep=""))
   if (e.mis$model$npar>0)
     mynames <- c(mynames,paste("p",1:e.mis$model$npar,sep=""))
  rownames(cc) <- mynames 

  
  mycc <- val$pattern.allcomp ## Position of complete-case model
  nmis <- with(val, as.numeric(table(mis.type)[pattern.allmis])) ## Number of completely missing observations
  if (length(nmis)>0 & length(mycc)>0) ## Any individuals with all missing?
    if (val$pattern.allmis<mycc)
      mycc <- mycc-1

  if (length(xfix)>0) {
    nrow <- length(vars(x))
    xpos <- lapply(xfix,function(y) which(regfix(x)$labels==y))
    colpos <- lapply(xpos, function(y) ceiling(y/nrow))
    rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
    myfix <- list(var=xfix, col=colpos, row=rowpos)
    for (i in 1:length(xfix)) 
      ##      for (j in 1:length(myfix$col[[i]])) 
      regfix(x, from=vars(x)[rowpos[[i]]],to=vars(x)[colpos[[i]]]) <-
        rep(colMeans(data[,xfix[i],drop=FALSE],na.rm=TRUE),length(rowpos[[i]]))
    x <-
      updatelvm(x,zeroones=TRUE,deriv=TRUE)
                ##                reindex(x,zeroones=TRUE,deriv=TRUE)
  }

  ord <- c()
  ordlist <- list()
  for (i in 1:nrow(val$patterns)) {
    ordlist <- c(ordlist, list(which(val$mis.type==i)))
    ord <- c(ord, ordlist[[i]])    
  }

  res <- with(val, list(coef=cc,
                        patterns=patterns, table=table(mis.type),
                        mis.type=mis.type,
                        order=ord,
                        orderlist=ordlist,
                        nmis=nmis,
                        allmis=pattern.allmis,
##                        cc=pattern.allcomp,
                        cc=mycc,
                        ncc=as.numeric(table(mis.type)[pattern.allcomp]),
                        multigroup=e.mis$model,
##                        mg0=e.mis$model0$lvm,
                        estimate=e.mis,
                        model=x,
                        model0=x0,
                        vcov=e.mis$vcov, opt=e.mis$opt,
                        control=control,
                        data=list(model.frame=data),
                        estimator=estimator,
                        call=cl
                        ))
  class(res) <- c("lvm.missing","lvmfit") #,"lvmfit")
  if ("lvmfit.randomslope"%in%class(e.mis))
    class(res) <- c(class(res),"lvmfit.randomslope")
  if (hessian & is.null(cluster)) {
    if (!silent)
      cat("Calculating asymptotic variance...\n")
    res$vcov <- solve(information(res$estimate,type="hessian"))
    cc[] <- coef(e.mis,level=0,vcov=res$vcov)
    res$coef <- cc
  }
  
  res <- edgelabels(res,type="est")
  return(res)
}

###}}} estimate.MAR.lvm
