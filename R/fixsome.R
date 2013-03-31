
##' @export
fixsome <- function(x, exo.fix=TRUE, measurement.fix=TRUE, S, mu, n, data, x0=FALSE, na.method="complete.obs", param=lava.options()$param,...) {

  if (is.null(param)) {
    param <- "none"
  } else {
    paramval <- c("hybrid","relative","none","absolute")
    param <- agrep(param,paramval,max.distance=0,value=TRUE)
  }
  
  if (is.character(measurement.fix)) {
    param <- measurement.fix
    measurement.fix <- TRUE
  }
  var.missing <- c()
  if (!missing(data) | !missing(S)) {
        
    if (!missing(data)) {
      dd <- procdata.lvm(x,data=data,na.method=na.method)
    } else {
      dd <- procdata.lvm(x, list(S=S,mu=mu,n=n))
    }
    S <- dd$S; mu <- dd$mu; n <- dd$n    
    var.missing <- setdiff(index(x)$manifest,colnames(S))
  } else { S <- NULL; mu <- NULL }
  
  if (measurement.fix & param!="none") {
    if (length(var.missing)>0) {## Convert to latent:
      new.lat <- setdiff(var.missing,latent(x))
      if (length(new.lat)>0)
      x <- latent(x, new.lat)
    }
    etas <- latent(x)
##    etas <- index(x)$latent
##    ys <- index(x)$endogenous
    ys <- endogenous(x)
    M <- x$M

    for (e in etas) { ## Makes sure that at least one arrow from latent variable is fixed (identification)
      ys. <- names(which(M[e,ys]==1))
      if (length(ys.)>0) {      
        if (tolower(param)=="absolute") {
          if (is.na(intercept(x)[[e]])) intercept(x,e) <- 0
          if (is.na(x$covfix[e,e]) & is.na(x$covpar[e,e])) covariance(x,e) <- 1
        } else {        
          if (param=="hybrid") {
            if (is.na(intercept(x)[[e]])) intercept(x,e) <- 0
###            if (all(is.na(x$fix[e, ys.]==1)) &
            if (all(is.na(x$fix[e, ]==1)) &
                is.na(x$covpar[e,e]) & is.na(x$covfix[e,e])) 
              regfix(x,from=e,to=ys.[1]) <- 1
          } else { ## relative
###           if (all(is.na(x$fix[e, ys.]==1)) &
            if (all(is.na(x$fix[e, ]==1)) &
                is.na(x$covpar[e,e]) & is.na(x$covfix[e,e])) 
              regfix(x,from=e,to=ys.[1]) <- 1
            if (!any(unlist(lapply(intercept(x)[ys.],is.numeric))) &
                ##is.na(intercept(x)[[ys.[1]]]) &
                is.na(intercept(x)[[e]]))
              intercept(x,ys.[1]) <- 0
          }
        }
      }
    }

    ## latintNA <- unlist(lapply(intfix(x)[latent(x)],is.na))
    ## if (length(latintNA)>0) {
    ##   if (any(latintNA))
    ##   intfix(x, latent(x)[which(latintNA)]) <- 0 ## For identifiality we fix mean of latent variables to zero unless already fixed
    ## }
  }

  if (is.null(S)) x0 <- TRUE
  if (exo.fix) {
    if (x0) {
      S0 <- diag(length(index(x)$manifest))
      mu0 <- rep(0,nrow(S0))      
    }
    else {
      S0 <- S
      mu0 <- mu
    }
##    exo <- exogenous(x);
    exo.idx <- index(x)$exo.obsidx;
    ##exo.idx_match(exo,manifest(x)); exo_all.idx <- match(exo, vars(x))
    exo_all.idx <- index(x)$exo.idx

    if (length(exo.idx)>0) {
      for (i in 1:length(exo.idx))
        for (j in 1:length(exo.idx)) {
          i. <- exo_all.idx[i]; j. <- exo_all.idx[j]
          myval <- S0[exo.idx[i],exo.idx[j]];          
          if (i.==j. & myval==0) {
            warning("Overparametrized model. Problem with '"%+%index(x)$vars[j.]%+%"'")
            myval <- 1
          }
          else if (is.na(myval) || is.nan(myval)) myval <- 0
          x$covfix[i.,j.] <- x$covfix[j.,i.] <- myval
        }
      x$mean[exo_all.idx] <- mu0[exo.idx]
    }
  }
  
  index(x) <- reindex(x)  
  return(x)
}
