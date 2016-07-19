cvCheck <- function (x, ...) {
  UseMethod("cvCheck", x)
}

#' @title Test the sensibility of the lvm estimate to the initialization points
#
#' @param object a lvm model
#' @param data a data frame
#' @param factor.vcov inflation factor for the variance when sampling the initialization points
#' @param n.init number of initialization points to be used
#' @param ncpus the number of CPU to be used
#' @param verbose should a progression bar be displayed?
#' @param ... additional arguments to be passed to estimate
#' 
#' @details 
#' Simulation is based on a multivariate truncated normal law (even though it is not satifying for the variance components)
#' 
#' @return a data frame/cvlvm object containing the convergence status (by default 0 indicates successful convergence, see ?optim), the value of the log-likelihood and the estimated parameters (in columns) for each initialization (in rows)
#' 
#' @examples #' 
#' m <- lvm(list(y~v1+v2+v3+v4,c(v1,v2,v3,v4)~x))
#' covariance(m) <- v1~v2+v3+v4
#' dd <- sim(m,10000) ## Simulate 10000 observations from model
#' e <- estimate(m, dd) ## Estimate parameters
#' 
#' system.time(
#' summary(cvCheck(m, dd))
#' )
#' 
#' system.time(
#' summary(cvCheck(m, dd, ncpus = 4))
#' )

cvCheck.lvm <- function(object, data, factor.vcov = 9, n.init = 100, ncpus = 1, verbose = TRUE, ...){
  
  require(lava)
  require(tmvtnorm)
  
  if(is.null(ncpus)){ ncpus <- parallel::detectCores()}
  
  dots <- list(...)
  if("control" %in% names(dots) == FALSE){
    dots$control <- list()
  }
  
  ## automatic intialisation
  test.W <- 0
  
  if(verbose){cat("* initialisation \n")}
  suppressWarnings(
    lvm.init <- estimate(object, data, ...)  
  )
  
  names.coef <- names(coef(lvm.init))
  n.coef <- length(names.coef)
  var.param <- grep(",",names.coef, fixed = TRUE)
  mean.param <- setdiff(1:n.coef, var.param)
  
  VCOV <- vcov(lvm.init)
  diag(VCOV) <-  diag(VCOV)*factor.vcov
  sample.start <- tmvtnorm:::rtmvnorm(n = n.init, mean = coef(lvm.init), sigma = VCOV,
                                      lower = c(rep(-Inf, length(mean.param)), rep(0, length = length(var.param))),
                                      algorithm = "gibbs"
  )
  
  ## grid
  warper <- function(x){
    if(verbose){
      if(ncpus == 1){
        utils:::setTxtProgressBar(.pb, x) 
      }else{
        snowfall::sfCat("*")
      }
    }
    snowfall::sfCat("*")
    dots$control$start <- sample.start[x,]
    
    suppressWarnings(
      resStart <- try(do.call("estimate", 
                              c(list(x = object, data = data), dots))
      )
    )
    
    if(class(resStart) != "try-error"){
      return(c(resStart$opt$convergence, as.numeric(logLik(resStart)), coef(resStart) ))
    }else{
      return(rep(NA, n.coef+2))
    }
  }
  
  if(ncpus == 1){
    if(verbose){.pb <- utils:::txtProgressBar(min = 0, max = n.init, style = 3)}
    Mres <- sapply(1:n.init, warper)
  }else{
    require(snowfall)
    sfInit( parallel=TRUE, cpus=ncpus )
    sfLibrary(lava) ; sfLibrary(snowfall)
    sfExportAll()
    
    Mres <- sfSapply(1:n.init, warper)
    
    sfStop()
  }
  if(verbose){cat("\n")}
  
  Mres <- cbind(Mres,
                c(lvm.init$opt$convergence, as.numeric(logLik(lvm.init)), coef(lvm.init) ))  
  df.resCV <- setNames(as.data.frame(t(Mres)), c("cv", "logLik",names.coef))
  
  ## export
  class(df.resCV) <- c("cvlvm","data.frame")
  return(df.resCV)
}

#' @title Summary function associated with cvCheck
#
#' @param object the output of cvCheck

summary.cvlvm <- function(object, threshold = NULL){
  
  n.init <- length(object$cv)
  pc.cv <- mean(object$cv==0, na.rm = TRUE)
  cat("Convergence rate: ",round(100*pc.cv,2),"% (",sum(object$cv==0, na.rm = TRUE)," over ",n.init,") \n", sep = "")
  
  if(pc.cv>0){
    
    index.cv <- which(object$cv==0)
    index.bestcv <- which.max(object$logLik)
    if(is.null(threshold)){
      threshold <- min(apply(object[index.cv, -(1:2)],2,median))/1e3 # median / 1e3 of the smaller parameter
    }
    dist.coef <- dist(na.omit(object[index.cv, -(1:2)]))
    
    hclust.res <- hclust(dist.coef, method="complete")
    n.clusters <- 1 + sum(hclust.res$height > threshold)
    quantile.cv <- apply(object[index.cv,-1], 2, function(x){ c(range = diff(range(x, na.rm = TRUE)), min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))})
    
    indexRed.bestcv <- which(index.bestcv == na.omit(cbind(object,index = 1:n.init)[index.cv, -(1:2)])$index)
    dist.coef <- sum(as.matrix(dist.coef)[indexRed.bestcv,] < threshold)
    
    cat("Summary of the estimates: \n", sep = "")
    print(quantile.cv)
    cat("Threshold: ",threshold," \n", sep = "")
    cat("Number of convergence points: ",n.clusters," \n", sep = "")
    cat("Convergence rate to the best convergence point: ",round(100*dist.coef/n.init,2),"% (",dist.coef," over ",n.init,") \n", sep = "")
    print(object[index.bestcv,])
    
  }else{
    quantile.cv <- NA
    n.clusters <- NA
  }
  
  return(invisible(list(quantile = quantile.cv,
                        pc.cv = pc.cv,
                        n.clusters = n.clusters)))
}



optimx1 <- function(start,objective,gradient,hessian,...) {
  require(optimx)
  
  nlminbcontrols <- c("eval.max","iter.max","trace","abs.tol","rel.tol","x.tol","step.min","optim.method","optimx.method")
  dots <- list(...)
  control <- list(...)$control
  control <- control[names(control)%in%nlminbcontrols]
  dots$control <- control
  if (length(dots$trace)>0 && dots$trace>0) cat("\n")
  mypar <- c(list(start=start,objective=objective,gradient=gradient,hessian=hessian),dots)
  mypar["debug"] <- NULL
  
  if(!is.null(mypar$control$optimx.method)){
    optimx.method <- mypar$control$optimx.method
    mypar$control$optimx.method <- NULL
    
    optimx.res <- optimx(par = mypar$start, 
                         fn = mypar$objective, 
                         gr = mypar$gradient,
                         method = optimx.method,
                         control = mypar$control)
    
  }else{
    mypar$control$all.methods <- TRUE
    
    optimx.res <- optimx(par = mypar$start, 
                         fn = mypar$objective, 
                         gr = mypar$gradient,
                         control = mypar$control)  
  }
  
  
  index.cv <- which(optimx.res$convcode == 0)
  if(length(index.cv) == 0){
    res <- list(par = rep(NA, length(mypar$start)),
                objective = NA,
                convergence = 1,
                evaluations = c("function" = NA, "gradient" = NA),
                message = "optimx has not converged \n"
    )
  }else if(length(index.cv) == 1){
    
    par_tempo <- unlist(lapply(1:length(mypar$start), function(x){optimx.res[[x]][index.cv]})) 
    names(par_tempo) <- names(mypar$start)
    
    res <- list(par = par_tempo,
                objective = optimx.res$value[index.cv],
                convergence = optimx.res$convcode[index.cv],
                evaluations = c("function" = optimx.res$fevals[index.cv], "gradient" = optimx.res$gevals[index.cv]),
                message = paste0("One algorithm has converged: ",attributes(optimx.res)$row.names[index.cv], " (value = ",optimx.res$value[index.cv],")\n")
    )
  }else{
    index.cvBest <- index.cv[which.max(optimx.res$value[index.cv])]
    par_tempo <- unlist(lapply(1:length(mypar$start), function(x){optimx.res[[x]][index.cvBest]})) 
    names(par_tempo) <- names(mypar$start)
    
    res <- list(par = par_tempo,
                objective = optimx.res$value[index.cvBest],
                convergence = optimx.res$convcode[index.cvBest],
                evaluations = c("function" = optimx.res$fevals[index.cvBest], "gradient" = optimx.res$gevals[index.cvBest]),
                message = paste0("Several algorithms have converged but only the best convergence point is retained \n",
                                 "method: ",paste(attributes(optimx.res)$row.names[index.cv], collapse = " | "), "\n",
                                 "value : ",paste(optimx.res$value[index.cv], collapse = " | "), "\n")
                
    )
  }
  
  return(res)
  
}
