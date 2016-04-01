#' @title Estimate a penalized lvm model
#
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param lambda1 L1 penalization parameter
#' @param lambda2 L2 penalization parameter
#' @param fn_penalty penalty function relative to the objective function
#' @param gn_penalty penalty function relative to the first derivative of the objective function
#' @param hn_penalty penalty function relative to the second derivative of the objective function
#' @param method optimize to be used
#' @param ... additional arguments to be passed to estimate
#' 
#' @details 
#' Remaining issues with ridge regresion + elastic net to check
#' 
#' @return a
#' 
#' @examples #' 
#' set.seed(10)
#' lvm.model <- lvm(list(y~X1+X2+X3+X4))
#' df.data <- sim(lvm.model,300)
#' plvm.model <- penalize(lvm.model)
#' 
#' ## unpenalized
#' lvm.fit <- estimate(lvm.model, df.data)
#' 
#' ## L1 penalisation
#'  
#' p1lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = 0.1,
#'                control = list(constrain = TRUE, iter.max = 1000))
#' p1lvm.fit
#' 
#' ls.p1lvm <- lapply(c(0,0.1,0.2,0.25,0.5), 
#'                  function(x){estimate(plvm.model,  data = df.data, lambda1 = x, control = list(constrain = TRUE, iter.max = 1000))})
#'                
#' lapply(ls.p1lvm,coef)
#' 
#' ## L2 penalisation
#' 
#' p2lvm.fit <- estimate(plvm.model,  data = df.data, lambda2 = 0.1, 
#'                control = list(constrain = TRUE, iter.max = 1000))
#' p2lvm.fit
#' 
#' ls.p2lvm <- lapply(c(0,0.1,0.2,0.25,0.5), 
#'                  function(x){estimate(plvm.model,  data = df.data, lambda2 = x, control = list(constrain = TRUE, iter.max = 1000))})
#'                
#' lapply(ls.p2lvm,coef)
#' 
#' ## elastic net
#' 
#' p12lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = 0.1, lambda2 = 0.1,
#'                control = list(constrain = TRUE, iter.max = 1000))
#' p12lvm.fit
#' 
#' ls.p12lvm <- lapply(c(0,0.1,0.2,0.25,0.5), 
#'                  function(x){estimate(plvm.model,  data = df.data, lambda1 = x, lambda2 = x, control = list(constrain = TRUE, iter.max = 1000))})
#'                
#' lapply(ls.p12lvm,coef)
#' 
#' 
`estimate.plvm` <- function(x, data, lambda1, lambda2, fn_penalty, gn_penalty, hn_penalty, 
                            regularizationPath = FALSE, fix.sigma = FALSE,
                            method = penalized_method.lvm, proxGrad.method = "ISTA", ...) {
  
  ## generic penalty
  penalty  <- x$penalty
  
  ## specific values
  if(!missing(lambda1)){
    penalty$lambda1 <- lambda1
  }
  if(!missing(lambda2)){
    penalty$lambda2 <- lambda2
  }
  if(!missing(fn_penalty)){
    penalty$fn_penalty <- fn_penalty
  }
  if(!missing(gn_penalty)){
    penalty$gn_penalty <- gn_penalty
  }
  if(!missing(hn_penalty)){
    penalty$hn_penalty <- hn_penalty
  }
  
  
  # optimizer
  if(proxGrad.method %in% c("ISTA","FISTA") == FALSE){
    stop("estimate.plvm: wrong specification of control$proxGrad.method \n",
         "only  \"ISTA\" and \"FISTA\" mehtods are available \n",
         "control$proxGrad.method: ",proxGrad.method,"\n")
  }  
  
  ## dots  :  cleaning [not completely checked - may be issues here]
  dots <- list(...)
  
  names.dots <- names(dots)
  if(any(names.dots %in% names(penalty))){
    fixedArgs <- setdiff(names(formals("estimate.plvm")),c("..."))
    dotsArgs <- setdiff(names.dots[names.dots %in% names(penalty)], fixedArgs)
    
    if(length(dotsArgs)>0){
      penalty[dotsArgs] <- dots[dotsArgs]
      dots[dotsArgs] <- NULL
    }
  }
  
  #### main
  res <- do.call(`estimate.lvm`, args = c(list(x = x, data = data, estimator = "penalized", 
                                               penalty = penalty, regularizationPath = regularizationPath, fix.sigma = fix.sigma,
                                               method = method, proxGrad.method = proxGrad.method), 
                                          dots)
  )
  
  #### export
  return(res)
}