#### 1- S3 Methods ####

# ##' @export
`penalize` <-
  function(x,...) UseMethod("penalize")

##' @export
"penalize<-" <- function (x, ..., value) {
  UseMethod("penalize<-", x)
}


#### 2- Penalty functions ####
`penalize.lvm` <- function(x, ...){
  
  dots <- list(...)
  if("value" %in% names(dots)){
    value <- dots$value
    dots$value <- NULL
    penalize(x, dots) <- dots$value
  }else{
    penalize(x, ...) <- NULL
  }
  
  return(x)
}
`penalize.plvm` <- `penalize.lvm`

`penalize<-.lvm` <- function(x, ..., value){
  
  ## add slot to object
  if("plvm" %in% class(x) == FALSE){
    x$penalty <- list(fn_penalty = NULL,
                      gn_penalty = NULL,
                      hn_penalty = NULL,
                      names.coef = NULL,
                      index.coef = NULL,
                      lambda1 = 0, 
                      lambda2 = 0,
                      trace.objective = FALSE, 
                      trace.gradient = FALSE, 
                      trace.hessian = FALSE)
    class(x) <- append("plvm", class(x))
  }
  
  ## main (call `penalty<-.plvm`)
  penalize(x, ...) <- value
  # `penalize<-.plvm`(x, ..., value = value)
  # do.call(`penalize<-.plvm`, args = list(x = x, ..., value = value))
 
  ## export
  return(x)
}

`penalize<-.plvm` <- function(x, lambda1, lambda2, fn_penalty, gn_penalty, hn_penalty, ..., value){
  
  ## coefficients
  if(!is.null(value)){
    
    if(any(value %in% coef(x) == FALSE)){
      stop("penalty<-.lvm: coefficients to be penalized do not match those of the model\n",
           "unknown coefficients: ",paste(value[value %in% coef(x) == FALSE], collapse = " "),"\n",
           "available coefficients: ",paste(coef(x)[coef(x) %in% value == FALSE], collapse = " "),"\n")
    }
    
    x$penalty$names.coef <- value
    x$penalty$index.coef <- which(coef(x) %in% names.coef)
  } else if(is.null(x$penalty$names.coef)){
    x$penalty$names.coef <- paste(x$index$endogenous ,x$exogenous , sep = "~")
    x$penalty$index.coef <- which(coef(x) %in% x$penalty$names.coef)
  } 
  
  ## penalization parameters
  if(!missing(lambda1)){
    x$penalty$lambda1 <- lambda1
  }
  
  if(!missing(lambda1)){
    x$penalty$lambda2 <- lambda2
  }
  
  ## objective function
  if(!missing(fn_penalty)){
    
    x$penalty$fn_penalty <- fn_penalty
    
  } else if(is.null(x$penalty$fn_penalty)){
    
    x$penalty$fn_penalty <-  fn_penalty <- function(coef, lambda1, lambda2){
      fn1 <- lambda1 * sum( abs(coef) )
      fn2 <- lambda2/2 * sum( coef^2 )
      
      return( fn1 + fn2 )
    }
    
  }
  
  ## gradient function
  if(!missing(gn_penalty)){
    
    x$penalty$gn_penalty <- gn_penalty
    
  } else if(is.null(x$penalty$gn_penalty)){
    
    x$penalty$gn_penalty <- function(coef, lambda2){
      return( lambda2 * abs(coef) )
    }
    
  }
  
  ## hessian function
  if(!missing(hn_penalty)){
    
    x$penalty$hn_penalty <- hn_penalty
    
  } else if(is.null(x$penalty$hn_penalty)){
    
    x$penalty$hn_penalty <- hn_penalty <- function(coef, lambda2){
      return( lambda2 * sign(coef) )
    }
    
  }
  
  ## dots
  dots <- list(...)
  names.dots <- names(dots)
  
  if(length(names.dots) > 0){
    
    
    if(any(names.dots %in% names(x$penalty) == FALSE)){
      fixedArgs <- c("names.coef", "index.coef",
                     setdiff(names(formals("penalty<-.lvm")),c("..."))
      )
      
      stop("penalty<-.lvm: some additional arguments are invalid\n",
           "invalid arguments: ",paste(names.dots[names.dots %in% names(x$penalty) == FALSE], collapse = " "),"\n",
           "valid additional arguments: ", paste(setdiff(names(x$penalty), fixedArgs), collapse = " "),"\n")
      
    }
    x$penalty[names(dots)] <- dots
  }
  
  ## export
  return(x)
}

#### 3- optim functions #### 

penalized_method.lvm <- "nlminb2"#lava:::gaussian_method.lvm

penalized_objective.lvm <- function(x, ...){  # proportional to the log likelihood
  
  obj.UP <- lava:::gaussian_objective.lvm(x,...) # log likelihood
  
  dots <- list(...)
  penalty <- dots$penalty
  
  if(penalty$lambda1>0 || penalty$lambda2>0){ # L1 or L2 penalty
    
    obj.P <- penalty$fn(coef = dots$p[penalty$index.coef],
                        lambda1 = penalty$lambda1,
                        lambda2 = penalty$lambda2)
    
  }else{
    obj.P <- 0
  }
  
  return( obj.UP + obj.P )
}

penalized_gradient.lvm <- function(x, ...){
  
  grad.UP <- lava:::gaussian_gradient.lvm(x,...)
  
  dots <- list(...)
  penalty <- dots$penalty
  
  if(penalty$lambda2>0){ # only L2 penalty
    res_tempo <- penalty$gn(coef = dots$p[penalty$index.coef],
                                 lambda2 = penalty$lambda2)
    
    grad.P <- rep(0, length(dots$p))
    grad.P[penalty$index.coef] <- res_tempo
    
  }else{
    grad.P <- rep(0, length(grad.UP))
  }

  return( grad.UP + grad.P )
  
}

penalized_hessian.lvm <- function(x, ...){ # second order partial derivative
  
  hess.UP <- lava:::gaussian_hessian.lvm(x, ...)
  
  dots <- list(...)
  penalty <- dots$penalty
  
  if(penalty$lambda2>0){ # only L2 penalty
    
    res_tempo <- penalty$hn(coef = dots$p[penalty$index.coef],
                                 lambda2 = penalty$lambda2)
    
    hess.P <- rep(0, length(dots$p))
    hess.P[penalty$index.coef] <- res_tempo
    hess.P <- diag(hess.P)
    
  }else{
    hess.P <- 0
  }
  
  return( hess.UP + hess.P )
  
}

penalized_logLik.lvm <- function(object, ...){ # log likelihood
  
  logLik.UP <- lava:::gaussian_logLik.lvm(object,...)
  
  dots <- list(...)
  penalty <- dots$penalty
  
  if(!is.null(penalty)){
    logLik.P <- penalty$fn(coef = dots$p[penalty$names.coef],
                                lambda = penalty$lambda, 
                                power = penalty$power)
    
  }else{
    logLik.P <- 0
  }
  
  # browser()
  return( logLik.UP + logLik.P )
}

