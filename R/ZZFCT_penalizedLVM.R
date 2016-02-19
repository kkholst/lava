#### A- lava functions #### 

#### > penalised_objective.lvm ####
penalised_method.lvm <- lava:::gaussian_method.lvm

penalised_objective.lvm <- function(x, ...){  # proportional to the log likelihood

  obj.UP <- lava:::gaussian_objective.lvm(x,...) # log likelihood
  
  dots <- list(...)
  penalisation <- dots$penalisation
  
  if(!is.null(penalisation)){
    
    obj.P <- penalisation$fn(coef = dots$p[penalisation$names.coef],
                             lambda1 = penalisation$lambda1, 
                             lambda2 = penalisation$lambda2)
   
  }else{
    obj.P <- 0
  }
 
  if(!is.null(penalisation) && penalisation$trace.objective){
    cat("> objectif (UP,P,all): ",obj.UP," , ",obj.P," , ",obj.UP+obj.P,"\n")
    print(dots$p)
    cat("\n")
  }
  
  return( obj.UP + obj.P )
}

#### > penalised_gradient.lvm ####

# penalised_gradient.lvm <- lava:::gaussian_gradient.lvm

penalised_gradient.lvm <- function(x, ...){
  
  grad.UP <- lava:::gaussian_gradient.lvm(x,...)

  dots <- list(...)
  penalisation <- dots$penalisation
  
  if(!is.null(penalisation)){
    res_tempo <- penalisation$gn(coef = dots$p[penalisation$names.coef],
                                 lambda1 = penalisation$lambda1, 
                                 lambda2 = penalisation$lambda2)

  grad.P <- setNames(rep(0, length(dots$p)), names(dots$p))
  grad.P[penalisation$names.coef] <- res_tempo
  ## correction suggested by FU 1998 Penalized Regressions: The Bridge Versus the Lasso
  if(!is.null(penalisation$correction) && penalisation$correction == TRUE){
  grad.P[abs(grad.UP) < dots$tol.grad_pen] <- 0 
  }
  
  }else{
    grad.P <- 0
  }
  
  if(!is.null(penalisation) && penalisation$trace.gradient){
    cat("gradient (UP,P,all): \n ",grad.UP,"\n ",grad.P,"\n ",grad.UP+grad.P,"\n \n")
  }
 
  return( grad.UP + grad.P )

}

#### > penalised_hessian.lvm ####

# penalised_hessian.lvm <- lava:::gaussian_hessian.lvm

penalised_hessian.lvm <- function(x, ...){ # second order partial derivative

  hess.UP <- lava:::gaussian_hessian.lvm(x, ...)
   
  dots <- list(...)
  penalisation <- dots$penalisation
  
  if(!is.null(penalisation)){
    
    res_tempo <- penalisation$hn(coef = dots$p[penalisation$names.coef],
                                 lambda1 = penalisation$lambda1, 
                                 lambda2 = penalisation$lambda2)
    
    hess.P <- setNames(rep(0, length(dots$p)), names(dots$p))
    hess.P[penalisation$names.coef] <- res_tempo
    hess.P <- diag(hess.P)
    
  }else{
    hess.P <- 0
  }
  
  if(!is.null(penalisation) && penalisation$trace.hessian){
    cat("hessian (UP,P,all): ",hess.UP," , ",hess.P," , ",hess.UP+hess.P,"\n \n")
  }
  
  # browser()
  return( hess.UP + hess.P )
  
}

#### > penalised_logLik.lvm ####

# penalised_logLik.lvm <- lava:::gaussian_logLik.lvm
  
penalised_logLik.lvm <- function(object, ...){ # log likelihood
 
  logLik.UP <- lava:::gaussian_logLik.lvm(object,...)
  
  dots <- list(...)
  penalisation <- dots$penalisation
  
  if(!is.null(penalisation)){
    logLik.P <- penalisation$fn(coef = dots$p[penalisation$names.coef],
                                lambda = penalisation$lambda, 
                                power = penalisation$power)
   
  }else{
    logLik.P <- 0
  }
  
  # browser()
  return( logLik.UP + logLik.P )
}

#### B- extra functions ####

#### > ls.fct_penalty ####

ls.fct_penalty <- function(lvm, lambda1 = 0, lambda2 = 0, names.coef = NULL,
                           trace.objective = FALSE, trace.gradient = FALSE, trace.hessian = FALSE){
  
  
  ## preparation
  if(is.null(names.coef)){
    #     formula.lvm <- formula(lvm)[[1]]
    #     names.coef <- c(all.vars(formula.lvm)[1],
    #                     paste(all.vars(formula.lvm)[1],all.vars(formula.lvm)[-1], sep = "~")
    #     )

    names.coef <- paste(lvm$index$endogenous ,lvm$exogenous , sep = "~")
    
  }
  
  ## main
  fn_penalty <- function(coef, lambda1, lambda2){
    fn1 <- lambda1 * sum( abs(coef) )
    fn2 <- lambda2/2 * sum( coef^2 )
    
    return( fn1 + fn2 )
  }
  
  gn_penalty <- function(coef, lambda1, lambda2){
    gn1 <- lambda1 * sign(coef)
    gn2 <- lambda2 * abs(coef)
    
    return( gn1 + gn2 )
  }
  
  hn_penalty <- function(coef, lambda1, lambda2){
    hn1 <- 0
    hn2 <- lambda2 * sign(coef)
    
    return( hn1 + hn2 )
    
  }
  
  ## export
  return(list(fn_penalty = fn_penalty,
              gn_penalty = gn_penalty,
              hn_penalty = hn_penalty,
              names.coef = names.coef,
              lambda1 = lambda1, 
              lambda2 = lambda2,
              trace.objective = trace.objective, 
              trace.gradient = trace.gradient, 
              trace.hessian = trace.hessian))
}


