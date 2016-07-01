#### OVERVIEW
# glmPath: Compute the regularization path for a penalized regression 
# nextNode_lvm: Find the next value of the regularization parameter

#' @title  Compute the regularization path for a penalized regression 
#' 
glmPath <- function(beta0, objective, hessian, gradient, 
                    indexPenalty, indexNuisance,
                    sd.X, base.lambda1, lambda2, group.lambda1,
                    control, iter.max = length(beta0)*3){
  
  p <- length(beta0)
  
  gradientPen <- function(beta, grad_Lv, ...){
    ifelse(beta == 0, -sign(grad_Lv), -sign(beta))  # because grad_Lv in lava is -grad(Lv)
  }
  
  #### initialisation
  ## the nuisance parameter is always fixed
  res <- initSigmaConstrain(beta0, constrain = control$constrain, indexNuisance = indexNuisance)
  beta0 <- res$start
  constrain <- res$constrain
  
  M.beta <- matrix(0, nrow = 1, ncol = p)
  colnames(M.beta) <- names(beta0)
  M.beta[1,] <- beta0
  
  V.lambda <- max( abs(-gradient(M.beta[1,]) * sd.X)[indexPenalty] )
  
  cv <- FALSE
  iter <- 1
  
  #### main loop
  while(iter < iter.max && cv == FALSE){
    if(control$regPath$trace){cat("*")}
    
    ## prediction step: next breakpoint assuming linear path and constant sigma
    resNode <- nextNode_lvm(hessian = hessian, gradient = gradient,  gradientPen = gradientPen,
                            beta = M.beta[iter,], lambda1 = V.lambda[iter] * base.lambda1, lambda2 = lambda2, 
                            indexPenalty = indexPenalty, indexNuisance = indexNuisance)
    
    newLambda <- V.lambda[iter] * (1 - resNode$gamma)
    
    if(newLambda < 0 || is.infinite(newLambda)){ ## cv, estimate for no penalization
      cv <- TRUE
      
      proxOperator <- function(x, step){
        control$proxOperator(x, step,
                             lambda1 = 0*base.lambda1, lambda2 = lambda2, test.penalty1 = group.lambda1, test.penalty2 = lambda2>0, expX = control$proxGrad$expX)
      }
      
      
      resNode$beta <- do.call("proxGrad",
                              list(start = M.beta[nrow(M.beta),], proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                                   constrain = constrain,
                                   step = control$proxGrad$step, BT.n = control$proxGrad$BT.n, BT.eta = control$proxGrad$BT.eta, trace = FALSE, force.descent = control$proxGrad$force.descent,
                                   iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, method = control$proxGrad$method))$par
      newLambda <- 0
      
    }else if(any(lambda2>0)){ ## correction step or update beta value with L2 penalization
      proxOperator <- function(x, step){  
        control$proxOperator(x, step,
                             lambda1 = newLambda*base.lambda1, lambda2 = lambda2, test.penalty1 = group.lambda1, test.penalty2 = lambda2>0, expX = control$proxGrad$expX)
      }
      
      resNode$beta <- do.call("proxGrad",
                              list(start = resNode$beta, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                                   constrain = constrain,
                                   step = control$proxGrad$step, BT.n = control$proxGrad$BT.n, BT.eta = control$proxGrad$BT.eta, trace = FALSE,  force.descent = control$proxGrad$force.descent,
                                   iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, method = control$proxGrad$method))$par
      
    } 
    
    ## update
    V.lambda <- c(V.lambda, newLambda)
    M.beta <- rbind(M.beta, resNode$beta)
    iter <- iter + 1
  }
  if(control$regPath$trace){cat("\n ")}
  
  #### export
  V.lambda <- unname(V.lambda)
  rownames(M.beta) <- NULL
  
  return(as.data.frame(cbind(lambda1.abs = V.lambda, 
                             lambda1 = NA, 
                             lambda2.abs = median(lambda2[indexPenalty]*sd.X[indexPenalty]^2), 
                             lambda2 = NA, 
                             M.beta)))
}

#' @title  Find the next value of the regularization parameter
nextNode_lvm <- function(hessian, gradient, gradientPen,
                         beta, lambda1, lambda2, indexPenalty, indexNuisance){
  
  ##
  hess_Lv <- -hessian(beta)  # because hess_Lv in lava is -hess(Lv)
  grad_Lv <- -attr(hess_Lv, "grad")  # because grad_Lv in lava is -grad(Lv) // gradient(beta)
  
  grad_Pen <- gradientPen(beta, grad_Lv = grad_Lv)
  
  set_A <- union(setdiff(which(beta!=0),indexNuisance),
                 setdiff(which(abs(grad_Lv)/lambda1 > 1-1e-04),indexNuisance)
  )
  
  grad_B.A <- try(solve(hess_Lv[set_A,set_A, drop = FALSE], grad_Pen[set_A, drop = FALSE] * lambda1[set_A, drop = FALSE]), silent = TRUE)
  if("try-error" %in% class(grad_B.A)){ # cannot invert because high dimensional 
    return(list(gamma = Inf, beta = rep(NA,length(beta))))
  }
  grad_Rho <- hess_Lv[indexPenalty,set_A, drop = FALSE] %*% grad_B.A / lambda1[indexPenalty, drop = FALSE]
  
  ## find next knot
  gamma1 <- -beta[set_A]/(grad_B.A*lambda1[set_A]) 
  gamma1 <- gamma1[set_A %in% indexPenalty] # remove active but non penalized coefficients
  gamma2Moins <- (1 - grad_Lv[indexPenalty, drop = FALSE]/lambda1[indexPenalty, drop = FALSE])/(1 + grad_Rho)
  gamma2Moins[indexPenalty %in% set_A] <- Inf
  gamma2Plus <- (1 + grad_Lv[indexPenalty, drop = FALSE]/lambda1[indexPenalty, drop = FALSE])/(1 - grad_Rho)
  gamma2Plus[indexPenalty %in% set_A] <- Inf
  
  # cat(paste(gamma1, collpase = " ")," | ",paste(gamma2Moins, collpase = " ")," | ",paste(gamma2Plus, collpase = " ")," | \n \n")
  
  test.min <- any(c(gamma1,gamma2Moins,gamma2Plus)>=0)
  if(test.min){
    gamma <- max(1e-04, min(c(gamma1,gamma2Moins,gamma2Plus)[c(gamma1,gamma2Moins,gamma2Plus)>=0]))
  }else{
    gamma <- Inf
  }
  
  ## update coef
  beta[set_A] <- beta[set_A] + gamma * grad_B.A
  
  return(list(gamma = gamma, beta = beta))
}


#' @title Prepare the data for the glmPath algorithm
prepareData_glmPath <- function(model, data, penalty, label = "_pen"){
  
  outcomes <- model$index$endogenous
  n.outcomes <- length(outcomes)
  exogeneous <- model$exogenous
  names.coef <- coef(model)
  n.coef <- length(names.coef)
  
  #### rename penalized variables 
  for(iter_Y in 1:n.outcomes ){
    Y_tempo <- outcomes[iter_Y]
    
    index_coef <- grep(Y_tempo, x = penalty$names.penaltyCoef, fixed = TRUE) # any penalize coefficient related to outcome Y_tempo
    
    if( length(index_coef)>0 ){
      ls.tempo <- lapply(exogeneous, grep, penalty$names.penaltyCoef[index_coef], fixed = TRUE) # which exogeneous varaible is related to outcome Y_tempo
      names.varData <- exogeneous[unlist(lapply(ls.tempo, function(t){length(t)>0}))] 
      data <- cbind(data,
                    setNames(data[, names.varData, drop = FALSE], paste0(names.varData, label, Y_tempo))
      )
      
      # remove link corresponding to a penalized coefficient
      for(iter_link in unlist(ls.tempo)){
        cancel(model) <- as.formula(penalty$names.penaltyCoef[index_coef][iter_link])
      }
      
      # add new link corresponding to the penalized coeffcient but renaæed
      for(iter_link in unlist(ls.tempo)){
        regression(model) <- as.formula(paste0(penalty$names.penaltyCoef[index_coef][iter_link], label, Y_tempo))
      }
      
      # update penalty
      penalty$names.penaltyCoef[index_coef][unlist(ls.tempo)] <- paste0(penalty$names.penaltyCoef[index_coef][unlist(ls.tempo)], label, Y_tempo)
    }
    
  }
  
  #### remove useless variables
  varCoef <- unlist(strsplit(coef(model), split = "~", fixed = TRUE))
  varCoef <- unlist(strsplit(varCoef, split = ",", fixed = TRUE))
  
  indexRM <- which(model$exogenous %in% unique(varCoef) == FALSE)
  if(length(indexRM)>0){
    data <- data[,-which(names(data) %in%  model$exogenous[indexRM])]
    kill(model) <- model$exogenous[indexRM]
  }
  
  #### orthogonalize the dataset
  names.coef <- coef(model)
  scaleLambda1 <- setNames(rep(0, n.coef),names.coef)
  scaleLambda2 <- setNames(rep(0, n.coef),names.coef)
  mu.X <- setNames(rep(0, n.coef),names.coef)
  sd.X <- setNames(rep(1, n.coef),names.coef)
  orthogonalizer <- setNames(vector("list", n.outcomes), outcomes)
  
  for(iter_Y in 1:n.outcomes ){
    Y_tempo <- outcomes[iter_Y]
    
    index_coef <- grep(Y_tempo, x = penalty$names.penaltyCoef, fixed = TRUE) # any penalize coefficient related to outcome Y_tempo
    if( length(index_coef)>0 ){
      
      # orthogonalize data relative to non penalize coefficients
      resOrtho <- orthoData_glmPath(model, name.Y = Y_tempo,
                                    allCoef = names.coef[grep(Y_tempo, names.coef, fixed = TRUE)], 
                                    penaltyCoef = penalty$names.penaltyCoef[index_coef], 
                                    data = data)
      
      # update results
      data[, colnames(resOrtho$orthogonalizer)] <- resOrtho$data[, colnames(resOrtho$orthogonalizer),drop = FALSE]
      scaleLambda1[names(resOrtho$lambda1)] <- resOrtho$lambda1
      scaleLambda2[names(resOrtho$lambda2)] <- resOrtho$lambda2
      mu.X[names(resOrtho$lambda1)] <- resOrtho$mu.X
      sd.X[names(resOrtho$lambda1)] <- resOrtho$sd.X
      orthogonalizer[[iter_Y]] <- resOrtho$orthogonalizer
    }
    
  }
  
  #### export
  penalty$sd.X <- sd.X
  penalty$scaleLambda1 <- scaleLambda1
  penalty$scaleLambda2 <- scaleLambda2
  
  return(list(model = model,
              penalty = penalty,
              data = data,
              orthogonalizer = orthogonalizer,
              mu.X = mu.X))
  
}

#' @title Orthogonalize the non-penalized variables relatively to the penalized variables for a regression model
orthoData_glmPath <- function(model, name.Y, allCoef, penaltyCoef, data){
  
  ## function
  extractVar <- function(names){
    names.formula <- grep("~", names ,fixed = TRUE)  
    if(length(names.formula)>0){
      names[names.formula] <- unlist(lapply(strsplit(names[names.formula], split = "~", fixed = TRUE),"[",2))
    }
    
    return(names)
  }
  
  ## preparation
  n <- nrow(data)
  n.coef <- length(allCoef)
  names.interceptCoef <- intersect(allCoef,coef(model)[model$index$parBelongsTo$mean])
  n.interceptCoef <- length(names.interceptCoef)
  names.covCoef <- intersect(allCoef,coef(model)[model$index$parBelongsTo$cov])
  if(length(model$latent)>0){
    names.latentCoef <- allCoef[sapply(names(model$latent), grep, x = allCoef, fixed = TRUE)]
  }else{
    names.latentCoef <- NULL
  }
  
  var.penalized <-  extractVar(penaltyCoef)
  var.unpenalized <- setdiff(extractVar(setdiff(allCoef,c(names.covCoef,names.latentCoef))), 
                             c(var.penalized))
  
  ## rebuild data
  X_tempo <- data[,setdiff(c(var.penalized,var.unpenalized),names.interceptCoef), drop = FALSE]
  
  if(n.interceptCoef>0){
    X_tempo <- cbind(matrix(1, nrow = n, ncol = n.interceptCoef),
                     X_tempo)
  }
  names(X_tempo)[1:n.interceptCoef] <- names.interceptCoef
  
  ## distinguish penalized from non penalized
  penalized <-  as.matrix(X_tempo[,var.penalized, drop = FALSE])
  unpenalized <-  as.matrix(X_tempo[,var.unpenalized, drop = FALSE])
  
  ## orthogonlize
  orthogonalizer <- solve(crossprod(unpenalized), crossprod(unpenalized, penalized))
  penalized <- penalized - unpenalized %*% orthogonalizer
  
  ## starting coefficients
  mu.X <- setNames(rep(0,n.coef), allCoef)
  lm.fitted <- lm.fit(y = data[[name.Y]], x = unpenalized)
  mu.X[names(mu.X) %in% penaltyCoef == FALSE] <- c(coef(lm.fitted), var(lm.fitted$residuals)) ### issue with the latent variable here!!
  
  ## scale
  sd.X <- setNames(rep(1,n.coef), allCoef)
  
  varNI.penalized <- setdiff(var.penalized, names.interceptCoef)
  index.penalized <- setdiff(which(names(mu.X) %in% penaltyCoef),
                             which(names(mu.X) %in% names.interceptCoef) )
  
  if(length(varNI.penalized)>0){
    sd.X[index.penalized] <- sqrt(apply(penalized[,varNI.penalized, drop = FALSE], 2, var)*(n-1)/n)
    penalized <- sweep(penalized[,varNI.penalized, drop = FALSE], MARGIN = 2, FUN = "/", STATS = sd.X[index.penalized])
  }
  
  varNI.unpenalized <- setdiff(var.unpenalized, names.interceptCoef)
  index.unpenalized <- setdiff(which(names(mu.X) %in% penaltyCoef == FALSE), 
                               which(names(mu.X) %in% c(names.interceptCoef, names.covCoef))
  )
  
  if(length(varNI.unpenalized)>0){
    sd.X[index.unpenalized] <- sqrt(apply(unpenalized[,varNI.unpenalized, drop = FALSE], 2, var)*(n-1)/n)
    unpenalized <- sweep(unpenalized[,varNI.unpenalized, drop = FALSE], MARGIN = 2, FUN = "/", STATS = sd.X[index.unpenalized])
  }
  mu.X <- mu.X * sd.X
  
  ## update initial dataset
  data[,var.penalized] <- penalized
  if(length(setdiff(var.unpenalized, names.interceptCoef))>0){
    data[,setdiff(var.unpenalized, names.interceptCoef)] <- unpenalized
  }
  
  ## lambda
  lambda1 <- setNames(rep(0,n.coef), allCoef)
  lambda1[which(names(mu.X) %in% penaltyCoef)] <- 1/sd.X[which(names(mu.X) %in% penaltyCoef)]
  lambda2 <- setNames(rep(0,n.coef), allCoef)
  lambda2[which(names(mu.X) %in% penaltyCoef)] <- 1/(sd.X[which(names(mu.X) %in% penaltyCoef)]^2)
  
  ## export
  return(list(data = data,
              orthogonalizer = orthogonalizer,
              sd.X = sd.X,
              mu.X = mu.X,
              lambda1 = lambda1,
              lambda2 = lambda2))
}

#' @title Cancel the effect of the orthogonalization on the estimated parameters
rescaleCoef_glmPath <- function(Mres, penalty, orthogonalizer){
  
  #### rescale
  names.rescale <- names(penalty$sd.X)
  Mres[,names.rescale] <- sweep(Mres[,names.rescale], MARGIN = 2, FUN = "/", STATS = penalty$sd.X)
  
  #### de-orthogonalize
  name.Y <- names(orthogonalizer)
  covCoef <- penalty$names.varCoef
  n.Y <- length(name.Y)
  
  for(iter_Y in 1:n.Y){
    Y_tempo <- name.Y[iter_Y]
    names.allCoef <- colnames(Mres)[grep(Y_tempo, x = colnames(Mres), fixed = TRUE)] # coefficients related to outcome Y_tempo
    names.penalizedCoef <- intersect(penalty$names.penaltyCoef, names.allCoef) # penalized coefficient related to outcome Y_tempo
    names.unpenalizedCoef <- setdiff(names.allCoef, c(names.penalizedCoef,covCoef)) # penalized coefficient related to outcome Y_tempo
    
    Mres[,names.unpenalizedCoef] <- Mres[,names.unpenalizedCoef] - Mres[,names.penalizedCoef] %*% t(orthogonalizer[[iter_Y]])
  }
  
  return(as.data.frame(Mres))
}