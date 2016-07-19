calcLambda <- function(model, seq_lambda, data.fit, data.test, trace = TRUE, ...){
  # if no test set then CV
  
  if("plvmfit" %in% class(model)){
    penPath <- penPath(model)
    model <- model$model
    seq_lambda <- penPath$lambda1.abs
  }else{
    penPath <- NULL
  }
  
  n.lambda <- length(seq_lambda)
  res.cv <- numeric(n.lambda)
  n.endogeneous <- length(endogenous(model))
  n.obs <- nrow(data.test)
  
  best.res <- Inf
  best.lambda <- NA
  best.subset <- NULL
  coef.penalty <- model$penalty$names.penaltyCoef
  
  if(trace){pb <- txtProgressBar(min = 0, max = n.lambda, style = 3)}
  for(iterLambda in 1:n.lambda){
    if(trace){setTxtProgressBar(pb, iterLambda)}
    
    if(is.null(penPath)){
      fitTempo <- estimate(model, data = data.fit, lambda1 = seq_lambda[iterLambda], ...)
      coef0_lambda <- coef0(fitTempo, tol = 1e-6, penalized = TRUE)
    }else{
      # fitTempo <- estimate(model, data = data.fit, lambda1 = 10,  control = list(trace = TRUE, constrain = TRUE))
      coef_lambda <- penPath[iterLambda,intersect(names(penPath),model$penalty$names.penaltyCoef)]
      coef0_lambda <- setNames(coef_lambda[coef_lambda == 0],colnames(coef_lambda == 0)[coef_lambda == 0])
    }
    
    model2 <- model
    if(length(coef0_lambda)>0){
      ToKill <- names(coef0_lambda)
      for(iterKill in ToKill){
        model2 <- rmLink.lvm(model2, as.formula(iterKill))  
      }
    }
    
    fitTempo2 <- lava:::estimate.lvm(model2, data = data.fit, control = list(constrain = TRUE))
    predY <- as.data.frame(predict(fitTempo2,
                                   data = data.test[,exogenous(fitTempo2), drop = FALSE]))
    
    res.cv[iterLambda] <- 1/(n.endogeneous*n.obs) * sum((data.test[,endogenous(fitTempo2)] - predY)^2)
    
    if(best.res>res.cv[iterLambda] ){
      best.res <- res.cv[iterLambda]
      best.lambda <- seq_lambda[iterLambda]
      best.subset <- names(coef(fitTempo2))
    }
  }
  if(trace){close(pb)}
  
  
  sd.index <- tail(which(res.cv <= min(res.cv)+sd(res.cv)),1)
  sd.lambda <- seq_lambda[sd.index]
  sd.subset <- coefN0(estimate(model, data = data.fit, lambda1 = sd.lambda))
  
  return(list(resCV = res.cv,
              lambda = best.lambda,
              subset = best.subset,
              lambda1sd = sd.lambda,
              subset1sd = sd.subset))
}
