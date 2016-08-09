calcLambda <- function(object, model, seq_lambda, data.fit, data.test, fit = "BIC", trace = TRUE, ...){
  # if no test set then CV
  
 if(fit %in% c("AIC","BIC","P_error") == FALSE){
   stop("fit must be in AIC BIC P_error \n")
 }
  if(all(c("plvmfit", "lvm") %in% class(object) == FALSE)){
    stop("model must be either a plvmfit or a lvm object \n")
  }
  
  if("plvmfit" %in% class(object)){
    regPath <- getPath(object, getLambda = "lambda1.abs", getCoef = "coef0")
    seq_lambda <- unlist(lapply(regPath, function(x){attr(x,"lambda1.abs")}))
    seq_row <- unlist(lapply(regPath, function(x){attr(x,"row")}))
    seq_coef <- lapply(regPath, function(x){as.character(x)})
  }
  n.lambda <- length(seq_lambda)
  res.cv <- numeric(n.lambda)
  n.endogeneous <- length(endogenous(model))
  
  best.res <- Inf
  best.lambda <- NA
  best.subset <- NULL
  coef.penalty <- model$penalty$name.coef
  
  if(trace){pb <- txtProgressBar(min = 0, max = n.lambda, style = 3)}
  for(iterLambda in 1:n.lambda){
    if(trace){setTxtProgressBar(pb, iterLambda)}
    
    
    #### define the variables to include in the model 
    if("plvmfit" %in% class(object) == FALSE){
      fitTempo <- estimate(model, data = data.fit, lambda1 = seq_lambda[iterLambda], ...)
      coef0_lambda <- coef0(fitTempo, tol = 1e-6, penalized = TRUE, value = FALSE)
    }else{
      coef0_lambda <- seq_coef[[iterLambda]]
    }
    
    #### form the reduced model
    model2 <- model
    if(length(coef0_lambda)>0){
      ToKill <- coef0_lambda
      for(iterKill in ToKill){
        model2 <- rmLink.lvm(model2, iterKill)  
      }
    }
    
    #### fit the reduced model
    fitTempo2 <- lava:::estimate.lvm(model2, data = data.fit, control = list(constrain = TRUE))
   
    #### gof criteria
    if(fit %in% c("AIC","BIC")){
      res.cv[iterLambda] <- gof(fitTempo2)[[fit]]
    }else{
      predY <- as.data.frame(predict(fitTempo2,
                                     data = data.test[,exogenous(fitTempo2), drop = FALSE]))
      
      res.cv[iterLambda] <- 1/(n.endogeneous*NROW(data.test)) * sum((data.test[,endogenous(fitTempo2)] - predY)^2)
    }
    
    #### storage
    if(best.res>res.cv[iterLambda] ){
      best.res <- res.cv[iterLambda]
      best.lambda <- seq_lambda[iterLambda]
      best.subset <- names(coef(fitTempo2))
      best.lvm <- fitTempo2
      if("plvmfit" %in% class(object)){attr(best.lambda,"row") <- seq_row[iterLambda]}
    }
  }
  if(trace){close(pb)}
  if("plvmfit" %in% class(object)){
    best.lvm$penalty <- object$penalty 
    best.lvm$penalty$lambda1 <- seq_lambda 
    attr(res.cv, "criterion") <- fit
    best.lvm$penalty$performance <- res.cv
    best.lvm$penalty$lambda1.best <- best.lambda
    best.lvm$regularizationPath <- object$regularizationPath
    class(best.lvm) <- append("plvmfit", class(best.lvm))
    return(best.lvm)
  }else{
    return(list(resCV = res.cv,
                lambda = best.lambda,
                subset = best.subset))
  }
}
