`extendModel` <-
  function(x,...) UseMethod("extendModel")

###

extendModel.lvm <- function(x, data, type, alpha = 0.05, covariance = TRUE, warn = TRUE, trace = TRUE, ...){
  
  match.arg(type, choices = c("all","modelsearch", "modelsearchLR"))
  
  if(type != "all"){
   
    ## handling factor variables
    test.factor <- unlist(lapply(data, function(x){is.factor(x) + is.character(x) > 0}))
    varFactor <- names(test.factor)[test.factor]
    ls.levels <- lapply(varFactor, function(x){levels(data[[x]])})
    names(ls.levels) <- varFactor
    
    lvmfit <- estimate(x, data = data, ...)
    cv <- FALSE
    
    while(cv == FALSE){ 
      if(trace){cat("*")}
      resSearch <- do.call(type, args = list(lvmfit, silent = TRUE))
     
      if(tail(p.adjust(resSearch$test[,"P-value"], method = "holm"), 1) < alpha){
        var1 <- tail(resSearch$var,1)[[1]][,1]
        var2 <- tail(resSearch$var,1)[[1]][,2]
        
        if(var1 %in% vars(x) == FALSE){
          var1 <- renameFactor(var1, ls.levels = ls.levels)
          if(warn && length(ls.levels[[var1]])>2){
            warning("extendModel.lvm: one of the levels of a factor variable reach the significance level \n",
                    "a link with the whole variable is added in the model \n")
          }
        }
        if(var2 %in% vars(x) == FALSE){
          var2 <- renameFactor(var2, ls.levels = ls.levels)
          if(warn && length(ls.levels[[var2]])>2){
            warning("extendModel.lvm: one of the levels of a factor variable reach the significance level \n",
                    "a link with the whole variable is added in the model \n")
          }
        }
        
        x <- addLink(x, var1, var2,
                     covariance = covariance)
        
        lvmfit <- estimate(x, data = data, ...)
      }else{
        cv <- TRUE
      }
      
      if(lvmfit$opt$iterations == lvmfit$control$iter.max){
        return(NA)
      }
    }
    if(trace){cat("\n")}
    
    return(lvmfit)
    
  }else{ # all
    
    newlinks <- findNewLink(x, rm.exoexo = TRUE)
    for(iterLink in 1:nrow(newlinks)){
     x <- addLink(x, newlinks[iterLink,1], newlinks[iterLink,2], covariance = covariance,
                   silent = TRUE)
    }
    
    return(x)
  }
  

}


`modelsearchLR` <- function(object, ...) UseMethod("modelsearchLR")

modelsearchLR.lvmfit <- function (object, silent = FALSE, ...){

  #### newlinks 
  restricted <- findNewLink(object$model, rm.exoexo = FALSE, output = "names")
  seq_i <- seq_len(NROW(restricted))
  
  #### initialisation
  M.test <- cbind("Test Statistic" = rep(NA,length(seq_i)),
                  "P-value" = rep(NA,length(seq_i))
  )
  ls.var <- list()
  
  if(silent == FALSE){pb <- utils::txtProgressBar(max = tail(seq_i,1), style = 3) }
  
  for (iterI in seq_i) {
    
    newmodel <- addLink(object$model, var1 = restricted[iterI,1], var2 = restricted[iterI,2],
                        covariance = all(restricted[iterI,1:2] %in% endogenous(object$model)))
    newcontrol <- object$control
    newcontrol$start <- coef(object)
    newcontrol$trace <- FALSE
    
    newfit <- tryCatch(estimate(newmodel, data = object$data$model.frame, control = newcontrol, 
                                missing = "lvm.missing" %in% class(object), ...),
                       error = function(x){NA},
                       finally = function(x){x})
    
    if("lvmfit" %in% class(newfit)){
      compareT <- compare(object,newfit)
      M.test[iterI,] <- c(compareT$statistic[[1]], compareT$p.value[[1]])
    }
    
    ls.var[[iterI]] <- matrix(c(restricted[iterI,1], restricted[iterI,2]), nrow = 1)
    if(silent == FALSE){ utils::setTxtProgressBar(pb, value = iterI) }
    
  }
  if(silent == FALSE){  close(pb) }
  
  #### reorder
  index.order <- order(M.test[,"P-value"], decreasing = TRUE)
  ls.var <- ls.var[index.order]
  M.test <- M.test[index.order,]
  
  #### export
  M.res <- cbind(apply(M.test, 2, signif, digit = 4), unlist(lapply(ls.var, paste, collapse = ",")))
  colnames(M.res) <- c("Score:", "S P(S>s)", "Index")
  rownames(M.res) <- rep("",nrow(M.res))
  
  output <- list(res = M.res,
                 test = M.test,
                 var = ls.var)
  class(output) <- "modelsearch"
  return(output)
}


renameFactor <- function(var, ls.levels, data, sep = ""){
  
  if(!missing(data)){
    test.factor <- unlist(lapply(data, function(x){is.factor(x) + is.character(x) > 0}))
    allvars <- names(test.factor)[test.factor]
    ls.levels <- lapply(allvars, function(x){levels(data[[x]])})
  }else {
    allvars <- names(ls.levels)
  }
 
  test <- unlist(lapply(1:length(allvars), function(x){any(paste(allvars[x], ls.levels[[x]], sep = sep) == var)}))
  
  return(allvars[unlist(test)])
}



