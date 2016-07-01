`extendModel` <-
  function(x,...) UseMethod("extendModel")

`extendNames` <-
  function(x,...) UseMethod("extendNames")

`coefVar` <-
  function(x,...) UseMethod("coefVar")

`coefCov` <-
  function(x,...) UseMethod("coefCov")

`addLink` <-
  function(x,...) UseMethod("addLink")

`setLink` <-
  function(x,...) UseMethod("setLink")

`rmLink` <-
  function(x,...) UseMethod("rmLink")

`rmLink` <-
  function(x,...) UseMethod("rmLink")

###

extendModel.lvm <- function(x, data, type, alpha = 0.05, covariance = TRUE, warn = TRUE, trace = TRUE, ...){
  
  match.arg(type, choices = c("all","modelsearch"))
  
  if(type == "modelsearch"){
   
    ## handling factor variables
    test.factor <- unlist(lapply(data, function(x){is.factor(x) + is.character(x) > 0}))
    varFactor <- names(test.factor)[test.factor]
    ls.levels <- lapply(varFactor, function(x){levels(data[[x]])})
    names(ls.levels) <- varFactor
    
    lvmfit <- estimate(x, data = data, ...)
    cv <- FALSE
    
    while(cv == FALSE){ 
      if(trace){cat("*")}
      resSearch <- modelsearch(lvmfit, silent = TRUE)

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
    
    newlinks <- extendNames(x, rm.exoexo = TRUE)
    for(iterLink in 1:nrow(newlinks)){
      x <- addLink(x, newlinks[iterLink,1], newlinks[iterLink,2], covariance = covariance,
                   silent = TRUE)
    }
    
    return(x)
  }
  

}

extendNames.lvm <- function(x, rm.exoexo){
  AP <- with(index(x), A + t(A) + P)
  
  restricted <- c()
  for (i in seq_len(ncol(AP) - 1)){
    for (j in seq(i + 1, nrow(AP))){
      test.exo <- (rownames(AP)[i] %in% exogenous(x)) + (colnames(AP)[j] %in% exogenous(x))
      if (AP[j, i] == 0 && (rm.exoexo == FALSE || test.exo!=2)){
          restricted <- rbind(restricted, c(i, j)) 
      }
    }
  }
 
  if (is.null(restricted)){
    return(NULL)
  }
  
  names.restricted <- cbind(rownames(AP)[restricted[,1]],
                            colnames(AP)[restricted[,2]])
  
  return(names.restricted)
  
}


addLink.lvm <- function(x, var1, var2 = NA, covariance, silent = FALSE){
 
  res <- initVar_link(var1, var2)
  var1 <- res$var1
  var2 <- res$var2
  
  if(var1 %in% vars(x) == FALSE){
    if(silent == FALSE){
      warning("addLink.lvm: var1 does not match any variable in x, no link is added \n",
              "var1: ",var1,"\n")
    }
  }
  
  ####
  if(is.na(var2)){
    
    intercept(x) <- as.formula(paste0("~", var1))
    
  }else{
    
    if(var1 == var2){
      if(silent == FALSE){
        warning("addLink.lvm: var1 equals var2, no link is added \n",
                "var1/2: ",var1,"\n")
      }
    }
    
    
    if(var2 %in% vars(x) == FALSE){
      if(silent == FALSE){
        warning("addLink.lvm: var2 does not match any variable in x, no link is added \n",
                "var2: ",var2,"\n")
      }
    }
    test.1 <- var1 %in% exogenous(x)
    test.2 <- var2 %in% exogenous(x)
    
    
    if(test.1 && test.2){
      if(silent == FALSE){
        warning("addLink.lvm: both variable are exogenous, no link is added \n",
                "var1: ",var1,"\n",
                "var2: ",var2,"\n")
      }
    }else if(test.1){
      regression(x) <- as.formula(paste(var2, var1,  sep = "~"))
    }else if(test.2){
      regression(x) <- as.formula(paste(var1, var2, sep = "~"))
    }else if(covariance){
      covariance(x) <- as.formula(paste(var1, var2, sep = "~"))  
    }
  }
  
  return(x)
}

setLink.lvm <- function(x, var1, var2 = NA, value, silent = FALSE){
  
  res <- initVar_link(var1, var2)
  var1 <- res$var1
  var2 <- res$var2
  
  #### set the link
  if(is.na(var2)){
    intercept(x, as.formula(paste0("~",var1))) <- value
  }else if(paste(var1, var2, sep = "~") %in% coef(x)){
    regression(x, as.formula(paste(var1,var2, sep = "~"))) <- value
  }else if(paste(var1,var2, sep = ",") %in% coef(x)){
    covariance(x, as.formula(paste(var1,var2, sep = "~"))) <- value
  }else if(paste(var2,var1, sep = ",") %in% coef(x)){
    covariance(x, as.formula(paste(var1,var2, sep = "~"))) <- value
  }else{
    if(silent == FALSE){
      warning("setLink.lvm: no link was found from var1 to var2, no link is set \n",
              "var1: ",var1,"\n",
              "var2: ",var2,"\n")
    }
  }
  
  return(x)
}

rmLink.lvm <- function(x, var1, var2 = NA, silent = FALSE){
  
  res <- initVar_link(var1, var2)
  var1 <- res$var1
  var2 <- res$var2
  
  #### remove the link
  if(is.na(var2)){
    cancel(x) <- as.formula(paste0("~",var1))
  }else if(paste(var1, var2, sep = "~") %in% coef(x)){
    cancel(x) <- as.formula(paste(var1,var2, sep = "~"))
  }else if(paste(var1,var2, sep = ",") %in% coef(x)){
    cancel(x) <-  as.formula(paste(var1,var2, sep = "~"))
  }else if(paste(var2,var1, sep = ",") %in% coef(x)){
    cancel(x) <-  as.formula(paste(var2,var1, sep = "~"))
  }else{
    if(silent == FALSE){
      warning("addLink.lvm: no link was found from var1 to var2, no link is removed \n",
              "var1: ",var1,"\n",
              "var2: ",var2,"\n")
    }
  }
  
  
  #### if unused variable remove it from the model
  if(length(grep(paste0("~",var1,"$|^",var1,"~|^",var1,"$"), x = coef(x), fixed = FALSE))==0){
      kill(x) <- as.formula(paste0(var1,"~1"))
  }
  if(!is.na(var2) && length(grep(paste0("~",var2,"$|^",var2,"~|^",var2,"$"), x = coef(x), fixed = FALSE))==0){
    kill(x) <- as.formula(paste0(var2,"~1"))
  }
  
  return(x)
}

initVar_link <- function(var1, var2){
  
  if(is.na(var2) && is.character(var1)){
    if(grepl(",",var1)==TRUE){var1 <- gsub(",","~", x = var1)}
    if(grepl("~",var1)==TRUE){var1 <- as.formula(var1)}
  }
  
  if(class(var1) == "formula"){
    var2 <- all.vars(var1)[2]
    var1 <- all.vars(var1)[1]
  }
  
  return(list(var1 = var1,
              var2 = var2))
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

coefVar.lvm <- function(x, index = FALSE){
  names.var <- paste(x$index$vars, x$index$vars, sep = ",")
  return(grep(paste(names.var, collapse = "|"), coef(x), value = index))
}

coefCov.lvm <- function(x, index = FALSE){
  names.cov <- setdiff(coef(x)[x$index$parBelongsTo$cov],
                       paste(x$index$vars, x$index$vars, sep = ",")
  )
  return(grep(paste(names.cov, collapse = "|"), coef(x), value = index))
}

