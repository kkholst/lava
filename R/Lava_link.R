`findNewLink` <-
  function(x,...) UseMethod("findNewLink")

`addLink` <-
  function(x,...) UseMethod("addLink")

`setLink` <-
  function(x,...) UseMethod("setLink")

`rmLink` <-
  function(x,...) UseMethod("rmLink")

#' @title Find the new possible links between variables (copied from lava::modelsearch)
findNewLink.lvm <- function(x, rm.exoexo, output = "names"){
  
  match.arg(output, choices = c("names","index"))
  
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
  
  if(output == "names"){
    names.restricted <- cbind(rownames(AP)[restricted[,1]],
                              colnames(AP)[restricted[,2]])
    
    return(names.restricted)
  }else{
    return(restricted)
  }
  
}

#' @title Add a new link between two variables in a lvm
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
    if(var1 %in% endogenous(x) && var1 %in% endogenous(x)){
      covariance <- TRUE
    }
    
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
    }else {
      if(var1 %in% endogenous(x) && var2 %in% latent(x)){
        regression(x) <- as.formula(paste(var1, var2, sep = "~"))  
      }else if(var2 %in% endogenous(x) && var1 %in% latent(x)){
        regression(x) <- as.formula(paste(var2, var1, sep = "~"))  
      }else{
        stop("addLink.lvm: unknow configuration \n")
      }
      
    }
  }
  
  return(x)
}

#' @title Affect a given value to a link between two variables in a lvm 
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

#' @title Remove a new link between two variables in a lvm model
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

#' @title Convert var1 and var2 from formula to character
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
