`print.plvm` <- function(x, ...) {
  
  ## normal display
  out <- capture.output(lava:::print.lvm(x))
  if(!is.null(x$penaltyNuclear$name.Y)){
    charY <- paste0(x$penaltyNuclear$name.Y," ~ ")
    indexEq <- grep(charY,out)
    out[indexEq] <- gsub(charY, replacement = paste0(charY,LCSseq(x$penaltyNuclear$name.X),"(image)+"),x = out[indexEq])
  }
  sapply(out, function(o){cat(o,"\n")})
  
  ## additional display - lasso
  if(!is.null(x$penalty$names.penaltyCoef)){
    if(x$penalty$lambda1>0 && x$penalty$lambda2>0){
      penaltyType <- "Elastic net"
    }else if(x$penalty$lambda1>0){
      penaltyType <- "Lasso"
    }else if(x$penalty$lambda1>0){
      penaltyType <- "Ridge"
    }else{
      penaltyType <- "None"
    }
    
    if(all(x$penalty$group.penaltyCoef<1)){
      cat("Penalty: ", penaltyType,"\n",
          "On     : ", paste(x$penalty$names.penaltyCoef, collapse = " "),"\n")
    }else{
      
      test.lasso <- (x$penalty$group.penaltyCoef<1)*(x$penalty$group.penaltyCoef>0)
      if(any(test.lasso==1)){
        cat("Penalty: ", penaltyType,"\n",
            "On     : ", paste(x$penalty$names.penaltyCoef[test.lasso==1], collapse = " "),"\n")
      }
      
      ls.penalty <- tapply(x$penalty$names.penaltyCoef[test.lasso!=1], x$penalty$group.penaltyCoef[test.lasso!=1],list)
      cat("Penalty: Grouped lasso \n")
      lapply(ls.penalty, function(x){cat("On     :",paste(x, collapse = " "),"\n")})
    }
    cat("\n")
  }
  if(!is.null(x$penaltyNuclear$name.Y)){
    penaltyType <- "Nuclear norm"
    cat("Penalty: ", penaltyType,"\n",
        "on     : ", LCSseq(x$penaltyNuclear$name.X)," (outcome: ",x$penaltyNuclear$name.Y,")\n",sep = "")
  }
  ## export
  invisible(x)
}





##' @export
`print.plvmfit` <- function(x,level=2,labels=FALSE, 
                            getCoef = "penalized", getLambda = "abs", rm.duplicated = TRUE,
                            ...) {
  if(!is.null(x$penalty$lambda1.best)){
    lava:::print.lvmfit(x)
    cat("\n Model selected using ",attr(x$penalty$performance,"criterion"),"criterion \n")
    cat("   range of lambda1: ",paste(range(x$penalty$lambda1), collapse = " "),"\n")
    cat("   best lambda1    : ",x$penalty$lambda1.best,"\n")
    
  }else if(x$penalty$regularizationPath == 0){
    
    Mtempo <- CoefMat(x,labels=labels,level=level,...) 
    ncol.M <- ncol(Mtempo)
    if(x$penalty$lambda1>0 || x$penalty$lambda2>0){
      Mtempo <- rbind(Mtempo, "Penalization:" = rep("", ncol.M))
    }
    if(x$penalty$lambda1>0){
      Mtempo <- rbind(Mtempo, "   L1 lambda (abs)" = c(x$penalty$lambda1.abs, rep("",ncol.M-1)))
    }
    if(x$penalty$lambda2>0){
      Mtempo <- rbind(Mtempo, "   L2 lambda (abs)" = c(x$penalty$lambda2.abs, rep("",ncol.M-1)))
    }
    
    print(Mtempo,quote=FALSE,right=TRUE)
    minSV <- attr(vcov(x),"minSV")
    if (!is.null(minSV) && minSV<1e-12) {
      warning("Small singular value: ", format(minSV))
    }
    pseudo <- attr(vcov(x),"pseudo")
    if (!is.null(pseudo) && pseudo) warning("Singular covariance matrix. Pseudo-inverse used.")
    
    
  }else {
    cat("Regularization path: \n")
    printPath <- getPath(x, rm.duplicated = rm.duplicated, getCoef = getCoef, getLambda = getLambda)
    print(printPath)
    diffRow <- nrow(getPath(x)) - nrow(printPath)
    if(diffRow>0){cat("[ omitted ",diffRow," rows ] \n",sep = "")}
    cat("estimated using ")
    switch(x$penalty$regularizationPath,
           "1" = cat("glmPath algorithm \n"),
           "2" = cat("EPSODE algorithm \n"))
  }


invisible(x)
}

##' @title get the common substring sequence in a vector of strings
LCSseq <- function(x){
  affixe <- strsplit(x[[1]], split = "")[[1]]
  
  for(iterX in 2:length(x)){
    affixe <- qualV::LCS(affixe, strsplit(x[[iterX]], split = "")[[1]])$LCS
  }
  
  return(affixe)
}