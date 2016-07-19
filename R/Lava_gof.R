#' @title Display the correlation between the residuals

`corDiag` <- function(object, ...) UseMethod("corDiag")

corDiag.lvmfit <- function(object, endogenous = NULL, adjust = "none"){
  
  if(is.null(endogenous)){endogenous <- endogenous(object)}
  butils:::validCharacter(endogenous, validValues = endogenous(object), validLength = NULL,
                          refuse.duplicates = TRUE, method = "lmDiag.lvmfit")
  
  fitted <- predict(object)
  residual <- matrix(NA, ncol = ncol(fitted), nrow = nrow(fitted))
  colnames(residual) <- endogenous
  
  for(iterEndo in endogenous){
    if(iterEndo %in% colnames(fitted) == FALSE){
      possibleVar <- paste0(iterEndo,"~",latent(object$model))
      test <- lapply(possibleVar, function(x){grep(x,coef(object$model))})
      index.latent <- which(unlist(lapply(test, function(x){any(!is.na(x))})))
      if(length(index.latent)>1){
        stop("lmDiag.lvmfit: cannot know how to reconstruct prediction from multiple latent variables \n")
      }
      residual[,iterEndo] <- predict(object)[,latent(object$model)[index.latent]] - object$data$model.frame[[iterEndo]]
    }else{
      residual[,iterEndo] <- predict(object)[,iterEndo] - object$data$model.frame[[iterEndo]]
    }
  }
  
  return(psych:::corr.test(residual, adjust = adjust))
  
}

#' @title Diagnostic for the model relating the LV and the endogenous variables
`lmDiag` <- function(object, ...) UseMethod("lmDiag")

lmDiag.lvmfit <- function(object, endogenous = NULL, which = c(0,2:3,5), digit = 3,
                          mfrow = c(2,2), mar = c(3,3,3,3), mgp = c(2,0.75,0), oma=c(0,0,1,0)){
  
  
  if(is.null(endogenous)){endogenous <- endogenous(object)}
  butils:::validCharacter(endogenous, validValues = endogenous(object), validLength = NULL,
                          refuse.duplicates = TRUE, method = "lmDiag.lvmfit")
  
  fitted <- predict(object)
  ls.lm <- list()
  
  
  
  
  for(iterEndo in endogenous){
    if(iterEndo %in% colnames(fitted) == FALSE){
      possibleVar <- paste0(iterEndo,"~",latent(object$model))
      test <- lapply(possibleVar, function(x){grep(x,coef(object$model))})
      index.latent <- which(unlist(lapply(test, function(x){any(!is.na(x))})))
      if(length(index.latent)>1){
        stop("lmDiag.lvmfit: cannot know how to reconstruct prediction from multiple latent variables \n")
      }
      df <- data.frame(observed = object$data$model.frame[[iterEndo]],
                       predicted = predict(object)[,latent(object$model)[index.latent]])
    }else{
      df <- data.frame(observed = object$data$model.frame[[iterEndo]],
                       predicted = predict(object)[,iterEndo])
    }
    
    
    ls.lm[[iterEndo]] <- lm(observed ~ predicted, data = df, subset = which(!is.na(df$observed)))
    
    par(mfrow = mfrow, mar = mar, mgp = mgp, oma = oma)
    for(iterPlot in which){
      if(iterPlot == 0){
        plot(df$predicted, df$observed, xlab = "Fitted values", ylab = "Observed values", main = "Observed vs Fitted")
        lowessTempo <- lowess(df$predicted[!is.na(df$observed)], df$observed[!is.na(df$observed)])
        points(lowessTempo$x, lowessTempo$y, type = "l", col = "red")
      }else if(iterPlot == 2){
        car:::qqPlot(ls.lm[[iterEndo]], main = "QQ plot")
      }else{
        plot(ls.lm[[iterEndo]], which = iterPlot, sub.caption = "")
      }
    }
    R2 <- summary(ls.lm[[iterEndo]])$r.squared
    if(R2 < 10^{-digit}){R2 <- paste0("<",10^{-digit})}else{R2 <- round(R2, digit)}
    title(paste0("Diagnostic plots for ",iterEndo," - R2=",R2), outer = TRUE)
    
  }
  names(ls.lm) <- endogenous
  
  
  #### export
  return(invisible(ls.lm))
  
}
