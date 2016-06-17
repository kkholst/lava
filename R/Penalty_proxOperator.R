proxL1 <- function(x, step, lambda, test.penalty){
  if(test.penalty){
    max(0, (abs(x) - lambda * step)) * sign(x)
  }else{
    x
  }
}

proxL2 <- function(x, step, lambda, test.penalty){
  if(test.penalty){
    1 / (1 + lambda * step) * x
  }else{
    x
  }
}

proxE2 <- function(x, step, lambda){ # adapted from Simon 2013
  max(0, 1 - sqrt(length(x)) * lambda * step/norm(x, type = "2"))*x
}



init.proxOperator <- function(lambda1, lambda2, group.penaltyCoef, regularizationPath){
  
  #### No penalty
  if(all(lambda1 == 0) && all(lambda2 == 0) &&  regularizationPath == FALSE){ 
    
    proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2, expX){x}
    
    objectivePenalty <- function(x, lambda1, lambda2, test.penalty1, test.penalty2, group.penaltyCoef, expX){
      0
    }
    
    return(list(proxOperator = proxOperator,
                objectivePenalty = objectivePenalty))
  }
  
  #### Grouped lasso 
  if(any(group.penaltyCoef>=1)){ 
    
    if(all(lambda2 == 0)){
      
      proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2, expX){ ## proximal operator
        #if(!is.null(expX)){x[expX] <- exp(x[expX])}
        
        levels.penalty <- setdiff(unique(test.penalty1),0)
        for(iter_group in 1:length(levels.penalty)){
          index_group <- which(test.penalty1==levels.penalty[iter_group])
          x[index_group] <- proxE2(x = x[index_group], step = step, lambda = mean(lambda1[index_group])) # normally lambda1 has the same value for each member of the group
        }
        return(x)
      }
      
      objectivePenalty <- function(x, lambda1, lambda2, test.penalty1, test.penalty2, group.penaltyCoef, expX){
        #if(!is.null(expX)){x[expX] <- exp(x[expX])}
        
        levels.penalty <- setdiff(unique(test.penalty1),0)
        res <- 0
        for(iter_group in 1:length(levels.penalty)){
          index_group <- which(test.penalty1==levels.penalty[iter_group])
          res <- res +  mean(lambda1[index_group]) * norm(x[index_group], type = "2") # normally lambda1 has the same value  for each member of the group
        }
        
        return(res)
        
      }
    }else{
      stop("init.proxOperator: proximal operator for group lasso penalty and ridge penalty not yet defined \n")
    }
    
    return(list(proxOperator = proxOperator,
                objectivePenalty = objectivePenalty))
    
  }
  
  #### Standard lasso-ridge penalty 
  if(all(lambda2 == 0)){ ## lasso penalty
    
    proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2, expX){
      #if(!is.null(expX)){x[expX] <- exp(x[expX])}
      mapply(proxL1, x = x, step = step, lambda = lambda1, test.penalty = test.penalty1)
    }
    
    objectivePenalty <- function(x, lambda1, lambda2, test.penalty1, test.penalty2, group.penaltyCoef, expX){
      #if(!is.null(expX)){x[expX] <- exp(x[expX])}
      sum(lambda1[test.penalty1>0] * abs(x[test.penalty1>0]))
    }
    
  }else if(all(lambda1 == 0) &&  regularizationPath == FALSE){ ## ridge penalty
    
    proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2, expX){
      # if(!is.null(expX)){x[expX] <- exp(x[expX])}
      mapply(proxL2, x = x, step = step, lambda = lambda2, test.penalty = test.penalty2)
    }
    
    objectivePenalty <- function(x, lambda1, lambda2, test.penalty1, test.penalty2, group.penaltyCoef, expX){
      #if(!is.null(expX)){x[expX] <- exp(x[expX])}
      sum(lambda2[test.penalty2>0]/2 *(x[test.penalty2>0])^2)
    }
    
  }else{ ## elastic net penalty
    
    proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2, expX){
      #if(!is.null(expX)){x[expX] <- exp(x[expX])}
      mapply(proxL2, 
             x = mapply(proxL1, x = x, step = step,  lambda = lambda1, test.penalty = test.penalty1),
             step = step, lambda = lambda2, test.penalty = test.penalty2)
    }
    
    objectivePenalty <- function(x, lambda1, lambda2, test.penalty1, test.penalty2, group.penaltyCoef, expX){
      #if(!is.null(expX)){x[expX] <- exp(x[expX])}
      sum(lambda1[test.penalty1>0] * abs(x[test.penalty1>0])) + sum(lambda2[test.penalty2>0]/2 * (x[test.penalty2>0])^2)
    }
    
  }
 
  return(list(proxOperator = proxOperator,
              objectivePenalty = objectivePenalty))
  
}