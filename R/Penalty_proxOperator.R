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


proxNuclear <-  function(x, step, lambda, nrow, ncol){
  eigen.Mx <- svd(matrix(x, nrow = nrow, ncol = ncol))
  n.eigen <- min(nrow, ncol)
  b <- mapply(proxL1, x = eigen.Mx$d, step = step, lambda = rep(lambda, n.eigen), test.penalty = rep(1, n.eigen))
  as.vector(eigen.Mx$u %*% diag(b) %*% t(eigen.Mx$v))
}


init.proxOperator <- function(coef,
                              lambda1, lambda2, group.coef, 
                              lambdaN, nrow, ncol, 
                              regularizationPath){
  
  proxOperator <- list()
  objective <- list()
  
  test.N <- !is.null(lambdaN) && lambdaN != 0
  test.L1 <- (!is.null(lambda1) && any(lambda1 > 0)) || (regularizationPath >= 0)
  if(test.N){test.L1 <- FALSE;test.L2 <- FALSE}
  test.L2 <- !is.null(lambda2) && any(lambda2 > 0)
  
  #### No penalty
  if(test.N == FALSE && test.L1 == FALSE &&  test.L2 == FALSE){ 
    return(list(proxOperator = function(x, ...){x},
                objectivePenalty = list(function(...){0})))
  }
  
  #### Nuclear norm penalty
  if(test.N){
    proxOperator$N <-  function(x, step, test.penaltyN, nrow, ncol, ...){
      x[test.penaltyN] <- proxNuclear(x = x[test.penaltyN], step = step, lambda = lambdaN, nrow = nrow, ncol = ncol)
      return(x)
    }
    objective$N <- function(x, step, test.penaltyN, nrow, ncol, ...){
      sum(svd(matrix(x[test.penaltyN], nrow = nrow, ncol = ncol))$d)
    }
  }else{
    objective$N <- function(...){0}
  }

  #### Lasso penalty
  if(test.L1){
    if(any(group.coef>=1)){ ## group lasso
      
      proxOperator$L1 <- function(x, step, lambda1, test.penalty1, expX, ...){
        
        levels.penalty <- setdiff(unique(test.penalty1),0)
        for(iter_group in 1:length(levels.penalty)){
          index_group <- which(test.penalty1==levels.penalty[iter_group])
          x[index_group] <- proxE2(x = x[index_group], step = step, lambda = mean(lambda1[index_group])) # normally lambda1 has the same value for each member of the group
        }
        return(x)
      }
      
      objective$L1 <- function(x, lambda1, lambda2, test.penalty1, test.penalty2, group.coef, expX){
        if(!is.null(expX)){x[expX] <- exp(x[expX])}
        
        levels.penalty <- setdiff(unique(test.penalty1),0)
        res <- 0
        for(iter_group in 1:length(levels.penalty)){
          index_group <- which(test.penalty1==levels.penalty[iter_group])
          res <- res +  mean(lambda1[index_group]) * norm(x[index_group], type = "2") # normally lambda1 has the same value  for each member of the group
        }
        return(res)
      }
      
    }else{ ## normal lasso
      proxOperator$L1 <- function(x, step, lambda1, test.penalty1, expX, ...){
        mapply(proxL1, x = x, step = step, lambda = lambda1, test.penalty = test.penalty1)
      }
    }
    
    objective$L1 <- function(x, lambda1, test.penalty1, group.coef, expX, ...){
      if(!is.null(expX)){x[expX] <- exp(x[expX])}
      sum(lambda1[test.penalty1>0] * abs(x[test.penalty1>0]))
    }
  }else{
    objective$L1 <- function(...){0}
  }
  
  #### Ridge penalization
  if(test.L2){
    proxOperator$L2 <- function(x, step, lambda2, test.penalty2, expX, ...){
      mapply(proxL2, x = x, step = step, lambda = lambda2, test.penalty = test.penalty2)
    }
    
    objective$L2 <- function(x, lambda1, test.penalty2, group.coef, expX, ...){
      if(!is.null(expX)){x[expX] <- exp(x[expX])}
      sum(lambda2[test.penalty2>0]/2 *(x[test.penalty2>0])^2)
    }
  }else{
    objective$L2 <- function(...){0}
  }
  
  ls.prox <- list(if(test.N){proxOperator$N},if(test.L1){proxOperator$L1},if(test.L2){proxOperator$L2})
  ls.prox <- ls.prox[!sapply(ls.prox,is.null)]
  
  ls.obj <- list(if(test.N){objective$N},if(test.L1){objective$L1},if(test.L2){objective$L2})
  ls.obj <- ls.obj[!sapply(ls.obj,is.null)]
  
  return(list(proxOperator = do.call(composeOperator, args = ls.prox),
              objectivePenalty = ls.obj
  ))
}


#### miscelaneous function ####
composeOperator <- function (...){ 
  ls.fct <- lapply(list(...), match.fun)
  
  newfct <- function(x, lambda1, lambda2, index.constrain = NULL, type.constrain, expX, ...) {
   
    if(length(index.constrain)>0){
      if(type.constrain){norm <- sum(exp(x[index.constrain]))
      }else{
        norm <- sum(x[index.constrain])
        if(norm < 0){stop("proxGrad: negative variance parameter - set constrain to TRUE in control \n")}
      }
      lambda1 <- lambda1/norm
      lambda2 <- lambda2/norm
    }
    
    for (f in ls.fct) {
      if(!is.null(expX)){x[expX] <- exp(x[expX])}
      x <- f(x, lambda1 = lambda1, lambda2 = lambda2, ...)
      if(!is.null(expX)){x[expX] <- log(x[expX])}
    }
    return(x)
  }
  
  return(newfct)
}
