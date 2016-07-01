#' @title Compute the degree of freedom for nuclear regularized regression
#' 
#' @references Zhou et al. Regularized Matrix Regression
#' 
dfNuclear <- function(B.LS, B.lambda, lambda, lambda2 = 0, tol = 1e-12){
  
  eigen.LS <- svd(B.LS)$d
  if(any(eigen.LS<0)){
    warning("dfNuclear: the ",if(lambda2>0){"regularized"}," least square solution has negative eigenvalues \n",
            "the formula for df may not be valid \n")
  }
  eigen.lambda <- svd(B.lambda)$d
  q <- length(eigen.lambda)
  p1 <- nrow(B.lambda)
  p2 <- ncol(B.lambda)
  seq_q <- which(eigen.lambda > tol)
  
  ####
  if(length(seq_q)==0){
    
    return(0)
    
  }else{
    
    #### upper bound
    #   rank.lambda <- length(seq_q)
    #   rank.lambda*(p1+p2) - rank.lambda^2
    
    #### explicit formula
    eigen.LS0 <- c(eigen.LS, rep(0, max(p1,p2) - q))
    eigen.lambda0 <- c(eigen.lambda, rep(0, max(p1,p2) - q))
    vec.lambda <- eigen.LS-eigen.lambda0
    
    vec.df <- sapply(seq_q, function(iter_q){
      #sumTempo <- eigen.LS[iter_q] * ((1 + lambda2) * eigen.LS0[iter_q] - lambda) / (eigen.LS0[iter_q]^2 - eigen.LS0^2)
      sumTempo <- eigen.LS[iter_q] * ((1 + lambda2) * eigen.LS0[iter_q] - vec.lambda[iter_q]) / (eigen.LS0[iter_q]^2 - eigen.LS0^2)
      sum1 <- sum(sumTempo[setdiff(1:p1,iter_q)])
      sum2 <- sum(sumTempo[setdiff(1:p2,iter_q)])
      
      return(1 + 1 / (1 + lambda2) * sum1 + 1 / (1 + lambda2) * sum2)
    })
    df <- sum(vec.df)
    
    return(df)
  }
}