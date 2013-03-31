
### Finds the unique commutation matrix K:
### \code{K%*%as.vector(A) = as.vector(t(A))}
## '
## ' @title Finds the unique commutation matrix
## ' @param m rows
## ' @param n columns
## ' @author Klaus K. Holst
commutation <- function(m, n=m) { 
  H <- function(i,j) { ## mxn-matrix with 1 at (i,j)
    Hij <- matrix(0, nrow=m, ncol=n)
    Hij[i,j] <- 1
    Hij
  }
  K <- matrix(0,m*n,m*n)  
  for (i in 1:m)
    for (j in 1:n)
      K <- K + H(i,j)%x%t(H(i,j))
  ## if (sparse)
  ##   return(as(K, "sparseMatrix"))
  K  
}
