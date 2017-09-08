##' Split data into folds
##'
##' @title Split data into folds
##' @param x Data or integer (size)
##' @param p Number of folds, or if a number between 0 and 1 is given two folds of size p and (1-p) will be returned
##' @param replace With or with-out replacement
##' @param return.index If TRUE index of folds are returned otherwise the actual data splits are returned (default)
##' @param ... additional arguments to lower level functions
##' @export
##' @examples
##' csplit(10,3)
##' csplit(iris[1:10,]) ## Split in two sets 1:(n/2) and (n/2+1):n
##' csplit(iris[1:10,],0.5) 
##' @author Klaus K. Holst
csplit <- function(x,p=NULL,replace=FALSE,return.index=FALSE,...) {
    if (length(x)==1 & is.numeric(x)) x <- seq(x)
    N <- NROW(x)    
    if (is.null(p)) { ## 
        k1 <- base::round(N/2)
        idx <- list(base::seq(k1),

                   base::seq(k1+1,N))
    } else {
        if (p<1) { ## two folds (size N*p and N*(1-p))
            idx1 <- base::sample(N,base::round(p*N),replace=replace)
            idx <- list(idx1,                       
                       base::sample(setdiff(seq(N),idx1),replace=replace))
        } else { ## Number of folds (equal size)
            idx <- split(sample(seq(N)), rep(seq(p), length=N))
        }
    }
    if (return.index)
        return(idx)
    if (!is.vector(x)) {
        return(lapply(idx,function(ii) x[ii,,drop=FALSE]))
    }
    return(lapply(idx,function(ii) x[ii]))
}