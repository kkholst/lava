##' Generates missing entries in data.frame/matrix
##'
##' @title Create random missing data
##' @param data data.frame
##' @param p Fraction of missing data in each column
##' @param cols Which columns (name or index) to alter
##' @param rowwise Should missing occur row-wise (either none or all selected columns are missing)
##' @param nafun (Optional) function to be applied on data.frame before return (e.g. \code{na.omit} to return complete-cases only)
##' @param seed Random seed
##' @return data.frame
##' @author Klaus K. Holst
##' @keywords utilities
##' @export
makemissing <- function(data,p=0.2,cols=seq_len(ncol(data)),rowwise=FALSE,nafun=function(x) x, seed=NULL) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }    
    p <- rep(p,length.out=length(cols))
    if (!rowwise)
        for (i in seq_along(cols)) {
            data[rbinom(nrow(data),1,p[i])==1,cols[i]] <- NA
        }
    else
        data[which(rbinom(nrow(data),1,p)==1),cols] <- NA
    return(nafun(data))
}
