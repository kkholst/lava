##' Create a Data Frame from All Combinations of Factors
##'
##' Simple wrapper of the 'expand.grid' function.  If x is a table
##' then a data frame is returned with one row pr individual
##' observation. 
##' @title Create a Data Frame from All Combinations of Factors
##' @param x Data.frame
##' @param ... vectors, factors or a list containing these
##' @author Klaus K. Holst
##' @export
##' @examples
##' dd <- Expand(iris, Sepal.Length=2:8, Species="virginica")
##' str(dd)
##'
##' T <- with(warpbreaks, table(wool, tension))
##' Expand(T)
Expand <- function(x,...) {
    if (inherits(x,"table")) {
        M <- as.data.frame(x)
        idx <- rep(seq(nrow(M)),M[,ncol(M)])
        return(M[idx,-ncol(M),drop=FALSE])                
    }
    dots <- list(...)
    nn <- names(dots)
    for (n in nn) {
        y <- dots[[n]]
        if (is.factor(x[1,n])) {
            dots[[n]] <- factor(y,levels=levels(x[1,n]))
        }
    }
    do.call("expand.grid",dots)
}
