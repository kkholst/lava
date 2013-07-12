##' Create a Data Frame from All Combinations of Factors
##'
##' Simple wrapper of the 'expand.grid' function
##' @title Create a Data Frame from All Combinations of Factors
##' @param x Data.frame
##' @param ... vectors, factors or a list containing these
##' @author Klaus K. Holst
##' @export
##' @examples
##' dd <- Expand(iris, Sepal.Length=2:8, Species="virginica")
##' str(dd)
Expand <- function(x,...) {
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
