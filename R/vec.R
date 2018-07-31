##' vec operator
##'
##' Convert array into vector
##' @title vec operator
##' @param x Array
##' @param matrix If TRUE a row vector (matrix) is returned
##' @param sep Seperator
##' @param ... Additional arguments
##' @author Klaus Holst
##' @export
vec <- function(x,matrix=FALSE,sep=".",...) {
    if (is.vector(x) && !is.list(x)) {
        res <- x
    } else if (is.list(x)) {
        res <- stats::setNames(unlist(x),names(x))
    } else {
        if (is.matrix(x) && is.null(rownames(x))) {
            nn <- colnames(x)
        } else {
            nn <- apply(expand.grid(dimnames(x)),1,function(x) paste(x,collapse=sep))
        }
        res <- as.vector(x); names(res) <- nn
    }
    if (matrix) return(cbind(res))
    return(res)
}
