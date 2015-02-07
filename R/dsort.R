
##' Sort data according to columns in data frame
##'
##' @title Sort data frame
##' @param data Data frame
##' @param x variable to order by
##' @param ... additional variables to order by
##' @return data.frame
##' @export
##' @examples
##' data(hubble)
##' dsort(hubble, "sigma")
##' dsort(hubble, hubble$sigma)
##' dsort(hubble,sigma,v)
dsort <- function(data,x,...) {
    e <- substitute(x)
    expr <- suppressWarnings(inherits(try(x,silent=TRUE),"try-error"))
    x1 <- if (expr) eval(e,envir=data) else {
        if (is.character(x) && length(x)<nrow(data)) lapply(x,function(z) data[,z]) else x
    }
    dots <- as.list(substitute(list(...))[-1])
    args <- lapply(dots, function(x) {
        e <- eval(x,envir=data)
        if (length(e)==1 && is.character(e)) e <- data[,e]
        e
    })
    if (!is.list(x1)) x1 <- list(x1)
    data[do.call("order",c(x1,args)),]
}
