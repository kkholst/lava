
##' Sort data according to columns in data frame
##'
##' @title Sort data frame
##' @param data Data frame
##' @param ... Arguments to order 
##' @return data.frame 
##' @export
##' @examples
##' data(hubble)
##' dsort(hubble, "sigma")
##' dsort(hubble, hubble$sigma)
##' dsort(hubble,sigma,v)
dsort <- function(data,...) {
    dots <- as.list(substitute(list(...))[-1])
    args <- lapply(dots, function(x) {        
        e <- eval(x,envir=data)
        if (length(e)==1 && is.character(e)) e <- data[,e]
        e
    })
    data[do.call("order",args),]    
}
