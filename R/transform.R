##' @export
"transform<-" <- function(x,...,value) UseMethod("transform<-")

##' @S3method transform<- lvm
"transform<-.lvm" <- function(x,formula,...,value) 
    transform(x,formula,value,...)


"transform.lvm" <- function(`_data`,formula,fun,...) {
    if (is.character(formula)) {
        y <- NULL; xx <- formula
    } else {
        y <- getoutcome(formula)
        xx <- attributes(y)$x
    }
    if (length(y)==0) {
        if (is.null(attributes(`_data`)$selftransform))
            attributes(`_data`)$selftransform <- list()
        attributes(`_data`)$selftransform[[xx]] <- fun
        return(`_data`)
    }
    addvar(`_data`) <- c(y)
    intercept(`_data`,y) <- 0; covariance(`_data`,y) <- 0
    if (is.null(attributes(`_data`)$transform))
        attributes(`_data`)$transform <- list()
    if (is.null(fun)) attributes(`_data`)$transform[y] <- NULL
    else
        attributes(`_data`)$transform[[y]] <- list(fun=fun,x=xx)  
    return(`_data`)
}
