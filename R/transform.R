##' @export
"transform<-" <- function(x,...,value) UseMethod("transform<-")

##' @S3method transform<- lvm
"transform<-.lvm" <- function(x,formula,...,value) {
    transform(x,formula,value,...)
}


"transform.lvm" <- function(`_data`,formula,fun,post=TRUE,...) {
    if (is.character(formula)) {
        y <- NULL; xx <- formula
    } else {
        y <- getoutcome(formula)
        xx <- attributes(y)$x
    }    
    if (length(xx)==0) { xx <- y; y <- NULL }
    if (length(y)==0) {
        if (post) {
            `_data`$constrainY[xx] <- NULL
            `_data`$constrain[xx] <- NULL
            if (is.null(attributes(`_data`)$selftransform))
                attributes(`_data`)$selftransform <- list()
            attributes(`_data`)$selftransform[[xx]] <- fun
            return(`_data`)
        }
        attributes(`_data`)$selftransform[xx] <- NULL
        constrain(`_data`,xx,y,...) <- fun
        return(`_data`)
    }

    attributes(`_data`)$selftransform[y] <- NULL
    addvar(`_data`) <- c(y)
    intercept(`_data`,y) <- 0; covariance(`_data`,y) <- 0
    if (is.null(attributes(`_data`)$transform))
        attributes(`_data`)$transform <- list()
    if (is.null(fun)) attributes(`_data`)$transform[y] <- NULL
    else
        attributes(`_data`)$transform[[y]] <- list(fun=fun,x=xx)  
    return(`_data`)
}
