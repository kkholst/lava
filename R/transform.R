##' @export
"transform<-" <- function(`_data`,...,value) UseMethod("transform<-")

##' @export
"transform<-.lvm" <- function(`_data`,formula=NULL,...,value) {
    transform(`_data`,formula,value,...)
}


"transform.lvm" <- function(`_data`,formula,fun,post=TRUE,y,x,...) {
    if (!missing(y) && !missing(x)) {
        xx <- x
    } else {
        if (is.character(formula)) {
            y <- NULL; xx <- formula
        } else {
            y <- getoutcome(formula)
            xx <- attributes(y)$x
        }
    }
    if (length(xx)==0) { xx <- y; y <- NULL }
    if (length(y)==0) {
        if (post) {
            `_data`$constrainY[xx] <- NULL
            `_data`$constrain[xx] <- NULL
            if (is.null(`_data`$attributes$selftransform))
                `_data`$attributes$selftransform <- list()
            `_data`$attributes$selftransform[[xx]] <- fun
            return(`_data`)
        }
        `_data`$attributes$selftransform[xx] <- NULL
        constrain(`_data`,xx,y,...) <- fun
        return(`_data`)
    }

    `_data`$attributes$selftransform[y] <- NULL
    addvar(`_data`) <- c(y)
    intercept(`_data`,y) <- 0; covariance(`_data`,y) <- 0
    if (is.null(`_data`$attributes$transform))
        `_data`$attributes$transform <- list()
    if (is.null(fun)) `_data`$attributes$transform[y] <- NULL
    else
        `_data`$attributes$transform[[y]] <- list(fun=fun,x=xx)
    return(`_data`)
}
