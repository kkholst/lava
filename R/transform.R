##' @export
"transform<-" <- function(`_data`,...,value) UseMethod("transform<-")

##' @export
"transform<-.lvm" <- function(`_data`,formula=NULL,...,value) {
    transform(`_data`,formula,value,...)
}

##' @export
print.transform.lvm <- function(x,...) {
    for (i in seq_along(x)) {
        cat("Variable: ", names(x)[i],"\n",sep="")
        cat("Transformation: (",paste0(x[[i]]$x,collapse=","),") -> ",sep="")
        print(x[[i]]$fun)
        cat("\n")
    }
    invisible(x)
}

##' @export
"transform.lvm" <- function(`_data`,formula,value,post=TRUE,y,x,...) {
    if (missing(formula)) {
        if (length(tr <- `_data`$attributes$transform)==0) {
            return(NULL)
        }        
        return(structure(`_data`$attributes$transform,class="transform.lvm"))
    }
    
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
            `_data`$attributes$selftransform[[xx]] <- value
            return(`_data`)
        }
        `_data`$attributes$selftransform[xx] <- NULL
        constrain(`_data`,xx,y,...) <- value
        return(`_data`)
    }
    
    
    `_data`$attributes$selftransform[y] <- NULL
    addvar(`_data`) <- y
    intercept(`_data`,y) <- 0; covariance(`_data`,y) <- 0
    if (is.null(`_data`$attributes$transform))
        `_data`$attributes$transform <- list()
    if (is.null(value)) `_data`$attributes$transform[y] <- NULL
    else {
        if (length(y)>1) {
            if (is.null(`_data`$attributes$multitransform))
                `_data`$attributes$multitransform <- list()
            `_data`$attributes$multitransform
            for (yi in y) {
                `_data`$attributes$transform[yi] <- NULL
            }
            rmidx <- c()
            for (i in seq_along(`_data`$attributes$multitransform)) {
                l <- `_data`$attributes$multitransform[[i]]
                if (any(y%in%letters)) rmidx <- c(rmidx,i)
            }
            if (length(rmidx)>0) `_data`$attributes$transform[rmidx] <- NULL            
            `_data`$attributes$multitransform <- c(`_data`$attributes$multitransform,                                                   
                                                   list(list(fun=value,y=y,x=xx)))
        } else {
            `_data`$attributes$transform[[y]] <- list(fun=value,x=xx)
        }
    }
    return(`_data`)
}


addhook("plothook_transform","plot.post.hooks")

plothook_transform <- function(x,...) {
    trans <- x$attributes$transform
    transnames <- names(trans)
    for (v in transnames) {
        xx <- trans[[v]][["x"]]
        if (length(xx)>0) {
            x <- regression(x,x=xx,y=v)
            edgelabels(x,from=xx,to=v,col="gray70") <- ""
        }        
    }
    return(x)
}

