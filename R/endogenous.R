##' @export
`endogenous` <-
    function(x,...) UseMethod("endogenous")

##' @export
`endogenous.lvmfit` <-
    function(x,...) {
        endogenous(Model(x),...)
    }

##' @export
`endogenous.lvm` <-
    function(x,top=FALSE,latent=FALSE,...) {
        observed <- manifest(x)
        if (latent) observed <- vars(x)
        if (top) {
            M <- x$M
            res <- c()
            for (i in observed)
                if (!any(M[i,]==1))
                    res <- c(res, i)
            return(res)
        }
        exo <- exogenous(x)
        return(setdiff(observed,exo))
    }

##' @export
endogenous.list <- function(x,...) {
    endolist <- c()
    for (i in seq_along(x)) {
        endolist <- c(endolist, endogenous(x[[i]]))
    }
    endolist <- unique(endolist)
    return(endolist)
}

##' @export
`endogenous.multigroup` <-
    function(x,...) {
        endogenous(Model(x))
    }

##' @export
`endogenous.lm` <-
    function(x,...) {
        getoutcome(formula(x))[1]
    }
