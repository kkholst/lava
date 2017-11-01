
##' @export
estimate.list <- function(x,...) {
    if (inherits(x[[1]],"lvm")) return(estimate.lvmlist(x,...))
    res <- lapply(x,function(x) estimate(x,...))
    class(res) <- c("estimate.list","list")
    res
}

coef.estimate.list <- function(object,...) {
    lapply(object,coef)
}


