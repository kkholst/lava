##' Wait for user input (keyboard or mouse)
##'
##' @title Wait for user input (keyboard or mouse)
##' @aliases waitclick
##' @author Klaus K. Holst
##' @export
##' @keywords iplot
wait <- function() {
    cat(gettext("\nPress <Enter> to continue..."))
    res <- try(scan("", what=0, quiet=TRUE, nlines=1), silent=TRUE)
}
waitclick <- function() if(is.null(locator(1))) invisible(NULL)

