###{{{ %+% concat operator

##' Concatenation operator
##'
##' For matrices a block-diagonal matrix is created. For latent variable models ('lvm' objects) \code{merge.lvm} is called. For all other data types he operator is a wrapper of \code{paste}.
##' @name concatoperator
##' @title Concatenation operator
##' @param x First object
##' @param y Second object
##' @author Klaus K. Holst
##' @export
`%+%` <- function(x,y) UseMethod("%+%",y)

##' @S3method %+% lvm
`%+%.lvm` <- function(x,y) merge(x,y)

##' @S3method %+% matrix
`%+%.matrix` <- function(x,y) blockdiag(x,y)

##' @S3method %+% default
`%+%.default` <- function(x,y) paste(x, y, sep="")

###}}}


##' Matching operator
##' 
##' Matching operator (x not in y) oposed to the \code{\%in\%}-operator (x in y).
##' @rdname op_match
##' @aliases %nin%
##' @param x vector
##' @param y vector of same type as \code{x}
##' @return A logical vector.
##' @author Klaus K. Holst
##' @seealso \code{\link{match}}
##' @keywords utilities misc
##' @examples
##' 
##' \dontrun{1:10 %nin% c(1,5,10)}
##'
##' @export
`%nin%` <-
function(x,y) {
  is.na(match(x,y))
}

