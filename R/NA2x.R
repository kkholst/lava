##' Convert to/from NA
##'
##' Convert vector to/from NA
##' @aliases NA2x x2NA
##' @param s The input vector (of arbitrary class)
##' @param x The elements to transform into \code{NA} resp. what to transform
##' \code{NA} into.
##' @return A vector with same dimension and class as \code{s}.
##' @author Klaus K. Holst
##' @keywords manip
##' @examples##' 
##' x2NA(1:10, 1:5)
##' NA2x(x2NA(c(1:10),5),5)##'
##' @export
NA2x <- function(s,x=0) { sapply(s, function(y) ifelse(is.na(y),x,y) ) }
##' @export
x2NA <- function(s,x=0) { sapply(s, function(y) ifelse(y%in%x,NA,y) ) }
