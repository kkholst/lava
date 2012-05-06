##' Extract model
##' 
##' Extract or replace model object
##' 
##' 
##' @aliases Model Model<- Model.lvmfit Model<-.lvmfit Model.multigroupfit
##' @usage
##' Model(x, ...)
##' Model(x, ...) <- value
##' @S3method Model
##' @method 
##' Model<-.multigroupfit
##' @param x Fitted model
##' @param value New model object (e.g. \code{lvm} or \code{multigroup})
##' @param \dots Additional arguments to be passed to the low level functions
##' @return Returns a model object (e.g. \code{lvm} or \code{multigroup})
##' @author Klaus K. Holst
##' @seealso \code{\link{Graph}}
##' @keywords models
##' @examples
##' 
##' m <- lvm(y~x)
##' e <- estimate(m, sim(m,100))
##' Model(e)
##'
##' @export
`Model` <- function(x,...) UseMethod("Model")
`Model.lvm` <- function(x,...) x
`Model.lvmfit` <- function(x,...) x$model
`Model.multigroup` <- function(x,...) x$lvm
`Model.multigroupfit` <- function(x,...) x$model

"Model<-" <- function(x,...,value) UseMethod("Model<-")
"Model<-.lvm" <- function(x,...,value) { x <- value; return(x) }
"Model<-.lvmfit" <- function(x,...,value) { x$model <- value; return(x) }
"Model<-.multigroup" <- function(x,...,value) { x$lvm <- value; return(x) }
"Model<-.multigroupfit" <- function(x,...,value) { x$model <- value; return(x) }
