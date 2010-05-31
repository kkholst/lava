
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
