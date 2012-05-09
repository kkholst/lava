##' Initialize new latent variable model
##' 
##' Function that constructs a new latent variable model object
##' 
##' @aliases lvm print.lvm summary.lvm
##' @param x Vector of variable names. Optional but gives control of the
##' sequence of appearance of the variables. The argument can be given as a
##' character vector or formula, e.g. \code{~y1+y2} is equivalent to
##' \code{c("y1","y2")}. Alternatively the argument can be a formula specifying
##' a linear model.
##' @param \dots Additional arguments to be passed to the low level functions
##' @param silent Logical variable which indicates whether messages are turned
##' on/off.
##' @return Returns an object of class \code{lvm}.
##' @author Klaus K. Holst
##' @seealso \code{\link{regression}}, \code{\link{covariance}},
##' \code{\link{intercept}}, ...
##' @keywords models regression
##' @export
##' @examples
##' 
##' m <- lvm() # Empty model
##' m1 <- lvm(y~x) # Simple linear regression
##' m2 <- lvm(~y1+y2) # Model with two independent variables (argument)
##' m3 <- lvm(list(c(y1,y2,y3)~u,u~x+z)) # SEM with three items
##' 
lvm <- function(x=NULL, ..., silent=lava.options()$silent) {
 
  m <- new("graphNEL", edgemode="directed"); C <- par <- fix <- numeric(); mu <- list()
  graphRenderInfo(m)$recipEdges <- "distinct"
  nodeDataDefaults(m, "fill") <- "white"
  nodeDataDefaults(m, "shape") <- "rectangle"
  nodeDataDefaults(m, "latent") <- FALSE
  nodeDataDefaults(m, "randomslope") <- FALSE
  nodeDataDefaults(m, "normal") <- TRUE
  nodeDataDefaults(m, "survival") <- FALSE
  nodeDataDefaults(m, "parameter") <- FALSE
  nodeDataDefaults(m, "categorical") <- FALSE
  nodeDataDefaults(m, "distribution") <- NA
  nodeDataDefaults(m, "label") <- expression(NA)
  myhooks <- gethook("init.hooks")
  for (f in myhooks) {
    m <- do.call(f, list(x=m))
  }        


  edgeDataDefaults(m, "lty") <- 1
  edgeDataDefaults(m, "lwd") <- 1
  edgeDataDefaults(m, "col") <- "black"
  edgeDataDefaults(m, "textCol") <- "black"
  edgeDataDefaults(m, "est") <- 0
  edgeDataDefaults(m, "arrowhead") <- "open"
  edgeDataDefaults(m, "dir") <- "forward"
  edgeDataDefaults(m, "cex") <- 1.5
  edgeDataDefaults(m, "label") <- expression()
  edgeDataDefaults(m, "futureinfo") <- list()

  res <- list(graph=m, par=par, cov=C, covpar=C, fix=fix, covfix=fix, mean=mu, index=NULL, exogenous=NA, constrain=list())
  class(res) <- "lvm"


  myvar <- NULL
  lvar <- x
##  lvar <- list(x,...)
  if (!is.list(lvar)) lvar <- list(x)
  for (myvar in lvar) {
    if (class(myvar)[1]=="formula") {
      ## if (length(getoutcome(myvar))>0) {
      ##   regression(res,...,silent=silent) <- myvar     
      ## } else {
      ##   myvar <- all.vars(myvar)
      ## }
      regression(res,...,silent=silent) <- myvar
    }
    if (is.character(myvar)) {
      res <- addvar(res, myvar, silent=silent)  }
  }
  if (!is.null(myvar)) {
    index(res) <- reindex(res,zeroones=TRUE) }

    
  return(res)
}

