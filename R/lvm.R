lvm <- function(var=NULL, silent=FALSE, ...) {
 
  x <- new("graphNEL", edgemode="directed"); C <- par <- fix <- numeric(); mu <- list()
  nodeDataDefaults(x, "fill") <- "white"
  nodeDataDefaults(x, "shape") <- "rectangle"
  nodeDataDefaults(x, "latent") <- FALSE
  nodeDataDefaults(x, "randomslope") <- FALSE
  nodeDataDefaults(x, "normal") <- TRUE
  nodeDataDefaults(x, "survival") <- FALSE
  nodeDataDefaults(x, "parameter") <- FALSE
  nodeDataDefaults(x, "categorical") <- FALSE
  nodeDataDefaults(x, "distribution") <- NA
  nodeDataDefaults(x, "label") <- expression(NA)
  myhooks <- gethook("init.hooks")
  for (f in myhooks) {
    x <- do.call(f, list(x=x))
  }        


  edgeDataDefaults(x, "lty") <- 1
  edgeDataDefaults(x, "lwd") <- 1
  edgeDataDefaults(x, "col") <- "black"
  edgeDataDefaults(x, "textCol") <- "black"
  edgeDataDefaults(x, "est") <- 0
  edgeDataDefaults(x, "arrowhead") <- "open"
  edgeDataDefaults(x, "dir") <- "forward"
  edgeDataDefaults(x, "cex") <- 1.5
  edgeDataDefaults(x, "label") <- expression()
  edgeDataDefaults(x, "futureinfo") <- list()

  res <- list(graph=x, par=par, cov=C, covpar=C, fix=fix, covfix=fix, mean=mu, index=NULL, exogenous=NA, constrain=list())
  class(res) <- "lvm"


  myvar <- NULL
  lvar <- var
  if (!is.list(lvar)) lvar <- list(var)
  for (myvar in lvar) {
    if (class(myvar)[1]=="formula") {
      if (length(getoutcome(myvar))>0) {
        regression(res,...,silent=silent) <- myvar     
      } else {
        myvar <- all.vars(myvar)
      }
    }
    if (is.character(myvar)) {
      res <- addvar(res, myvar, silent=silent)  }
  }
  if (!is.null(myvar)) {
    index(res) <- reindex(res,zeroones=TRUE) }

    
  return(res)
}

