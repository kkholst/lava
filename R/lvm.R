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

