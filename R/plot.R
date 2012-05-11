###{{{ plot.lvm



##' Plot path diagram
##' 
##' Plot the path diagram of a SEM
##' 
##' 
##' @aliases plot.lvmfit
##' @param x Model object
##' @param diag Logical argument indicating whether to visualize variance
##' parameters (i.e. diagonal of variance matrix)
##' @param cor Logical argument indicating whether to visualize correlation
##' parameters
##' @param labels Logical argument indiciating whether to add labels to plot
##' (Unnamed parameters will be labeled p1,p2,...)
##' @param intercept Logical argument indiciating whether to add intercept
##' labels (current version: not used))
##' @param addcolor Logical argument indiciating whether to add colors to plot
##' (overrides \code{nodecolor} calls)
##' @param plain if TRUE strip plot of colors and boxes
##' @param cex Fontsize of node labels
##' @param fontsize1 Fontsize of edge labels
##' @param noplot if TRUE then return \code{graphNEL} object only
##' @param graph Graph attributes (Rgraphviz)
##' @param attrs Attributes (Rgraphviz)
##' @param unexpr if TRUE remove expressions from labels
##' @param addstyle Logical argument indicating whether additional style should
##' automatically be added to the plot (e.g. dashed lines to double-headed
##' arrows)
##' @param Rgraphviz if FALSE igraph is used for graphics
##' @param init Reinitialize graph (for internal use)
##' @param \dots Additional arguments to be passed to the low level functions
##' @author Klaus K. Holst
##' @keywords hplot regression
##' @examples
##' 
##' \dontrun{
##' example(estimate)
##' plot(e)
##' }
##'
##' @S3method plot lvm
##' @method plot lvm
`plot.lvm` <-
  function(x,diag=FALSE,cor=TRUE,labels=FALSE,intercept=FALSE,addcolor=TRUE,plain=FALSE,cex,fontsize1=10,noplot=FALSE,graph=list(rankdir="BT"),
         attrs=list(graph=graph),
           unexpr=FALSE,
           addstyle=TRUE,Rgraphviz=lava.options()$Rgraphviz,init=TRUE,
           ...) {
  index(x) <- reindex(x)
  if (length(index(x)$vars)<2) stop("Not available for models with fewer than two variables")

  if (!Rgraphviz || (!require("Rgraphviz"))) {
    if (!require("igraph"))
      stop("package 'Rgraphviz' or 'igraph' not available")
    g <- igraph.lvm(x,...)
    if (noplot) return(g)
    plot(g)
    return(invisible(g))
  } 
  if (init) {
    g <- finalize(x,diag=diag,cor=cor,addcolor=addcolor,intercept=intercept,plain=plain,cex=cex,fontsize1=fontsize1,unexpr=unexpr,addstyle=addstyle)
  } else {
    g <- Graph(x)
  }
  if  (labels) {
    AP <- matrices(x,paste("p",seq_len(index(x)$npar),sep=""))
    mylab <- AP$P; mylab[AP$A!="0"] <- AP$A[AP$A!="0"]
    mylab[!is.na(x$par)] <- x$par[!is.na(x$par)]
    mylab[!is.na(x$covpar)] <- x$covpar[!is.na(x$covpar)]
    g <- edgelabels(g, lab=mylab)
  }
  if (lava.options()$debug) {
    plot(g)
  } else {

##    graphRenderInfo(g)$recipEdges <- "distinct"
    .savedOpt <- options(warn=-1) ## Temporarily disable warnings as renderGraph comes with a stupid warning when labels are given as "expression"
    dots <- list(...)
    dots$attrs <- attrs
    dots$x <- g
    dots$recipEdges <- "distinct"
    if (is.null(dots$layoutType) & all(index(x)$A==0))
      dots$layoutType <- "circo"
    g <- do.call("layoutGraph", dots)
    if (noplot)
      return(g)    
    res <- tryCatch(renderGraph(g),error=function(e) NULL)
    options(.savedOpt)
  }
  ## if (!is.null(legend)) {
  ##   op <- par(xpd=TRUE)
  ##   legend(legend, c("Exogenous","Endogenous","Latent","Time to event"),
  ##          pt.cex=1.5, pch=15, lty=0, col=cols[1:4], cex=0.8)
  ##   par(op)
  ## }
  invisible(g)
}

###}}} plot.lvm

###{{{ plot.lvmfit

##' @S3method plot lvmfit
`plot.lvmfit` <-
  function(x,diag=TRUE,cor=TRUE,type,noplot=FALSE,fontsize1=5,...) {
    .savedOpt <- options(warn=-1) ## Temporarily disable warnings as renderGraph comes with a stupid warning when labels are given as "expression"
    g <- Graph(x)
    newgraph <- FALSE
    if (is.null(g)) {
      newgraph <- TRUE
      Graph(x) <- finalize(Model(x), diag=TRUE, cor=FALSE, fontsize1=fontsize1, ...)
    }
    if(noplot) return(Graph(x))
    ##  browser()
    ##    newgraph <- FALSE
    ##cat("Setting up graph...\n")
    if (newgraph) {
      if (missing(type))
        type <- "est"
      x <- edgelabels(x, type=type, diag=diag, cor=cor, fontsize1=fontsize1, ...)
    } else {
      if (!missing(type)) {
        x <- edgelabels(x, type=type, diag=diag, cor=cor, fontsize1=fontsize1, ...)
      }
    }
    g <- Graph(x)
    var <- rownames(covariance(Model(x))$rel)
    if (!cor) {
       delta <- 1
      for (r in 1:(nrow(covariance(Model(x))$rel)-delta) ) {
        for (s in (r+delta):ncol(covariance(Model(x))$rel) ) {
          if (covariance(Model(x))$rel[r,s]==1) {
            g <- removeEdge(var[r],var[s], g)
            g <- removeEdge(var[s],var[r], g)
          }
        }
      }
    }
    if (!diag) {
      for (r in 1:(nrow(covariance(Model(x))$rel)) ) {
        if (isAdjacent(g,var[r],var[r]))
          g <- removeEdge(var[r],var[r],g)
      }
    } 
    m <- Model(x); Graph(m) <- g
    g <- plot(m, diag=diag, cor=cor, fontsize1=fontsize1, init=FALSE, ...)
    options(.savedOpt)
    invisible(g)    
  }

###}}} plot.lvmfit

###{{{ plot.multigroup

##' @S3method plot multigroup
plot.multigroup <- function(x,diag=TRUE,labels=TRUE,...) {
  k <- x$ngroup
  for (i in 1:k)
    plot(x$lvm[[i]],diag=diag,labels=labels, ...)
}

##' @S3method plot multigroupfit
plot.multigroupfit <- function(x,...) {
  plot(Model(x),...)
}

###}}}

###{{{ igraph.lvm

##' @export
igraph.lvm <- function(x,layout=igraph::layout.kamada.kawai,...) {
  require("igraph")
  oC <- covariance(x)$rel
  for (i in 1:(nrow(oC)-1))
    for (j in (i+1):nrow(oC)) {
      if (oC[i,j]!=0) {
        x <- regression(x,vars(x)[i],vars(x)[j])
        x <- regression(x,vars(x)[j],vars(x)[i])
      }
    }  
  g <- igraph.from.graphNEL(Graph(x))  
  V(g)$color <- "lightblue"
  V(g)$label <- vars(x)
  for (i in match(latent(x),V(g)$name)) {
    V(g)$shape[i] <- "circle"
    V(g)$color[i] <- "green"
  }
  endo <- index(x)$endogenous
  for (i in match(endo,V(g)$name)) {
    V(g)$color[i] <- "orange"
  }
  E(g)$label <- as.list(rep("",length(E(g))))
  oE <- edgelabels(x)
  for (i in 1:length(E(g))) {
    st <- as.character(oE[i])
    if (length(st)>0)
      E(g)$label[[i]] <- st
  }
  g$layout <- layout(g)
  return(g)  
}

###}}} igraph.lvm
