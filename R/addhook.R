
##' Set global options for \code{lava}
##'
##' Extract and set global parameters of \code{lava}. In particular optimization
##' parameters for the \code{estimate} function.
##'
##' \itemize{
##'   \item \code{param}: 'relative' (factor loading and variance of one
##' endogenous variables in each measurement model are fixed to one), 'absolute'
##' (mean and variance of latent variables are set to 0 and 1, respectively),
##' 'hybrid' (intercept of latent variables is fixed to 0, and factor loading of
##' at least one endogenous variable in each measurement model is fixed to 1),
##' 'none' (no constraints are added)
##'   \item \code{layout}: One of 'dot','fdp','circo','twopi','neato','osage'
##'   \item \code{messages}: Set to 0 to disable various output messages
##'   \item ...  }
##'
##' see \code{control} parameter of the \code{estimate} function.
##'
##' @param \dots Arguments
##' @return \code{list} of parameters
##' @author Klaus K. Holst
##' @keywords models
##' @examples
##'
##' \dontrun{
##' lava.options(iter.max=100,messages=0)
##' }
##'
##' @export
lava.options <- function(...) {
    dots <- list(...)
    newopt <- curopt <- get("options",envir=lava.env)
    if (length(dots)==0)
        return(curopt)
    if (length(dots)==1 && is.list(dots[[1]]) && is.null(names(dots))) {
        dots <- dots[[1]]
    }
    idx <- which(names(dots)!="")
    newopt[names(dots)[idx]] <- dots[idx]
    assign("options",newopt,envir=lava.env)
    invisible(curopt)
}

##' @export
gethook <- function(hook="estimate.hooks",...) {
    get(hook,envir=lava.env)
}

##' @export
addhook <- function(x,hook="estimate.hooks",...) {
    newhooks <- unique(c(gethook(hook),x))
    assign(hook,newhooks,envir=lava.env)
    invisible(newhooks)
}

versioncheck <- function(pkg="lava",geq,sep=".",...) {
    xyz <- tryCatch(
        char2num(strsplit(as.character(packageVersion(pkg)),sep,fixed=TRUE)[[1]]),
        error=function(x) NULL)
    if (is.null(xyz)) return(FALSE)
    if (missing(geq)) return(xyz)
    for (i in seq(min(length(xyz),length(geq)))) {
        if (xyz[i]>geq[i]) return(TRUE)
        if (xyz[i]<geq[i]) return(FALSE)
    }
    if (length(xyz)>=length(geq)) return(TRUE)
    return(FALSE)
}

lava.env <- new.env()
assign("init.hooks",c(),envir=lava.env)
assign("remove.hooks",c(),envir=lava.env)
assign("estimate.hooks",c(),envir=lava.env)
assign("color.hooks",c(),envir=lava.env)
assign("sim.hooks",c(),envir=lava.env)
assign("post.hooks",c(),envir=lava.env)
assign("print.hooks",c(),envir=lava.env)
assign("plot.post.hooks",c(),envir=lava.env)
assign("plot.hooks",c(),envir=lava.env)
assign("options", list(
                      trace=0,
                      tol=1e-6,
                      gamma=1,
                      backtrack="wolfe",
                      ngamma=0,
                      iter.max=300,
                      eval.max=250,
                      constrain=FALSE,
                      allow.negative.variance=FALSE,
                      progressbarstyle=3,
                      itol=1e-16,
                      cluster.index=versioncheck("mets",c(0,2,7)),
                      tobit=versioncheck("lava.tobit",c(0,5)),
                      Dmethod="simple", ##"Richardson"
                      messages=ifelse(interactive(), 1, 0),
                      parallel=TRUE,
                      param="relative",
                      sparse=FALSE,
                      test=TRUE,
                      coef.names=FALSE,
                      constrain=TRUE,
                      graph.proc="beautify",
                      regex=FALSE,
                      min.weight=1e-3,
                      exogenous=TRUE,
                      plot.engine="Rgraphviz",
                      node.color=c(exogenous="lightblue",endogenous="orange",
                                   latent="yellowgreen",transform="lightgray"),
                      edgecolor=FALSE,
                      layout="dot",
                      ## symbols=c("<-","<->"),
                      symbols=c("~","~~"),
                      devel=FALSE,
                      debug=FALSE), envir=lava.env)
