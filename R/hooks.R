##' Set global options for \code{lava}
##' 
##' Extract and set global parameters of \code{lava}. In particular optimization
##' parameters for the \code{estimate} function.
##' 
##' \itemize{ \item \code{param}: 'relative' (factor loading and variance of one
##' endogenous variables in each measurement model are fixed to one), 'absolute'
##' (mean and variance of latent variables are set to 0 and 1, respectively),
##' 'hybrid' (intercept of latent variables is fixed to 0, and factor loading of
##' at least one endogenous variable in each measurement model is fixed to 1),
##' 'none' (no constraints are added) \item \code{silent}: Set to \code{FALSE}
##' to disable various output messages \item ...  } see \code{control} parameter
##' of the \code{estimate} function.
##' 
##' @param \dots Arguments
##' @return \code{list} of parameters
##' @author Klaus K. Holst
##' @keywords models
##' @examples
##' 
##' \dontrun{
##' lava.options(iter.max=100,silent=TRUE)
##' }
##'
##' @export
lava.options <- function(...) {
  dots <- list(...)
  curopt <- get("options",envir=lava.env)
  if (length(dots)==0) 
    return(curopt)
  unlockBinding("options",lava.env)
  idx <- which(names(dots)!="")
  curopt[names(dots)[idx]] <- dots[idx]
  assign("options",curopt, envir=lava.env)  
  lockBinding("options", lava.env)  
}

##' @export
gethook <- function(hook="estimate.hooks",...) {
  get(hook,envir=lava.env)
}

##' @export
addhook <- function(x,hook="estimate.hooks",...) {
  ##  env <- as.environment("package:lava")
  unlockBinding(hook,lava.env)
  assign(hook,unique(c(gethook(hook),x)), envir=lava.env)
  lockBinding(hook, lava.env)
}

lava.env <- environment()
assign("init.hooks",c(),lava.env)
lockBinding("init.hooks", lava.env)
assign("estimate.hooks",c(),lava.env)
lockBinding("estimate.hooks", lava.env)
assign("color.hooks",c(),lava.env)
lockBinding("color.hooks", lava.env)
assign("sim.hooks",c(),lava.env)
lockBinding("sim.hooks", lava.env)
assign("post.hooks",c(),lava.env)
lockBinding("post.hooks", lava.env)
assign("options",
       list(
            trace=0,
            iter.max=300,
            eval.max=250,
            constrain=FALSE,
            silent=TRUE,            
            itol=1e-9,
            Dmethod="simple", ##Richardson"
            parallel=TRUE,
            param="relative",
            constrain=TRUE,
            exogenous=TRUE,
            Rgraphviz=TRUE,
            debug=FALSE
            ),
       lava.env)

