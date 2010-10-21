lava.env <- environment()
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
            constrain=FALSE,
            silent=FALSE,
            Dmethod="simple" ##Richardson"
            ),
       lava.env)
lockBinding("post.hooks", lava.env)

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
gethook <- function(hook="estimate.hooks",...) {
  get(hook,envir=lava.env)
}
addhook <- function(x,hook="estimate.hooks",...) {
  ##  env <- as.environment("package:lava")
  unlockBinding(hook,lava.env)
  assign(hook,c(gethook(hook),x), envir=lava.env)
  lockBinding(hook, lava.env)
}
