lava.env <- environment()
assign("estimate.hooks",c(),lava.env)
lockBinding("estimate.hooks", lava.env)
assign("color.hooks",c(),lava.env)
lockBinding("color.hooks", lava.env)

gethook <- function(hook="estimate.hooks",...) {
  get(hook,envir=lava.env)
}
addhook <- function(x,hook="estimate.hooks",...) {
  ##  env <- as.environment("package:lava")
  unlockBinding(hook,lava.env)
  assign(hook,c(gethook(hook),x), envir=lava.env)
  lockBinding(hook, lava.env)
}
