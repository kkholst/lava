`vars` <-
function(x,...) UseMethod("vars")

`vars.graph` <-
  function(x,...) {
    nodes(x)
  }

`vars.lvm` <-
  function(x,...) {
    nodes(Graph(x))
  }

`vars.lvmfit` <-
  function(x,...) {
    vars(Model(x),...)
  }

vars.list <- function(x,...) {
  varlist <- c()
  for (i in 1:length(x)) {
    varlist <- c(varlist, vars(x[[i]]))
  }
  varlist <- unique(varlist)
  return(varlist)
}
