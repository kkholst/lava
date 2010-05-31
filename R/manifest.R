
`manifest` <-
function(x,...) UseMethod("manifest")

`manifest.lvm` <-
function(x,...) {
  if (length(vars(x))>0) 
    setdiff(vars(x),latent(x))
  else
    NULL
}
`manifest.lvmfit` <-
function(x,...) {
  manifest(Model(x))
}

manifest.list <- function(x,...) {
  manifestlist <- c()
  for (i in 1:length(x)) {
    manifestlist <- c(manifestlist, manifest(x[[i]]))
  }
  endolist <- unique(manifestlist)
  return(manifestlist)
}

`manifest.multigroup` <-
function(x,...) {
  manifest(Model(x))
}

