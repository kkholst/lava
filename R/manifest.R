##' @export
`manifest` <-
function(x,...) UseMethod("manifest")

##' @S3method manifest lvm
`manifest.lvm` <-
function(x,...) {
  if (length(vars(x))>0) 
    setdiff(vars(x),latent(x))
  else
    NULL
}

##' @S3method manifest lvmfit
`manifest.lvmfit` <-
function(x,...) {
  manifest(Model(x))
}

##' @S3method manifest list
manifest.list <- function(x,...) {
  manifestlist <- c()
  for (i in 1:length(x)) {
    manifestlist <- c(manifestlist, manifest(x[[i]]))
  }
  endolist <- unique(manifestlist)
  return(manifestlist)
}

##' @S3method manifest multigroup
`manifest.multigroup` <-
function(x,...) {
  manifest(Model(x))
}

