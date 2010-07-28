
'.onLoad' <- function(lib, pkg="lava")
  {
    desc <- packageDescription(pkg)
    cat("Loading '", desc$Package, "' package...\n", sep="")
    cat("Version    : ", desc$Version, "\n", sep="")
    cat("Overview: help(package=", desc$Package, ")\n", sep="");
    ##lavalogo()
  }
