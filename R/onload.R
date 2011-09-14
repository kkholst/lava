
'.onLoad' <- function(lib, pkg="lava")
  {
    ##    addhook("heavytail.sim.hook","sim.hooks")
    addhook("heavytail.init.hook","init.hooks")
    addhook("cluster.post.hook","post.hooks")

    
    desc <- packageDescription(pkg)
    cat("Loading '", desc$Package, "' package...\n", sep="")
    cat("Version    : ", desc$Version, "\n", sep="")
    cat("Overview: help(package=", desc$Package, ")\n", sep="");
    ##lavalogo()
  }
