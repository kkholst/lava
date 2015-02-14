'.onAttach' <- function(libname, pkgname="lava")
  {
    addhook("heavytail.init.hook","init.hooks")
    addhook("glm.estimate.hook","estimate.hooks")
    addhook("cluster.post.hook","post.hooks")
    addhook("ordinal.sim.hook","sim.hooks")
    addhook("color.ordinal","color.hooks")

    desc <- utils::packageDescription(pkgname)
    ## packageStartupMessage("Loading '", desc$Package, "' package...\n",
    ##                       "\tVersion: ", desc$Version, "\n",
    ##                       "\tOverview: help(package=", desc$Package, ")");
    packageStartupMessage(desc$Package, " version ",desc$Version);

  }
