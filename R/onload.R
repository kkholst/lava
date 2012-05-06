'.onAttach' <- function(lib, pkg="lava")
  {
    addhook("heavytail.init.hook","init.hooks")
    addhook("glm.estimate.hook","estimate.hooks")
    addhook("cluster.post.hook","post.hooks")
    
    desc <- packageDescription(pkg)
    packageStartupMessage("\nLoading '", desc$Package, "' package...\n",
                          "Version    : ", desc$Version, "\n",
                          "Overview: help(package=", desc$Package, ")\n");
  }
