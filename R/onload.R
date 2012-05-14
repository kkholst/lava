'.onAttach' <- function(lib, pkg="lava")
  {
    options(lava=list(
              init.hooks=c(),
              estimate.hooks=c(),
              color.hooks=c(),
              sim.hooks=c(),
              post.hooks=c(),
              options=list(
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
            )))
    
    addhook("heavytail.init.hook","init.hooks")
    addhook("glm.estimate.hook","estimate.hooks")
    addhook("cluster.post.hook","post.hooks")
    
    desc <- packageDescription(pkg)
    packageStartupMessage("\nLoading '", desc$Package, "' package...\n",
                          "Version    : ", desc$Version, "\n",
                          "Overview: help(package=", desc$Package, ")\n");
  }
