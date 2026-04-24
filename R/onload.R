'.onLoad' <- function(libname, pkgname="lava") {
    addhook("heavytail_init_hook","init.hooks")
    addhook("glm_estimate_hook","estimate.hooks")
    addhook("ordinal_estimate_hook","estimate.hooks")
    addhook("cluster_post_hook","post.hooks")
    addhook("ordinal_sim_hook","sim.hooks")
    addhook("color_ordinal","color.hooks")
    addhook("ordinal_remove_hook","remove.hooks")
    lava.options(cluster.index = packagecheck("mets"))
}

'.onAttach' <- function(libname, pkgname="lava") {
    # desc <- utils::packageDescription(pkgname)
    # packageStartupMessage(desc$Package, " version ",desc$Version)
}
