##' @export
`measurement` <-
    function(x, ...) {
        M <- x$M
        latent.idx <- match(latent(x),vars(x))
        obs.idx <- match(manifest(x),vars(x))
        if (length(latent.idx)==0)
            return(NULL)

        measurementmodels <- c()
        for (i in seq_along(latent.idx)) {
            ii <- latent.idx[i]

            relation <- M[ii,obs.idx]==1
            byNodes <- names(relation)[relation]
            newnodes <- c(latent(x)[i],byNodes)
            lvm1 <- subset(x,newnodes)
            measurementmodels <- c(measurementmodels, list(lvm1))
        }

        measurementmodels
    }
