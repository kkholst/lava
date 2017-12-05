
##' Weighted K-means
##'
##' Weighted K-means via Lloyd's algorithm
##' @param x Data (or formula)
##' @param mu Initial centers (or number centers chosen randomly among x)
##' @param data optional data frmae
##' @param weights Optional weights
##' @param iter.max Max number of iterations
##' @param n.start Number of restarts
##' @param ... Additional arguments to lower level functions
##' @export 
##' @author Klaus K. Holst
##' 
km <- function(x, mu, data, weights=rep(1,NROW(x)), iter.max=20, n.start=5, ...) { ## Lloyd's algorithm
    if (inherits(x, "formula")) x <- stats::model.matrix(x,data=data)
    x <- cbind(x)
    random.start <- TRUE
    if (is.list(mu)) {
        random.start <- FALSE
        n.start=1
        K <- length(mu)
    } else {
        K <- mu
    }
    sswmin <- Inf
    mus <- ssws <- clmin <- mumin <- NULL
    cl0 <- rep(1,NROW(x))
    for (k in seq(n.start)) {
        if (random.start) {
            mu <- lapply(sample(NROW(x),K), function(i) cbind(x)[i,,drop=TRUE])
        }
        mus <- c(mus, list(mu))
        for (i in seq(iter.max)) {
            d <- Reduce(cbind,lapply(mu, function(m) weights*colSums((t(x)-m)*(t(x)-m))))
            cl <- apply(d,1,which.min)
            for (j in seq_along(mu)) {
                idx <- which(cl==j)
                if (length(idx)) {
                    mu[[j]] <- colSums(cbind(apply(x[idx,,drop=FALSE],2,
                                            function(x) x*weights[idx])))/sum(weights[idx])
                }

            }
            if (sum(cl0-cl)==0L) break; # No change in assigment
        }
        ssw <- sum(d[cbind(seq(NROW(d)),cl)])
        ssws <- c(ssws,ssw)
        if (ssw < sswmin) {
            sswmin <- ssw
            clmin <- cl
            mumin <- mu
        }       
    }
    mu <- structure(mu,class="by",dim=K,dimnames=list(class=seq(K)))
    withinclusterss <- as.vector(by(d[cbind(seq(NROW(d)),cl)],cl,sum))
    return(list(cluster=cl,
           center=mu, ssw=withinclusterss)) ##, ssws,mu=Reduce(cbind,mus)))
}
