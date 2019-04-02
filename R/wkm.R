## kmeans++
## Assign centre of first cluster randomly from observations.
## 1) Calculate min distances of observations to current clusters, D
## 2) Sample new cluster centre from distribution p(x) = D^2(x)/sum(D^2)
## Repeat 1-2 until all clusters are assigned
kmpp <- function(y,k=2) {
    Dist <- function(y1,y2) sum(y1-y2)^2
    n <- NROW(y)
    ii <- numeric(k)
    u <- runif(k)
    ii[1] <- sample(n,1)
    D <- matrix(0,n,k-1)
    for (i in seq_len(k-1)+1) { 
        D[,i-1] <- apply(y,1,function(x) Dist(x,y[ii[i-1],]))
        D2 <- apply(D[,seq(i-1),drop=FALSE],1,min)
        pdist <- cumsum(D2/sum(D2))
        ii[i] <- mets::fast.approx(pdist,u[i])
    }
    return(ii)
}


##' Weighted K-means
##'
##' Weighted K-means via Lloyd's algorithm
##' @param x Data (or formula)
##' @param mu Initial centers (or number centers chosen randomly among x)
##' @param data optional data frmae
##' @param weights Optional weights
##' @param iter.max Max number of iterations
##' @param n.start Number of restarts
##' @param init method to create initial centres (default kmeans++)
##' @param ... Additional arguments to lower level functions
##' @export
##' @author Klaus K. Holst
##'
wkm <- function(x, mu, data, weights=rep(1,NROW(x)), iter.max=20, n.start=5, init="kmpp", ...) { ## Lloyd's algorithm
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
    mus <- ssws <- NULL
    cl0 <- rep(1,NROW(x))
    for (k in seq(n.start)) {
        if (random.start) {
            if (!exists(init)) { ## Random select centres
                idx <- sample(NROW(x),K)
            } else {
                idx <- do.call(init, list(x, K))
            }            
            mu <- lapply(idx, function(i) cbind(x)[i,,drop=TRUE])
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
           center=mu, ssw=withinclusterss))
}
