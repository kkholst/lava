##' @export
Inverse <- function(X,tol=lava.options()$itol,det=TRUE,names=!chol,chol=FALSE,symmetric=FALSE) {
    n <- NROW(X)
    if (n==1L) {
        res <- 1/X
        if (det) attributes(res)$det <- X
        if (chol) attributes(res)$chol <- X
        return(res)
    }
    if (chol) {
        L <- chol(X)
        res <- chol2inv(L)
        if (det) attributes(res)$det <- prod(diag(L)^2)
        if (chol) attributes(res)$chol <- X
    } else {
        if(symmetric){
            decomp <- eigen(X, symmetric = TRUE)
            D <- decomp$values
            U <- decomp$vectors
            V <- decomp$vectors
        }else{
            X.svd <- svd(X)
            U <- X.svd$u
            V <- X.svd$v
            D <- X.svd$d
        }
        id0 <- numeric(n)
        idx <- which(abs(D)>tol)
        id0[idx] <- 1/D[idx]
        res <- V%*%diag(id0,nrow=length(id0))%*%t(U)

        if (det)
            attributes(res)$det <- prod(D[D>tol])
        attributes(res)$pseudo <- (length(idx)<n)
        attributes(res)$minSV <- min(D)

    }
    if (names && !is.null(colnames(X))) dimnames(res) <- list(colnames(X),colnames(X))
    return(res)
}
