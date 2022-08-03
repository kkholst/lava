##################################################
## Cohen's kappa
##################################################

##' @export
kappa.multinomial <- function(z,all=FALSE,...) {
    pp <- length(coef(z))
    if ((length(z$levels)!=2) || !(identical(z$levels[[1]],z$levels[[2]])))
        stop("Expected square table and same factor levels in rows and columns")
    k <- length(z$levels[[1]])
    zeros <- rbind(rep(0,pp))
    A0 <- zeros; A0[diag(z$position)] <- 1
    A <- matrix(0,ncol=pp,nrow=2*k)
    for (i in seq(k)) A[i,z$position[i,]] <- 1
    for (i in seq(k)) A[i+k,z$position[,i]] <- 1
    b <- estimate(z,function(p) as.vector(rbind(A0,A)%*%p),IC=TRUE)
    b2 <- estimate(b,function(p) c(p[1],sum(p[seq(k)+1]*p[seq(k)+k+1])),IC=TRUE)
    if (!all) {
        return(estimate(b2,function(p) list(kappa=(p[1]-p[2])/(1-p[2])),IC=TRUE,...))
    }
    estimate(b2,function(p) list(kappa=(p[1]-p[2])/(1-p[2]),agree=p[1], independence=p[2]),IC=TRUE,...)
}

##' @export
kappa.table <- function(z,...) {
    kappa(multinomial(Expand(z)),...)
}

##' @export
kappa.data.frame <- function(z,...) {
    kappa(multinomial(z),...)
}
