##' @export
stack.estimate <- function(x,y,Ix,iIy,weight,dweight,...) {    
    iid1 <- iid(x)
    iid2 <- iid(y)
    if (missing(iIy)) {
        iIy <- attributes(iid2)$bread
    }
    if (!missing(dweight)) {
        u2 <- apply(Inverse(iIy)%*%t(iid2),1,function(x) x/weight)
        Ix <- colSums(
            rbind(rep(1,ncol(dweight)))%x%u2*
            rbind(rep(1,ncol(u2)))%x%dweight
            )
        Ix <- matrix(Ix,byrow=TRUE,ncol=ncol(iid1))
    }
    ii <- iid(merge(x,y))
    iid1. <- ii[,seq_along(coef(x)),drop=FALSE]
    iid2. <- ii[,length(coef(x))+seq_along(coef(y)),drop=FALSE] 
    iid3 <- t(iIy%*%(Ix%*%t(iid1.)))
    estimate(coef=c(coef(x),coef(y)),iid=cbind(iid1.,iid2.-iid3))
}
