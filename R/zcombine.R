
excoef <- function(x,digits=2,p.digits=3,format=FALSE,fun,...) {
    cc <- coef(summary(x))
    res <- round(cbind(cc[,1:3,drop=FALSE],confint(x)),max(1,digits))
    pvalround <- round(cc[,4], max(1, p.digits))
    if (format) {
        res <- base::format(res,digits=digits)
        pval <- format(pvalround,p.digits=p.digits)
    } else {
        ## res <- format(res)
        pval <- format(pvalround)
    }
    pval <- paste("p=",pvalround,sep="")
    pval[which(pvalround<10^(-p.digits))] <- paste("p<0.",paste(rep("0",p.digits-1),collapse=""),"1",sep="")
    res <- cbind(res,pval)
    res2 <- apply(res,1,function(x)
                  paste(x[1]," [", x[4], ";",x[5],"], ",x[6],sep=""))
    names(res2) <- names(coef(x))
    if (!missing(fun)) {
        res2 <- c(res2,fun(x))
    }         
    res2    
}

##' Report estimates across different models
##'
##' @title Report estimates across different models
##' @param x list of model objects
##' @param ... additional arguments to lower level functions
##' @author Klaus Kähler Holst
##' @examples
##' data(serotonin)
##' m1 <- lm(cau ~ age*gene1 + age*gene2,data=serotonin)
##' m2 <- lm(cau ~ age + gene1,data=serotonin)
##' m3 <- lm(cau ~ age*gene2,data=serotonin)
##' 
##' Combine(list(A=m1,B=m2,C=m3),fun=function(x) c(R2=format(summary(x)$r.squared,digits=2)))
##' @export
Combine <- function(x,...) {
    ll <- lapply(x,excoef,...)
    nn <- lapply(ll,names)
    n0 <- unique(unlist(nn,use.names=FALSE))
    res <- matrix(NA,ncol=length(ll),nrow=length(n0))
    colnames(res) <- seq(length(ll))
    rownames(res) <- n0
    for (i in seq(length(ll))) {
        res[match(names(ll[[i]]),n0),i] <- ll[[i]]
    }
    colnames(res) <- names(ll)
    return(res)
}

