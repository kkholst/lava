##' @export
predict.lvm.mixture <- function(object,p=coef(object,full=TRUE),model="normal",predict.fun=NULL,...) {
    p0 <- coef(object,full=FALSE)
    pp <- p[seq_along(p0)]
    pr <- p[length(p0)+seq(length(p)-length(p0))];
    if (length(pr)<object$k) pr <- c(pr,1-sum(pr))
    myp <- modelPar(object$multigroup,p=pp)$p
    logff <- sapply(seq(object$k), function(j) (logLik(object$multigroup$lvm[[j]],p=myp[[j]],data=object$data,indiv=TRUE,model=model)))
    logplogff <- t(apply(logff,1, function(y) y+log(pr)))
    zmax <- apply(logplogff,1,max)
    logsumpff <- log(rowSums(exp(logplogff-zmax)))+zmax
    aji <- apply(logplogff,2,function(x) exp(x-logsumpff))
    gamma <- exp(apply(logplogff,2,function(y) y - logsumpff)) ## Posterior class probabilities, conditional mean
    M <- 0; V <- 0
    x <- lava::vars(object$model)
    for (i in seq(object$k)) {
        m <- Model(object$multigroup)[[i]]
        P <- predict(m,data=object$data,p=myp[[i]],x=x)
        if (!is.null(predict.fun)) {
            M <- M+gamma[,i]*predict.fun(P,as.vector(attributes(P)$cond.var), ...)
        } else {
            M <- M+gamma[,i]*P
            V <- V+gamma[,i]*as.vector(attributes(P)$cond.var) + gamma[,i]*P^2
        }
    }
    ## Var(X) = E[(X-m)^2] = EX^2 - m^2 
    ##        = sum(pi*Xi^2) - [sum(pi*mi)]^2 = sum(pi*sigmai^2) + sum(pi*mi^2) - [sum(pi*mi)]^2
    if (!is.null(predict.fun)) return(M)
    V <- V-M^2    
    structure(M,cond.var=V)
}
