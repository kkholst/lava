##' Conformal predicions
##'
##' @title Conformal prediction
##' @param object Model object (lm, glm or similar with predict method) or formula (lm)
##' @param data data.frame
##' @param newdata New data.frame to make predictions for
##' @param alpha Level of prediction interval
##' @param mad Conditional model (formula) for the MAD (locally-weighted CP)
##' @param ... Additional arguments to lower level functions
##' @return data.frame with fitted (fit), lower (lwr) and upper (upr) predictions bands.
##' @examples
##' set.seed(123)
##' n <- 200
##' x <- seq(0,6,length.out=n)
##' delta <- 3
##' ss <- exp(-1+1.5*cos((x-delta)))
##' ee <- rnorm(n,sd=ss)
##' y <- (x-delta)+3*cos(x+4.5-delta)+ee
##' d <- data.frame(y=y,x=x)
##'
##' newd <- data.frame(x=seq(0,6,length.out=50))
##' cc <- confpred(lm(y~splines::ns(x,knots=c(1,3,5)),data=d), data=d, newdata=newd)
##' if (interactive()) {
##' plot(y~x,pch=16,col=lava::Col("black"),ylim=c(-10,10),xlab="X",ylab="Y")
##' with(cc,
##'      lava::confband(newd$x,lwr,upr,fit,
##'         lwd=3,polygon=TRUE,col=Col("blue"),border=FALSE))
##' }
##' @export
confpred <- function(object,data,newdata=data,alpha=0.05,mad,...) { ## Split algorithm
    if (inherits(object,"formula")) {
        object <- do.call("lm",list(object,data=data,...))
    }
    dd <- csplit(data,0.5)
    muhat.new <- predict(object,newdata=newdata) ## New predictions
    muhat.1 <- predict(object,newdata=dd[[1]])      ## Training
    R1 <- abs(dd[[1]][,1]-muhat.1)
    muhat.2 <- predict(object,newdata=dd[[2]])   ## Ranking
    R2 <- abs(dd[[2]][,1]-muhat.2)
    if (missing(mad)) mad <- formula(object)
    if (is.null(mad)) { 
        mad.new <- rep(1,nrow(newdata))
    } else { ## Locally-weighted conformal ffinference
        if (names(dd[[2]])[1] %ni% names(newdata)) {
            newdata <- cbind(0,newdata); names(newdata)[1] <- names(dd[[2]])[1]
        }
        X0 <- model.matrix(mad,data=newdata)
        if (inherits(mad,"formula")) { 
            X2 <- model.matrix(mad,dd[[2]])            
            mad.obj <- stats::lm.fit(x=X2,y=R2)
            mad2 <- X2%*%mad.obj$coefficients
            mad.new <- X0%*%mad.obj$coefficients
        } else {
            mad.obj <- do.call(mad,list(y=R2,x=dd[[2]]))
            mad2 <- predict(mad.obj,newdata=dd[[2]])
            mad.new <- predict(mad.obj,newdata=newdata)
        }
        R2 <- R2/mad2
    }
    k <- ceiling((nrow(data)/2+1)*(1-alpha))
    if (k==0) k <- 1
    if (k>length(R2)) k <- length(R2)
    q <- sort(R2)[k] ## 1-alpha quantile
    lo <- muhat.new - q*mad.new
    up <- muhat.new + q*mad.new
    data.frame(fit=muhat.new,lwr=lo,upr=up)
}
