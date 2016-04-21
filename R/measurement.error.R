##' Two-stage (non-linear) measurement error
##'
##' Two-stage measurement error
##' @param model1 Stage 1 model
##' @param formula Formula specifying observed covariates in stage 2 model
##' @param data data.frame
##' @param predictfun Predictions to be used in stage 2
##' @param id1 Optional id-vector of stage 1
##' @param id2 Optional id-vector of stage 2
##' @param ... Additional arguments to lower level functions
##' @seealso stack.estimate
##' @export
##' @examples
##' m <- lvm(c(y1,y2,y3)~u,c(y3,y4,y5)~v,u~~v,c(u,v)~x)
##' transform(m,u2~u) <- function(x) x^2
##' transform(m,uv~u+v) <- prod
##' regression(m) <- z~u2+u+v+uv+x
##' set.seed(1)
##' d <- sim(m,1000,p=c("u,u"=1))
##' 
##' ## Stage 1
##' m1 <- lvm(c(y1[0:s],y2[0:s],y3[0:s])~1*u,c(y3[0:s],y4[0:s],y5[0:s])~1*v,u~b*x,u~~v)
##' latent(m1) <- ~u+v
##' e1 <- estimate(m1,d)
##' 
##' pp <- function(mu,var,data,...) {
##'     cbind(u=mu[,"u"],u2=mu[,"u"]^2+var["u","u"],v=mu[,"v"],uv=mu[,"u"]*mu[,"v"]+var["u","v"])
##' }
##' (e <- measurement.error(e1, z~1+x, data=d, predictfun=pp))
##' 
##' ## uu <- seq(-1,1,length.out=100)
##' ## pp <- estimate(e,function(p,...) p["(Intercept)"]+p["u"]*uu+p["u2"]*uu^2)$coefmat
##' if (interactive()) {
##'     plot(e,intercept=TRUE,vline=0)
##' 
##'     f <- function(p) p[1]+p["u"]*u+p["u2"]*u^2
##'     u <- seq(-1,1,length.out=100)
##'     plot(e, f, data=data.frame(u), ylim=c(-.5,2.5))
##' }
measurement.error <- function(model1, formula, data=parent.frame(), predictfun=function(mu,var,data,...) mu[,1]^2+var[1], id1, id2, ...) {
    if (!inherits(model1,c("lvmfit","lvm.mixture"))) stop("Expected lava object ('lvmfit','lvm.mixture',...)")
    if (missing(formula)) stop("formula needed for stage two (right-hand side additional covariates)")
    p1 <- coef(model1,full=TRUE)
    uhat <- function(p=p1) {
        P <- predict(model1,p=p,x=manifest(model1))
        cbind(predictfun(P,attributes(P)$cond.var,data))
    }
    if (missing(id1)) id1 <- seq(nrow(model.frame(model1)))
    if (missing(id2)) id2 <- seq(nrow(model.frame(model1)))
    if (!inherits(model1,"estimate"))
        e1 <- estimate(NULL,coef=p1,id=id1,iid=iid(model1))
    u <- uhat()
    X0 <- model.matrix(formula, data)
    Y <- model.frame(formula,data)[,1]
    X <- cbind(X0,u)
    stage.two <- lm(Y~-1+X)
    names(stage.two$coefficients) <- colnames(X)
    if (!inherits(stage.two,"estimate"))
        e2 <- estimate(stage.two, id=id2)
    U <- function(alpha=p1,beta=coef(stage.two)) {
        X <- cbind(X0,uhat(alpha))
        r <- (Y-X%*%beta)/summary(stage.two)$sigma^2
        apply(X,2,function(x) sum(x*r))
    }
    Ia <- -numDeriv::jacobian(function(p) U(p),p1)
    stacked <- stack(e1,e2,Ia)
    res <- c(stacked,list(naive=e2,lm=coef(summary(stage.two)),fun=predictfun))
    ## res <- list(estimate=stacked, naive=e2, lm=coef(summary(stage.two)),
    ##             fun=predictfun)
    structure(res,class=c("measurement.error","estimate"))
}



##' Two-stage estimator (non-linear SEM)
##'
##' Two-stage estimator for non-linear structural equation models
##' @export
##' @param model1 Stage 1 measurement model
##' @param model2 Stage 2 SEM
##' @param data data.frame
##' @param predictfun Prediction of latent variable
##' @param id1 Optional id-variable (stage 1 model) 
##' @param id2 Optional id-variable (stage 2 model)
##' @param ... Additional arguments to lower level functions
##' @examples
##' m <- lvm(c(x1,x2,x3)~f1,f1~z,
##'          c(y1,y2,y3)~f2,f2~f1+z)
##' latent(m) <- ~f1+f2
##' d <- simulate(m,100,p=c("f2,f2"=2,"f1,f1"=0.5),seed=1)
##' 
##' ## Full MLE
##' ee <- estimate(m,d)
##' 
##' ## Manual two-stage
##' \dontrun{
##' m1 <- lvm(c(x1,x2,x3)~f1,f1~z); latent(m1) <- ~f1
##' e1 <- estimate(m1,d)
##' pp1 <- predict(e1,f1~x1+x2+x3)
##' 
##' d$u1 <- pp1[,]
##' d$u2 <- pp1[,]^2+attr(pp1,"cond.var")
##' m2 <- lvm(c(y1,y2,y3)~eta,c(y1,eta)~u1+u2+z); latent(m2) <- ~eta
##' e2 <- estimate(m2,d)
##' }
##' 
##' ## Two-stage
##' m1 <- lvm(c(x1,x2,x3)~f1,f1~z); latent(m1) <- ~f1
##' m2 <- lvm(c(y1,y2,y3)~eta,c(y1,eta)~u1+u2+z); latent(m2) <- ~eta
##' pred <- function(mu,var,data,...)
##'     cbind("u1"=mu[,1],"u2"=mu[,1]^2+var[1])
##' (mm <- twostage(m1,m2,data=d,predictfun=pred))
##' 
##' if (interactive()) {
##'     pf <- function(p) p["eta"]+p["eta~u1"]*u + p["eta~u2"]*u^2
##'     plot(mm,f=pf,data=data.frame(u=seq(-2,2,length.out=100)),lwd=2)
##' }
twostage <- function(model1, model2, data=parent.frame(),
                         predictfun=function(mu,var,data,...)
                             cbind("u1"=mu[,1],"u2"=mu[,1]^2+var[1]),
                     id1,id2, ...) {
    if (inherits(model1,"lvm")) {
        model1 <- estimate(model1,data=data,...)
    }    
    if (!inherits(model1,c("estimate","lvmfit","lvm.mixture"))) stop("Expected lava object ('estimate','lvmfit','lvm.mixture',...)")
    if (!inherits(model2,c("lvm"))) stop("Expected lava object ('lvm',...)")
    p1 <- coef(model1)
    uhat <- function(p=p1) {
        P <- predictlvm(model1,p=p,data=model.frame(model1))
        cbind(predictfun(P$mean,P$var,model.frame(model1)))
    }
    pp <- uhat()
    newd <- data
    newd[,colnames(pp)] <- pp
    model2 <- estimate(model2,data=newd,...)
    p2 <- coef(model2)    
    if (missing(id1)) id1 <- seq(nrow(model.frame(model1)))
    if (missing(id2)) id2 <- seq(nrow(model.frame(model2)))
    if (!inherits(model1,"estimate")) {
        e1 <- estimate(NULL,coef=p1,id=id1,iid=iid(model1))
    }
    e2 <- estimate(model2, id=id2)
    U <- function(alpha=p1,beta=p2) {
        pp <- uhat(alpha)
        newd <- model.frame(model2)
        newd[,colnames(pp)] <- pp
        score(model2,p=beta,data=newd)
    }
    Ia <- -numDeriv::jacobian(function(p) U(p),p1)
    stacked <- stack(e1,e2,Ia)

    coef <- model2$coef    
    res <- model2
    res[names(stacked)] <- stacked
    res$estimator <- "generic"
    cc <- stacked$coefmat[,c(1,2)];
    cc <- cbind(cc,cc[,1]/cc[,2],stacked$coefmat[,5])
    coef[,] <- cc
    res$coef <- coef
    res$fun <- predictfun
        
    structure(res,class=c("measurement.error","lvmfit","estimate"))
}

