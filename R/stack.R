##' Stack estimating equations
##'
##' Stack estimating equations
##' @param x Model 1
##' @param model2 Model 2 
##' @param D1u Derivative of score of model 2 w.r.t. parameter vector of model 1
##' @param inv.D2u Inverse of deri
##' @param weight weight (vector or function)
##' @param dweight derivative of weight wrt parameters of model 1
##' @param U Optional score function (model 2) as function of all parameters
##' @param k Debug argument
##' @param keep1 If TRUE only parameters of model 2 i s returned
##' @param ... Additional arguments to lower level functions
##' @aliases stack.estimate
##' @examples
##' 2+2
##' @export
stack.estimate <- function(x,model2,D1u,inv.D2u,weight,dweight,U,k=1,keep1=FALSE,...) {
    iid1 <- iid(x)
    iid2 <- iid(model2)
    if (missing(inv.D2u)) {
        inv.D2u <- -attributes(iid2)$bread
    }
    if (is.null(inv.D2u)) stop("Need derivative of second stage score")
    if (!missing(U)) {
        D1u <- numDeriv::jacobian(U,coef(x))
    }
    if (!missing(weight) && is.function(weight)) {
        dweight <- numDeriv::jacobian(weight,coef(x))
        weight <- weight(coef(x))
    }
    if (!missing(dweight)) {
        D2u <- Inverse(inv.D2u)
        u2 <- -iid2%*%D2u ## Score of stage two equation derived from estimated influence function
        ## Derivative of score wrt first set of parameters (weight-parameters)
        D1u <- crossprod(apply(u2,2,function(x) -x/weight),dweight)
    }
    ii <- iid(merge(x,model2))
    iid1. <- ii[,seq_along(coef(x)),drop=FALSE]
    iid2. <- ii[,length(coef(x))+seq_along(coef(model2)),drop=FALSE]
    iid3 <- t(inv.D2u%*%(D1u%*%t(iid1.)))
    if (!keep1) return(estimate(coef=coef(model2),iid=cbind(iid2.+k*iid3)))
    estimate(coef=c(coef(x),coef(model2)),iid=cbind(iid1.,iid2. + k*iid3))
}

##' @export
stack.glm <- function(x,model2,...) {
    stack(estimate(x),estimate(model2),...)
}


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
##' (e <- measurement.error(e1, z~1+x, data=d,predictfun=pp)[[1]])
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
    list(estimate=stacked, naive=e2, lm=coef(summary(stage.two)))
}
