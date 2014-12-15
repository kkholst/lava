##' Stack estimating equations
##'
##' Stack estimating equations
##' @param model1 Model 1
##' @param model2 Model 2 
##' @param D1u Derivative of score of model 2 w.r.t. parameter vector of model 1
##' @param inv.D2u Inverse of deri
##' @param weight weight (vector or function)
##' @param dweight derivative of weight wrt parameters of model 1
##' @param U Optional score function (model 2) as function of all parameters
##' @param k Debug argument
##' @param keep.x If TRUE only parameters of model 2 i s returned
##' @param ... Additional arguments to lower level functions
##' @aliases stack.estimate
##' @examples
##' 2+2
##' @export
stack.estimate <- function(model1,model2,D1u,inv.D2u,weight,dweight,U,k=1,keep1=FALSE,...) {
    iid1 <- iid(model1)
    iid2 <- iid(model2)
    if (missing(inv.D2u)) {
        inv.D2u <- -attributes(iid2)$bread
    }
    if (is.null(inv.D2u)) stop("Need derivative of second stage score")
    if (!missing(U)) {
        D1u <- numDeriv::jacobian(U,coef(model1))
    }
    if (!missing(weight) && is.function(weight)) {
        dweight <- numDeriv::jacobian(weight,coef(model1))
        weight <- weight(coef(model1))
    }
    if (!missing(dweight)) {
        D2u <- Inverse(inv.D2u)
        u2 <- -iid2%*%D2u ## Score of stage two equation derived from estimated influence function
        ## Derivative of score wrt first set of parameters (weight-parameters)
        D1u <- crossprod(apply(u2,2,function(x) -x/weight),dweight)
    }
    ii <- iid(merge(model1,model2))
    iid1. <- ii[,seq_along(coef(model1)),drop=FALSE]
    iid2. <- ii[,length(coef(model1))+seq_along(coef(model2)),drop=FALSE]
    iid3 <- t(inv.D2u%*%(D1u%*%t(iid1.)))
    if (!keep1) return(estimate(coef=coef(model2),iid=cbind(iid2.+k*iid3)))
    estimate(coef=c(coef(model),coef(model2)),iid=cbind(iid1.,iid2. + k*iid3))
}

##' @export
stack.glm <- function(model1,model2,...) {
    stack(estimate(model1),estimate(model2),...)
}


##' @export
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
