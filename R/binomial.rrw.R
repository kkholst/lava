##' Define constant risk difference or relative risk association for binary exposure
##'
##' Set up model as defined in Richardson, Robins and Wang (2017).
##' @param x model
##' @param response response variable (character or formula)
##' @param exposure exposure variable (character or formula)
##' @param target.model variable defining the linear predictor for the target model
##' @param nuisance.model variable defining the linear predictor for the nuisance model
##' @param exposure.model model for exposure (default binomial logit link)
##' @param ... additional arguments to lower level functions
##' @aliases binomial.rd binomial.rr
##' @export
binomial.rd <- function(x,response,exposure,
                 target.model,nuisance.model,
                 exposure.model=binomial.lvm(),...) {
    binomial.rrw(x,response,exposure,target.model,nuisance.model,exposure.model,type="rd",...)
}

##' @export
binomial.rr <- function(x,response,exposure,
                 target.model,nuisance.model,
                 exposure.model=binomial.lvm(),...) {
    binomial.rrw(x,response,exposure,target.model,nuisance.model,exposure.model,type="rr",...)
}


binomial.rrw <- function(x, response, exposure,
                 target.model, nuisance.model,
                 exposure.model=binomial.lvm(), type="rd", ...) {
    if (inherits(response,"formula")) {
        vars <- all.vars(response)
        if (length(vars)==1L) {
            response <- vars
        } else {
            yf <- getoutcome(response, sep="|")
            exposure  <- attr(yf,"x")[[1]]
            if (length(attr(yf,"x"))>1)
                target.model  <- attr(yf,"x")[[2]]
            if (length(attr(yf,"x"))>2)
                nuisance.model  <- attr(yf,"x")[[3]]
            response <- yf[1]
        }
    }
    if (inherits(exposure,"formula")) exposure <- all.vars(exposure)
    if (inherits(target.model,"formula")) target.model <- all.vars(target.model)
    if (inherits(nuisance.model,"formula")) nuisance.model <- all.vars(nuisance.model)

    if (type=="rd") {
        val <- list(list(input=c(exposure,target.model,nuisance.model),
                         fun=simulate_binomial_rd,
                         type="Binomial regression (exposure | risk-difference | odds-product)"))
    } else {
        val <- list(list(input=c(exposure,target.model,nuisance.model),
                         fun=simulate_binomial_rr,
                         type="Binomial regression (exposure | relative-risk | odds-product)"))
    }
    if (is.null(distribution(x)[[exposure]]))
        distribution(x, exposure) <- binomial.lvm(link="logit")
    covariance(x,c(target.model,nuisance.model)) <- 0
    distribution(x,exposure) <- exposure.model
    names(val) <- response
    x$attributes$multiple.inputs <- val
    return(x)
}


simulate_binomial_rd <- function(x,data,inputs,...) {
    exposure <- data[,inputs[1]]
    lp1 <- data[,inputs[2]]
    lp2 <- data[,inputs[3]]
    rd <- tanh(lp1)
    op <- exp(lp2)
    pp <- RD_OP(rd,op)
    pp <- pp[,1]*(1-exposure) + pp[,2]*exposure
    y <- rbinom(NROW(data), 1, pp)
}

simulate_binomial_rr <- function(x,data,inputs,...) {
    exposure <- data[,inputs[1]]
    lp1 <- data[,inputs[2]]
    lp2 <- data[,inputs[3]]
    rr <- exp(lp1)
    op <- exp(lp2)
    pp <- RR_OP(rr,op)
    pp <- pp[,1]*(1-exposure) + pp[,2]*exposure
    y <- rbinom(NROW(data), 1, pp)
}


##' @export
Identical <- function(x,y=1,tolerance = .Machine$double.eps^0.5) {
    Mod(x-y)<tolerance
}

RD_OP <- function(rd,op) {
    a <- op-1
    b <- -op*(2-rd)-rd
    p0 <- (-b - sqrt(b^2 - 4*op*(1-rd)*a))/(2*a)
    op1 <- which(sapply(op,function(x) Identical(x,1)))
    if (length(op1)>0)
        p0[op1] = 0.5 * (1 - rd[op1])
    p1 <- p0 + rd
    cbind(p0,p1)
}


RR_OP <- function(rr,op) {
    b <- op*(1+rr)
    a <- rr*(1-op)
    p0 <- (-b + sqrt(b^2 + 4*a*op))/(2*a)
    op1 <- which(sapply(op,function(x) Identical(x,1)))
    if (length(op1)>0)
        p0[op1] = 1/(1+rr)
    p1 <- p0*rr
    cbind(p0,p1)
}
