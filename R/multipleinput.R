##' Define risk difference association for binary exposure
##'
##' Set up model (risk difference) as defined in Richardson, Robins and Wang (2017).
##' @param x model
##' @param response response variable
##' @param exposure exposure variable
##' @param target.model design matrix for target model (relative difference)
##' @param nuisance.model design matrix for nuisance model (odds product)
##' @param exposure.model design matrix for propensity model
##' @param ... additional arguments to lower level functions
##' @export 
binomial.rd <- function(x,response,exposure,target.model,nuisance.model,exposure.model=binomial.lvm(),...) {    
    val <- list(list(input=c(exposure,target.model,nuisance.model),
                     fun=simulate.binomial.rd,
                     type="Binomial regression (exposure | risk-difference | odds-product)"))
    covariance(x,c(target.model,nuisance.model)) <- 0
    distribution(x,exposure) <- exposure.model
    names(val) <- response
    x$attributes$multiple.inputs <- val
    return(x)    
}

simulate.binomial.rd <- function(x,data,inputs,...) {
    exposure <- data[,inputs[1]]    
    lp1 <- data[,inputs[2]]
    lp2 <- data[,inputs[3]]
    rd <- tanh(lp1)
    op <- exp(lp2)
    pp <- RD_OP(rd,op)
    pp <- pp[,1]*(1-exposure) + pp[,2]*exposure
    y <- rbinom(NROW(data), 1, pp)
}


##' @export
Identical <- function(x,y=1,tolerance = .Machine$double.eps^0.5) {
    Mod(x-y)<tolerance
}

RD_OP <- function(rd,op) { # RD,OP
    a <- op-1
    b <- -op*(2-rd)-rd
    p0 <- (-b - sqrt(b^2 - 4*op*(1-rd)*a))/(2*a)
    op1 <- which(sapply(op,function(x) Identical(x,1)))
    if (length(op1)>0)   
        p0[op1] = 0.5 * (1 - rd[op1])
    p1 <- p0 + rd
    cbind(p0,p1)
}

simulate.multiple.inputs <- function(x,data,...) {
    minp <- x$attributes$multiple.inputs
    if (length(minp)>0) {
        for (i in seq_along(minp)) {            
            outcome <- names(minp[i])
            inp <- minp[[i]]$input
            fun <- minp[[i]]$fun
            data[,outcome] <- fun(x, data, inp)
        }
    }
    return(data)
}

addhook("simulate.multiple.inputs","sim.hooks")

printhook.multiple.inputs <- function(x,...) {
    minp <- x$attributes$multiple.inputs
    if (length(minp)>0) {
        outcomes <- names(minp)
        for (i in seq_along(minp)) {
            cat(minp[[i]]$type, ":\n\n")
            st <- paste0(outcomes[i]," ~ ", paste0(minp[[i]]$input,collapse=" | "))
            cat("  ", st, "\n")
            cat("\n")
        }
    }
    return(NULL)
}

addhook("printhook.multiple.inputs","print.hooks")
