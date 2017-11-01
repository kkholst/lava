##' @export
estimate.formula <- function(x,family=stats::gaussian, data, weights, ...,
                     model="glm", lvm=FALSE) {
    cl <- match.call()
    if (lvm) {
        cl[[1]] <- quote(estimate0)
        return(eval(cl,envir=parent.frame()))
    }
    if (missing(data)) {
        data <- environment(x)
        cl$data <- quote(data)
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("x", "data", "weights", "subset",
                "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf$na.action <- na.pass0
    names(mf)[which(names(mf)=="x")] <- "formula"
    mf <- eval(mf, parent.frame())
    idx <- attr(na.omit(mf),"na.action")
    weights <- as.vector(model.weights(mf))
    if (length(idx)>0) {
        if (is.null(weights)) weights <- rep(1,nrow(mf))
        if (length(idx)>0) weights[idx] <- 0
        cl$weights <- quote(weights)
        cl$na.action <- quote(na.pass0)
    }
    argsModelObj <- names(formals(model))
    dots <- list(...)    
    idx <- which(names(dots) %ni% argsModelObj)
    rmarg <- c("model","raw","lvm")
    if (length(idx)>0) rmarg <- c(names(dots)[idx],rmarg)
    cl[rmarg] <- NULL
    names(cl)[names(cl)=="x"] <- "formula"
    cl[[1]] <- as.name(model)
    fit <- eval(cl)
    if (length(idx)==0) return(fit)
    optarg <- NULL
    if (length(idx)>0) {
        optarg <- dots[idx]
    }    
    do.call(estimate, c(list(fit),optarg))
}


na.pass0 <- function(object,...) {
    ## Fill in "zeros" in the design matrix where we have missing data
    df <- na.omit(object) 
    idx <- attr(df,"na.action")
    for (i in seq_len(ncol(object))) {
        if (is.logical(object[0,i])) 
            object[idx,i] <- FALSE
        else if (is.character(object[0,i]))
            object[idx,i] <- df[1,i]
        else if (is.factor(object[0,i]))
                object[idx,i] <- levels(object[0,i])[1]
            else object[idx,i] <- 0                
    }
    structure(object,na.action=idx)
}

estimate0 <- function(x,data=parent.frame(),pred.norm=c(),unstruct=FALSE,silent=TRUE,id=NULL,distribution=NULL,estimator="gaussian",...) {
    cl <- match.call()
    formulaId <- union(Specials(x,c("cluster")),Specials(x,c("id")))
    formulaSt <- paste0("~.-cluster(",formulaId,")-id(",formulaId,")")
    if (!is.null(formulaId)) {
        id <- formulaId
        x <- update(x,as.formula(formulaSt))
    }
    if (!is.null(id))
        x <- update(x,as.formula(paste(".~.+",id)))
    varnames <- all.vars(x)
    mf <- model.frame(x,data)
    mt <- attr(mf, "terms")
    yvar <- names(mf)[1]
    y <- mf[,yvar]
    opt <- options(na.action="na.pass")
    mm <- model.matrix(x,data)
    options(opt)
    covars <- colnames(mm)
    covars <- unlist(lapply(covars, function(x) gsub("[^a-zA-Z0-9._]","",x)))
    colnames(mm) <- covars

    if (attr(terms(x),"intercept")==1) {
        covars <- covars[-1]
        it <- c()
    } else {
        it <- "0"
    }

    if (!is.null(id)) covars <- setdiff(covars,id)
    model <- lvm(toformula(yvar,c(it,covars)),silent=TRUE)
    if (!is.null(distribution)) {
        lava::distribution(model,yvar) <- distribution
        estimator <- "glm"
    }
    mydata <- na.omit(as.data.frame(cbind(data.frame(y),mm))); names(mydata)[1] <- yvar
    exogenous(model) <- setdiff(covars,pred.norm)
    if (unstruct) {
        model <- covariance(model,pred.norm,pairwise=TRUE)
    }
    estimate(model,mydata,silent=silent,id=id,estimator=estimator,...)
}

