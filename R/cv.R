rmse1 <- function(fit,data,response=NULL,...) {
    yhat <- predict(fit,newdata=data,...)
    if (is.null(response)) response <- endogenous(fit)
    y <- data[,response]
    c(RMSE=mean(as.matrix(y-yhat)^2))
}

##' Cross-validation
##'
##' Generic cross-validation function
##' @title Cross-validation
##' @param modelList List of fitting functions or models
##' @param data data.frame
##' @param K Number of folds (default 5, 0 splits in 1:n/2, n/2:n with last part used for testing)
##' @param rep Number of repetitions (default 1)
##' @param perf Performance measure (default RMSE)
##' @param seed Optional random seed
##' @param mc.cores Number of cores used for parallel computations
##' @param shared function applied to each fold with results send to each model
##' @param ... Additional arguments parsed to models in modelList and perf
##' @author Klaus K. Holst
##' @examples
##' f0 <- function(data,...) lm(...,data)
##' f1 <- function(data,...) lm(Sepal.Length~Species,data)
##' f2 <- function(data,...) lm(Sepal.Length~Species+Petal.Length,data)
##' x <- cv(list(m0=f0,m1=f1,m2=f2),rep=10, data=iris, formula=Sepal.Length~.)
##' x2 <- cv(list(f0(iris),f1(iris),f2(iris)),rep=10, data=iris)
##' @export
cv <- function(modelList, data, K=5, rep=1, perf, seed=NULL, mc.cores=1, shared=NULL, ...) {
    if (is.vector(data)) data <- cbind(data)
    if (missing(perf)) perf <- rmse1
    if (!is.list(modelList)) modelList <- list(modelList)
    nam <- names(modelList)
    if (is.null(nam)) nam <- paste0("model",seq_along(modelList))
    args0 <- list(...)
    args <- args0
    if (!is.null(shared)) {
        sharedres <- shared(data,...)
        args <- c(args, sharedres)
    }
    ## Models run on full data:
    if (is.function(modelList[[1]])) {
        fit0 <- lapply(modelList, function(f) do.call(f,c(list(data),args)))
    } else {
        fit0 <- modelList
    }
    ## In-sample predictive performance:
    perf0 <- lapply(fit0, function(fit) do.call(perf,c(list(fit,data=data),args)))
    namPerf <- names(perf0[[1]])
    names(fit0) <- names(perf0) <- nam
    n <- NROW(data)
    M <- length(perf0)      # Number of models
    P <- length(perf0[[1]]) # Number of performance measures
    if (!is.null(seed)) {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1)
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    if (K==0) {
        rep <- 1
        K <- 1
        folds <- list(csplit(seq(n)))
    } else {
        folds <- foldr(n,K,rep)
    }
    arg <- expand.grid(R=seq(rep),K=seq(K)) #,M=seq_along(modelList))
    dim <- c(rep,K,M,P)
    PerfArr <- array(0,dim)
    dimnames(PerfArr) <- list(NULL,NULL,nam,namPerf)

    ff <- function(i) {
        R <- arg[i,1]
        k <- arg[i,2]
        fold <- folds[[R]]
        dtest <- data[fold[[k]],,drop=FALSE]
        dtrain <- data[unlist(fold[-k]),,drop=FALSE]
        args <- args0
        if (!is.null(shared)) {
            sharedres <- shared(dtrain,...)
            args <- c(args, sharedres)
        }
        if (is.function(modelList[[1]])) {
            fits <- lapply(modelList, function(f) do.call(f,c(list(dtrain),args)))
        } else {
            fits <- lapply(modelList, function(m) do.call(update,c(list(m,data=dtrain),args)))
        }
        perfs <- lapply(fits, function(fit) do.call(perf,c(list(fit,data=dtest),args)))
        do.call(rbind,perfs)
    }

    if (mc.cores>1) {
        val <- parallel::mcmapply(ff,seq(nrow(arg)),SIMPLIFY=FALSE,mc.cores=mc.cores)
    } else {
        val <- mapply(ff,seq(nrow(arg)),SIMPLIFY=FALSE)
    }
    for (i in seq(nrow(arg))) {
        R <- arg[i,1]
        k <- arg[i,2]
        PerfArr[R,k,,] <- val[[i]]
    }

    structure(list(cv=PerfArr,
                   call=match.call(),
                   names=nam,
                   rep=rep, folds=K,
                   fit=fit0),
              class="CrossValidated")
}

##' @export
print.CrossValidated <- function(x,...) {
    res <- coef(x)
    print(res,quote=FALSE)
}

coef.CrossValidated <- function(object,...) {
    res <- apply(object$cv,3:4,function(x) mean(x))
    if (length(object$names)==nrow(res)) rownames(res) <- object$names
    res
}

