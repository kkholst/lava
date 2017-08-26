rmse1 <- function(fit,data,response=NULL,...) {
    yhat <- predict(fit,newdata=data)
    if (is.null(response)) response <- endogenous(fit)
    y <- data[,response]
    c(RMSE=mean(as.matrix(y-yhat)^2))
}

##' Cross-validation
##'
##' Generic cross-validation function
##' @title Cross-validation
##' @param modelList List of fitting functions
##' @param data data.frame
##' @param K Number of folds (default 5)
##' @param rep Number of repetitions (default 1)
##' @param perf Performance measure (default RMSE)
##' @param seed Optional random seed
##' @param ... Additional arguments parsed to models in modelList and perf
##' @export 
##' @author Klaus K. Holst
##' @examples
##' f0 <- function(data,...) lm(...,data)
##' f1 <- function(data,...) lm(Sepal.Length~Species,data)
##' f2 <- function(data,...) lm(Sepal.Length~Species+Petal.Length,data)
##' x <- cv(list(model0=f0,model1=f1,model2=f2),rep=10, data=iris, formula=Sepal.Length~.)
##' @export 
cv <- function(modelList, data, K=5, rep=1, perf, seed=NULL, mc.cores=1, ...) {
    if (missing(perf)) perf <- rmse1
    if (!is.list(modelList)) modelList <- list(modelList)
    nam <- names(modelList)
    args <- list(...)
    ## Models run on full data:
    fit0 <- lapply(modelList, function(f) do.call(f,c(list(data),args)))
    ## In-sample predictive performance:
    perf0 <- lapply(fit0, function(fit) do.call(perf,c(list(fit,data=data),args)))
    namPerf <- names(perf0[[1]])
    names(fit0) <- names(perf0) <- nam
    n <- nrow(data)
    M <- length(perf0)      # Number of models
    P <- length(perf0[[1]]) # Number of performance measures
    if (!is.null(seed)) {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))        
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }


    
    nam <- list(NULL,NULL,nam,namPerf)
    dim <- c(rep,K,M,P)
    PerfArr <- array(0,dim)
    dimnames(PerfArr) <- nam

    folds <- foldsr(n,K,rep)
    arg <- expand.grid(R=seq(rep),K=seq(K)) #,M=seq_along(modelList))

    ff <- function(i) {
        R <- arg[i,1]
        k <- arg[i,2]
        fold <- folds[[R]]
        dtest <- data[fold[[k]],]
        dtrain <- data[unlist(fold[-k]),]
        fits <- lapply(modelList, function(f) do.call(f,c(list(dtrain),args)))
        perfs <- lapply(fits, function(fit) do.call(perf,c(list(fit,data=dtest),args)))
        Reduce(rbind,perfs)        
    }
    if (mc.cores>1) {
        val <- parallel::mcmapply(ff,seq(nrow(arg)),SIMPLIFY=FALSE)
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
                   rep=rep, folds=K,
                   fit=fit0),
              class="CrossValidated")
}

##' @export
print.CrossValidated <- function(x,...) {
    ##print(drop(x$cv))
    print(apply(x$cv,3:4,function(x) mean(x)),quote=FALSE)
}

