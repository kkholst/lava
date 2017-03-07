##' @export
parsedesign <- function(coef,x,...,regex=FALSE) {
    if (!is.vector(coef)) coef <- stats::coef(coef)
    if (is.numeric(coef) && !is.null(names(coef))) coef <- names(coef)
    dots <- substitute(list(...))[-1]
    expr <- suppressWarnings(inherits(try(x,silent=TRUE),"try-error"))
    if (expr) {
        ee <- c(deparse(substitute(x)), sapply(dots, deparse))
    } else {
        ee <- c(deparse(x), sapply(dots, function(x) deparse(x)))
    }
    res <- c()
    for (e in ee) {
        e0 <- gsub(" ","",e)
        ff <- strsplit(e0,'\"')[[1]]
        Val0 <- rbind(rep(0,length(coef)))
        Val <- c()
        for (i in seq(length(ff)/2)) {
            val0 <- gsub("[*()]","",ff[2*(i-1)+1])            
            suppressWarnings(val <- as.numeric(val0))
            if (is.na(val)) {
                val <- switch(val0,"-"=-1,1)
            }
            par0 <- ff[2*i]
            par0int <- suppressWarnings(as.integer(par0))
            if (!regex) par0 <- glob2rx(par0)
            if (is.na(par0int)) par0int <- grep(par0,coef)
            if (length(par0int)>1) {
                for (k in seq_along(par0int)) {
                    if (par0int[k]<=length(Val0)) {
                        Val00 <- Val0
                        Val00[par0int[k]] <- val
                        Val <- rbind(Val,Val00)
                    }
                }
            } else {
                if (length(par0int)>0 && par0int<=length(Val0)) Val0[par0int] <- val
            }
        }
        if (length(Val)==0) Val <- Val0
        if (any(Val!=0)) res <- rbind(res,Val)
    }
    res
}
