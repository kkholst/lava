sumsplit <- function(x,...) {
    plus <- strsplit(x,"\\+")[[1]]
    spl <- unlist(lapply(plus, function(x) {
        val <- strsplit(x,"\\-")[[1]]
        val[-1] <- paste0("-",val[-1])
        setdiff(val,"")
    }))
    res <- c()
    for (st in spl) {
        st <- gsub(" ","",st)
        st0 <- gsub("^[-0-9\\*]*","",st)
        val <- gsub("\\*","",regmatches(st,gregexpr("^[-0-9\\*]*",st))[[1]])
        if (val=="") val <- "1"
        val <- switch(val,"-"=-1,val)
        res <- c(res,val,st0)    
    }
    return(res)
}

##' @export
parsedesign <- function(coef,x,...,regex=FALSE,diff=TRUE) {
    if (!is.vector(coef)) coef <- stats::coef(coef)
    if (is.numeric(coef) && !is.null(names(coef))) coef <- names(coef)    
    dots <- lapply(substitute(list(...)),function(x) x)[-1]
    expr <- suppressWarnings(inherits(try(x,silent=TRUE),"try-error"))
    if (expr) {
        ee <- c(deparse(substitute(x)), unlist(lapply(dots, deparse)))
    } else {
        ee <- c(deparse(x), sapply(dots, function(x) deparse(x)))
    }
    if (!expr && is.numeric(x)) {
        return(do.call(contr, list(c(list(x),list(...)),n=length(coef),diff=diff)))
    }
    res <- c()
    diff <- rep(diff,length.out=length(ee))
    count <- 0
    for (e in ee) {
        count <- count+1
        diff0 <- FALSE
        Val <- rbind(rep(0,length(coef)))
        if (grepl('\"',e)) {        
            diff0 <- diff[count] && grepl("^c\\(",e)
            e0 <- gsub(" |\\)$|^c\\(","",e)
            ff <- strsplit(e0,'\"')[[1L]]
        } else {
            ff <- sumsplit(e)
        }
        for (i in seq(length(ff)/2)) {
            val0 <- gsub("[*()]","",ff[2*(i-1)+1])            
            val <- char2num(val0)
            if (is.na(val)) {
                val <- switch(val0,"-"=-1,1)
            }
            par0 <- ff[2*i]
            par0int <- as.integer(char2num(par0))
            if (!regex) par0 <- glob2rx(par0)
            if (is.na(par0int)) par0int <- grep(par0,coef)
            if (length(par0int)>1) {
                diff0 <- diff[count]
                for (k in seq_along(par0int)) {
                    if (par0int[k]<=length(Val)) {
                        if (diff[count]) {
                            Val[par0int[k]] <- val
                        } else {
                            Val0 <- Val; Val0[] <- 0
                            Val0[par0int[k]] <- val
                            res <- rbind(res,Val0)
                        }
                    }
                }
            } else {
                if (length(par0int)>0 && par0int<=length(Val)) Val[par0int] <- val
            }
        }
        if (diff0) {
            n <- sum(Val!=0)
            if (n>1) {
                Val0 <- Val
                ii <- which(Val0!=0)
                Val <- matrix(0,nrow=n-1,ncol=length(Val))
                for (i in seq(n-1)) {
                    Val[i,ii[c(1,i+1)]] <- Val0[ii[c(1,i+1)]]*c(1,-1)
                }
            }
        }
        if (any(Val!=0)) res <- rbind(res,Val)
    }
    res
}
