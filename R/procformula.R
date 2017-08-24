##' @export
procformula <- function(object=NULL,value,exo=lava.options()$exogenous,...) {

    ## Split into reponse and covariates by ~ disregarding expressions in parantheses
    ## '(?!...)' Negative lookahead assertion
    regex <- "~(?![^\\(].*\\))"
    yx <- lapply(strsplit(as.character(value),regex,perl=TRUE),function(x) gsub(" ","",x))[-1]
    yx <- lapply(yx,function(x) gsub("\n","",x))
    iscovar <- FALSE
    if (length(yx)==1) {
        lhs <- NULL; xidx <- 1
    } else {
            lhs <- yx[1]; xidx <- 2
            if (yx[[xidx]][1]=="") {
                yx[[xidx]] <- yx[[xidx]][-1]
                iscovar <- TRUE
            }
    }
    ##Check for link function
    invlink <- NULL
    if (xidx==2) {
        if (length(grep("[a-zA-Z0-9_]*\\(.*\\)$",yx[[xidx]]))>0) { ## rhs of the form F(x+y)
            invlink <- strsplit(yx[[xidx]],"\\(.*\\)")[[1]][1]
                if (invlink%in%c("f","v","I","") ||
                    grepl("\\+",invlink))
                { ## Reserved for setting linear constraints
                    invlink <- NULL
                } else {
                    yx[[xidx]] <- gsub(paste0(invlink,"\\(|\\)$"),"",yx[[xidx]])
                }
            }
    }

    ## Handling constraints with negative coefficients
    ## while not tampering with formulas like y~f(x,-2)
    st <- yx[[xidx]]
    st <- gsub("\\-","\\+\\-",gsub("\\+\\-","\\-",st)) ## Convert - to +- (to allow for splitting on '+')
    ##gsub("[^,]\\-","\\+\\-",st) ## Convert back any - not starting with ','
    st <- gsub(",\\+",",",st) ## Remove + inside 'f' and 'v' constraints
    st <- gsub("^\\+","",st) ## Remove leading plus
    yx[[xidx]] <- st
    
    ## Match '+' but not when preceeded by ( ... )
    X <- strsplit(yx[[xidx]],"\\+(?![^\\(]*\\))", perl=TRUE)[[1]]
    
    ##regex <- "(?!(\\(*))[\\(\\)]"
    regex <- "[\\(\\)]"
    ## Keep squares brackets and |(...) statements
    ## Extract variables from expressions like
    ## f(x,b) -> x,b  and  2*x -> 2,cx
    ## but avoid to tamper with transformation expressions:
    ## a~(x*b)
    res <- lapply(X,decomp.specials,regex,pattern2="\\*",pattern.ignore="~",reverse=TRUE,perl=TRUE)
    ##OLD:
    ##res <- lapply(X,decomp.specials,pattern2="[*]",reverse=TRUE)
    xx <- unlist(lapply(res, function(x) x[1]))

    xxf <- lapply(as.list(xx),function(x) decomp.specials(x,NULL,pattern2="\\[|~",perl=TRUE))
    xs <- unlist(lapply(xxf,function(x) x[1]))

    ## Alter intercepts?    
    intpos <- which(vapply(xs,function(x) grepl("^[\\+\\-]*[\\.|0-9]+$",x), 0)==1)
    ## Match '(int)'
    intpos0 <- which(vapply(X,function(x) grepl("^\\([\\+\\-]*[\\.|0-9]+\\)$",x),0)==1)

    yy <- ys <- NULL
    if (length(lhs)>0) {
        yy <- decomp.specials(lhs)
        yyf <- lapply(yy,function(y) decomp.specials(y,NULL,pattern2="[",fixed=TRUE,perl=FALSE))
        ys <- unlist(lapply(yyf,function(x) x[1]))
    }
    
    notexo <- c()
    if (!is.null(object)) {
      if (length(lhs)>0) {
        object <- addvar(object,ys,reindex=FALSE,...)
        notexo <- ys
        ## Add link-transformation
        if (!is.null(invlink)) {
            if (invlink=="") {
                object <- transform(object,ys,NULL,post=FALSE)
                covariance(object,ys) <- NA
            } else {
                ff <- function(x) {};  body(ff) <- parse(text=paste0(invlink,"(x)"))
                object <- transform(object,ys,ff,post=FALSE)
                covariance(object,ys) <- 0
            }
        }
      }

      if (length(intpos>0)) {
          xs[intpos[1]] <- gsub("\\+","",xs[intpos[1]])
          if (xs[intpos[1]]==1 && (!length(intpos0)>0) ) {
              xs[intpos[1]] <- NA
          }
          intercept(object,ys) <- char2num(xs[intpos[1]])
          xs <- xs[-intpos]
          res[intpos] <- NULL
      }

        object <- addvar(object,xs,reindex=FALSE ,...)
        exolist <- c()

        for (i in seq_len(length(xs))) {

            ## Extract transformation statements: var~(expr)
            xf0 <- strsplit(xx[[i]],"~")[[1]]
            if (length(xf0)>1) {
                myexpr <- xf0[2]
                ftr <- toformula(y="",x=paste0("-1+I(",myexpr,")"))
                xtr <- all.vars(ftr)
                xf0 <- xf0[1]
                transform(object, y=xf0, x=xtr) <- function(x) {
                    structure(model.matrix(ftr,as.data.frame(x)),dimnames=list(NULL,xf0))
                }
            }

            xf <- unlist(strsplit(xf0,"[\\[\\]]",perl=TRUE))
            if (length(xf)>1) {
                xpar <- strsplit(xf[2],":")[[1]]
                if (length(xpar)>1) {
                    val <- ifelse(xpar[2]=="NA",NA,xpar[2])
                    valn <- char2num(val)
                    covariance(object,xs[i]) <- ifelse(is.na(valn),val,valn)
                }
                val <- ifelse(xpar[1]=="NA",NA,xpar[1])
                valn <- char2num(val)
                if (is.na(val) || val!=".") {
                    intercept(object,xs[i]) <- ifelse(is.na(valn),val,valn)
                    notexo <- c(notexo,xs[i])
                }
            } else { exolist <- c(exolist,xs[i]) }
        }

        for (i in seq_len(length(ys))) {
            y <- ys[i]
            yf <- unlist(strsplit(yy[i],"[\\[\\]]",perl=TRUE))
            if (length(yf)>1) {
                ypar <- strsplit(yf[2],":")[[1]]
                if (length(ypar)>1) {
                    val <- ifelse(ypar[2]=="NA",NA,ypar[2])
                    valn <- char2num(val)
                    covariance(object,y) <- ifelse(is.na(valn),val,valn)
                }
                val <- ifelse(ypar[1]=="NA",NA,ypar[1])
                valn <- char2num(val)
                if (is.na(val) || val!=".") {
                    intercept(object,y) <- ifelse(is.na(valn),val,valn)
                }
            }
        }

        curvar <- index(object)$var
        if (exo) {
            oldexo <- exogenous(object)
            newexo <- setdiff(exolist,c(notexo,curvar,ys))
            exogenous(object) <- union(newexo,setdiff(oldexo,notexo))
        }
    }

    return(list(object=object,
                yx=yx,
                X=X,
                ys=ys,
                xx=xx,
                xs=xs,
                yy=yy,
                ys=ys,
                res=res,
                notexo=notexo,
                intpos=intpos,
                invlink=invlink,
                lhs=lhs,
                iscovar=iscovar))
}
