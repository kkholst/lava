##' Backdoor criterion
##'
##' Check backdoor criterion of a lvm object
##' @param object lvm object
##' @param f formula. Conditioning, z, set can be given as y~x|z
##' @param cond Vector of variables to conditon on
##' @param ... Additional arguments to lower level functions
##' @examples
##' m <- lvm(y~c2,c2~c1,x~c1,m1~x,y~m1, v1~c3, x~c3,v1~y,
##'          x~z1, z2~z1, z2~z3, y~z3+z2+g1+g2+g3)
##' ll <- backdoor(m, y~x)
##' backdoor(m, y~x|c1+z1+g1)
##' @export
backdoor <- function(object, f, cond, ...) {
    y <- getoutcome(f, sep = "|")
    x <- attr(y, "x")
    if (length(x) > 1) {
        cond <- all.vars(x[[2]])
    }
    x <- all.vars(x[[1]])
    nod <- vars(object)
    des <- descendants(object, x)
    ch <- children(object, x)
    g0 <- cancel(object, toformula(x, ch))

    if (!base::missing(cond)) {
        return(dsep(g0, c(y, x), cond = cond) && !any(cond %in% des))
    }
    
    cset <- base::setdiff(nod, c(des, x, y)) ## possible conditioning set
    pp <- path(g0,y~x,all=TRUE) ## All backdoor paths
    M <- adjMat(g0)
    Collider <- function(vec) {
        M[vec[2],vec[1]] & M[vec[2],vec[3]]
    }
    blockList <- collideList <- c()
    for (i in seq_along(pp)) {
        p0 <- pp[[i]]
        blocks <- c()
        collide <- c()
        for (j in seq(length(p0)-2)) {
            if (Collider(p0[0:2 + j])) {
                collide <- c(collide,p0[1+j])
            } else {
                blocks <- c(blocks,p0[1+j])
            }
        }
        blockList <- c(blockList,list(blocks))
        collideList <- c(collideList,list(collide))
    }
    res <- list(blockList)
    ## Paths with colliders:
    col <- unlist(lapply(collideList,function(x) !is.null(x)))
    if (length(col)>0) col <- which(col)
    ## List of variables which are not on path between x and y:
    optional <- setdiff(cset,c(unlist(collideList),unlist(blockList)))
    callrecurs <- function(col,res=list()) {
        if (length(col)==0) return(res)
        blockList0 <- blockList
        blockList0[col] <- NULL
        blockList0 <- lapply(blockList0, function(x) setdiff(x,unlist(collideList[col])))
        if (!any(unlist(lapply(blockList0,is.null)))) {
            res <- c(res, list(blockList0))
        }            
        for (i in seq_along(col)) {
            col0 <- col[-i]
            if (length(col0)>0)
                res <- callrecurs(col0,res)
        }
        return(res)
    }
    if (length(col)>0)
        res <- c(res,callrecurs(col))
    ## Any element can be included from 'optional' For a given element
    ## in 'include' at least one element in each member of the list
    ## must be included in the conditioning set.
    return(list(optional=optional, include=res))
}
