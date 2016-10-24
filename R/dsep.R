##' @export
`dsep` <-
  function(object,...) UseMethod("dsep")

##' Check d-separation criterion
##'
##' Check for conditional independence (d-separation)
##' @export
##' @aliases dsep dsep.lvm
##' @param object lvm object
##' @param x Variables for which to check for conditional independence
##' @param cond Conditioning set
##' @param return.graph If TRUE the moralized ancestral graph with the
##'     conditioning set removed is returned
##' @param ... Additional arguments to lower level functions
##' @details The argument 'x' can be given as a formula, e.g.  x~y|z+v
##'     or ~x+y|z+v With everything on the rhs of the bar defining the
##'     variables on which to condition on.
##' @examples
##' m <- lvm(x5 ~ x4+x3, x4~x3+x1, x3~x2, x2~x1)
##' if (interactive()) {
##' plot(m,layoutType='neato')
##' }
##' dsep(m,x5~x1|x2+x4)
##' dsep(m,x5~x1|x3+x4)
##' dsep(m,~x1+x2+x3|x4)
##' 
dsep.lvm <- function(object,x,cond=NULL,return.graph=FALSE,...) {
    if (inherits(x,"formula")) {
        xf <- getoutcome(x,sep="|")
        xx <- attr(xf,"x")
        if (length(xx)==0) stop("Not a valid formula")
        x <- c(xf,all.vars(xx[[1]]))
        if (length(xx)>1) {
            cond <- all.vars(xx[[2]])
        }
    }
    if (inherits(cond,"formula")) {
        cond <- all.vars(cond)
    }
    nod <- vars(object)
    x <- intersect(x,nod)
    cond <- intersect(cond,nod)
    V <- c(x,cond)
    ## Ancenstral graph
    keep <- c(V,ancestors(object,V))
    del <- setdiff(nod,keep)
    if (length(del)>0) object <- rmvar(object,del)
    ## moralized graph
    man <- object
    for (v in V) {
        pa <- parents(object,v)
        if (length(pa)>1)  
            man$M[pa,pa] <- 1
        ## for (i in seq(length(pa)-1)) {
        ##     for (j in seq(length(pa)-1)+1) {
        ##         man$M[i,j]
        ##         man <- regression(man,from=pa[i],to=pa[j])
        ##     }
        ## }
    }    
    man.sel <- rmvar(man,cond)
    ## with(man.sel, solve(diag(nrow=nrow(M))-M))
    ii <- match(x,vars(man.sel))
    A <- with(man.sel, (t(M)+M)>0)
    dsep <- c()
    for (i in ii) {
        conn <- DFS(A,i)
        i0 <- setdiff(ii,i)
        dsep <- c(dsep,!any(i0%in%conn))
    }
    res <- all(dsep)
    attr(man.sel,"dsep") <- res
    if (return.graph) return(man.sel)
    return(res)
}
