##' Missing value generator
##'
##' This function adds a binary variable to a given \code{lvm} model and also
##' a variable which is equal to the original variable where the binary variable is
##' 
##' @title Missing value generator
##' @param object  \code{lvm}-object.
##' @param formula The right hand side specifies the name of a latent variable
##' which is not always observed. The left hand side specifies the name of
##' a new variable which is equal to the latent variable but has missing values.
##' @param Rformula The left hand side specifies the name of the missing value
##' indicator
##' @param ... Passed to binomial.lvm.
##' @return lvm object
##' @examples
##' library(lava)
##' set.seed(17)
##' m <- lvm(y0~x01+x02+x03)
##' m <- Missing(m,formula=x1~x01,Rformula=R1~0.3*x02+-0.7*x01,p=0.4)
##' sim(m,10)
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
Missing <- function(object,formula,Rformula,...){
    name <- all.vars(Rformula)[1]
    distribution(object,name) <- binomial.lvm(...)
    transform(object,update(formula,paste(".~.+",name))) <- function(u){
        out <- u[,1]
        out[u[,2]==1] <- NA
        out
    }
    regression(object) <- Rformula
    object
}
