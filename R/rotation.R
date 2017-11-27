##' Performs a rotation in the plane
##'
##' @title Performs a rotation in the plane
##' @aliases rotate2 rot2D rot3D
##' @param x Matrix to be rotated (2 times n)
##' @param theta Rotation in radians
##' @return Returns a matrix of the same dimension as \code{x}
##' @author Klaus K. Holst
##' @export
##' @examples
##' rotate2(cbind(c(1,2),c(2,1)))
##' @keywords hplot
`rotate2` <-
function(x,theta=pi) {
  R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), byrow=TRUE, ncol=2)
  x%*%R
}

## clockwise rotation 2d:
##' @export
rot2D <- function(theta) {
  matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2)
}

##' @export
rot3D <- function(x=0,y=0,z=0) {
  Rx <- function() {
    R2 <- rot2D(x)
    R <- diag(3)
    R[2:3,2:3] <- R2
    return(R)
  }
  Ry <- function() {
    R2 <- rot2D(y)
    R <- diag(3)
    R[c(1,3),c(1,3)] <- R2
    return(R)
  }
  Rz <- function() {
    R2 <- rot2D(z)
    R <- diag(3)
    R[1:2,1:2] <- R2
    return(R)
  } 
  res <- diag(3)
  if (x!=0) res <- res%*%Rx()
  if (y!=0) res <- res%*%Ry()
  if (z!=0) res <- res%*%Rz()
  return(res)
}
