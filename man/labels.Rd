\name{labels}
\alias{labels}
\alias{labels<-}
\alias{labels<-.default}
\alias{labels.lvm}
\alias{labels.lvmfit}
\alias{labels.graphNEL}
\alias{edgelabels}
\alias{edgelabels<-}
\alias{edgelabels<-.lvm}
\alias{nodecolor}
\alias{nodecolor<-}
\alias{nodecolor<-.default}



\title{ Define labels of graph }
\description{
  Alters labels of nodes and edges in the graph of a latent variable model
}
\usage{
\method{labels}{default}(object, ...) <- value
\method{edgelabels}{lvm}(object, to, ...) <- value
\method{nodecolor}{default}(object, var=vars(object),
border, labcol, shape, lwd, ...) <- value
}
\arguments{
  \item{object}{\code{lvm}-object.}
  \item{value}{node label/edge label/color}
  \item{to}{Formula specifying outcomes and predictors defining relevant
    edges.}
  \item{\dots}{ Additional arguments (\code{lwd}, \code{cex},
    \code{col}, \code{labcol}), \code{border}. }
  \item{var}{Formula or character vector specifying the nodes/variables to alter.}
  \item{border}{Colors of borders}
  \item{labcol}{Text label colors}
  \item{shape}{Shape of node}
  \item{lwd}{Line width of border}
}
%%\details{
%%  Notice that for a \code{lvmfit}-object, adding labels with edgelabels,
%%  will result in added 
%%}
\author{ Klaus K. Holst }
\examples{

m <- lvm(y~x+z)
labels(m) <- c(y=expression(psi), z=expression(zeta))
nodecolor(m,~y+z+x,border=c("white","white","black"), labcol="white", lwd=c(1,1,5)) <- c("orange","indianred","lightgreen")
edgelabels(m,y~z+x, cex=c(2,3), col=c("orange","black"),labcol="darkblue",
lwd=c(3,1)) <- expression(phi,rho)
plot(m,addcolor=FALSE)

}
\keyword{ graphs }
\keyword{ aplot }
