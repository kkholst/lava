\name{colorbar}
\alias{colorbar}
\title{Add color-bar to plot}
\usage{
colorbar(clut = Col(rev(rainbow(11, start = 0, end = 0.69)), 0.5),
  x.range = c(-0.5, 0.5), y.range = c(-0.1, 0.1), values = seq(clut),
  digits = 2, label.offset, srt = 45, cex = 0.5, border = NA,
  position = 1, direction = c("horizontal", "vertical"), ...)
}
\arguments{
  \item{clut}{Color look-up table}

  \item{0.5)}{}

  \item{x.range}{x range}

  \item{y.range}{y range}

  \item{values}{label values}

  \item{digits}{number of digits}

  \item{label.offset}{label offset}

  \item{srt}{rotation of labels}

  \item{cex}{text size}

  \item{border}{border of color bar rectangles}

  \item{position}{Label position left/bottom (1) or
  top/right (2) or no text (0)}

  \item{direction}{horizontal or vertical color bars}

  \item{...}{}

  \item{\dots}{additional low level arguments (i.e. parsed
  to \code{text})}
}
\description{
Add color-bar to plot
}
\examples{
\dontrun{
plotNeuro(x,roi=R,mm=-18,range=5)
colorbar(clut <- Col(rev(rainbow(11,start=0,end=0.69)),0.5),x=c(-40,40),y.range=c(84,90),values=c(-5:5))

colorbar(clut <- Col(rev(rainbow(11,start=0,end=0.69)),0.5),x=c(-10,10),y.range=c(-100,50),values=c(-5:5),direction="vertical",border=1)
}
}

