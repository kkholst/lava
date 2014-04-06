##' Convert PDF file to print quality png (default 300 dpi)
##'
##' Access to ghostscript program 'gs' is needed
##' @title Convert pdf to raster format
##' @param files Vector of (pdf-)filenames to process
##' @param dpi DPI
##' @param resolution Resolution of raster image file
##' @param gsopt Optional ghostscript arguments
##' @param format Raster format (e.g. png, jpg, tif, ...)
##' @param \dots Additional arguments
##' @seealso \code{dev.copy2pdf}, \code{printdev}
##' @export
##' @author Klaus K. Holst
##' @keywords iplot
pdfconvert <- function(files, dpi=300, resolution=1024, gsopt, format="png", ...) {
  if (missing(gsopt))
    gsopt <- "-dSAFTER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -dTextAlphaBits=4"
  cmd1 <- paste("gs -r",dpi," -dBackgroundColor='16#ffffff'",sep="")
  cmd2 <- paste("mogrify -resize ", resolution,sep="")
  for (f in files) {
    f0 <- strsplit(f,".pdf")[1]
    f.out <- paste(f0,format,sep=".")
    f.pdf <- paste(f0,"pdf",sep=".")
    mycmd1 <- paste(cmd1, " ", gsopt, " -sOutputFile=", f.out, " > /dev/null ", f.pdf, sep="")
    ##                   " && ", cmd2, " ", f.png, sep="")
    mycmd2 <- paste(cmd2, " ", f.out, sep="")
    cat(f.pdf)
    system(mycmd1)
    cat(" -> ")
    system(mycmd2)
    cat(f.out, "\n")
  }
}
