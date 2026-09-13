# Add color-bar to plot

Add color-bar to plot

## Usage

``` r
colorbar(
  clut = Col(rev(rainbow(11, start = 0, end = 0.69)), alpha),
  x.range = c(-0.5, 0.5),
  y.range = c(-0.1, 0.1),
  values = seq(clut),
  digits = 2,
  label.offset,
  srt = 45,
  cex = 0.5,
  border = NA,
  alpha = 0.5,
  position = 1,
  direction = c("horizontal", "vertical"),
  ...
)
```

## Arguments

- clut:

  Color look-up table

- x.range:

  x range

- y.range:

  y range

- values:

  label values

- digits:

  number of digits

- label.offset:

  label offset

- srt:

  rotation of labels

- cex:

  text size

- border:

  border of color bar rectangles

- alpha:

  Alpha (transparency) level 0-1

- position:

  Label position left/bottom (1) or top/right (2) or no text (0)

- direction:

  horizontal or vertical color bars

- ...:

  additional low level arguments (i.e. parsed to `text`)

## Examples

``` r
if (FALSE) { # \dontrun{
plotNeuro(x,roi=R,mm=-18,range=5)
colorbar(clut=Col(rev(rainbow(11,start=0,end=0.69)),0.5),
         x=c(-40,40),y.range=c(84,90),values=c(-5:5))

colorbar(clut=Col(rev(rainbow(11,start=0,end=0.69)),0.5),
         x=c(-10,10),y.range=c(-100,50),values=c(-5:5),
         direction="vertical",border=1)
} # }
```
