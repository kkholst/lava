# Identify points on plot

Extension of the `identify` function

## Usage

``` r
# Default S3 method
click(x, y=NULL, label=TRUE, n=length(x), pch=19, col="orange", cex=3, ...)
idplot(x, y ,..., id=list(), return.data=FALSE)
```

## Arguments

- x:

  X coordinates

- ...:

  Additional arguments parsed to `plot` function

- y:

  Y coordinates

- label:

  Should labels be added?

- n:

  Max number of inputs to expect

- pch:

  Symbol

- col:

  Colour

- cex:

  Size

- id:

  List of arguments parsed to `click` function

- return.data:

  Boolean indicating if selected points should be returned

## Details

For the usual 'X11' device the identification process is terminated by
pressing any mouse button other than the first. For the 'quartz' device
the process is terminated by pressing either the pop-up menu equivalent
(usually second mouse button or 'Ctrl'-click) or the 'ESC' key.

## See also

`idplot`, `identify`

## Author

Klaus K. Holst

## Examples

``` r
if (interactive()) {
    n <- 10; x <- seq(n); y <- runif(n)
    plot(y ~ x); click(x,y)

    data(iris)
    l <- lm(Sepal.Length ~ Sepal.Width*Species,iris)
    res <- plotConf(l,var2="Species")## ylim=c(6,8), xlim=c(2.5,3.3))
    with(res, click(x,y))

    with(iris, idplot(Sepal.Length,Petal.Length))
}
```
