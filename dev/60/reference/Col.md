# Generate a transparent RGB color

This function transforms a standard color (e.g. "red") into an
transparent RGB-color (i.e. alpha-blend\<1).

## Usage

``` r
Col(col, alpha = 0.2, locate = 0)
```

## Arguments

- col:

  color (numeric or character)

- alpha:

  degree of transparency (0-1)

- locate:

  Choose colour (with mouse)

## Value

A character vector with elements of 7 or 9 characters, `#` followed by
the red, blue, green and optionally alpha values in hexadecimal (after
rescaling to '0 ... 255').

## Details

This only works for certain graphics devices (Cairo-X11 (x11 as of
R\>=2.7), quartz, pdf, ...).

## Author

Klaus K. Holst

## Examples

``` r
plot(runif(1000),cex=runif(1000,0,4),
     col=Col(c("darkblue","orange"),0.5),pch=16)
```
