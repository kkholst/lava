# Check d-separation criterion

Check for conditional independence (d-separation)

## Usage

``` r
# S3 method for class 'lvm'
dsep(object, x, cond = NULL, return.graph = FALSE, ...)
```

## Arguments

- object:

  lvm object

- x:

  Variables for which to check for conditional independence

- cond:

  Conditioning set

- return.graph:

  If TRUE the moralized ancestral graph with the conditioning set
  removed is returned

- ...:

  Additional arguments to lower level functions

## Details

The argument 'x' can be given as a formula, e.g. x~y\|z+v or ~x+y\|z+v
With everything on the rhs of the bar defining the variables on which to
condition on.

## Examples

``` r
m <- lvm(x5 ~ x4+x3, x4~x3+x1, x3~x2, x2~x1)
if (interactive()) {
plot(m,layoutType='neato')
}
dsep(m,x5~x1|x2+x4)
#> [1] FALSE
dsep(m,x5~x1|x3+x4)
#> [1] TRUE
dsep(m,~x1+x2+x3|x4)
#> [1] FALSE
```
