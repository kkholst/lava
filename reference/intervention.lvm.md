# Define intervention

Define intervention in a \`lvm\` object

## Usage

``` r
# S3 method for class 'lvm'
intervention(object, to, value, dist = none.lvm(), ...)
```

## Arguments

- object:

  lvm object

- to:

  String defining variable or formula

- value:

  function defining intervention

- dist:

  Distribution

- ...:

  Additional arguments to lower level functions

## See also

regression lvm sim

## Examples

``` r
m <- lvm(y ~ a + x, a ~ x)
distribution(m, ~a+y) <- binomial.lvm()
mm <- intervention(m, "a", value=3)
sim(mm, 10)
#>    y a          x
#> 1  1 3  0.4577942
#> 2  1 3  1.6463162
#> 3  1 3  1.3117039
#> 4  1 3  0.1758473
#> 5  1 3 -0.7003017
#> 6  1 3  0.8740681
#> 7  0 3 -0.7909800
#> 8  1 3  0.5525168
#> 9  1 3  0.3305292
#> 10 1 3 -0.8173789
mm <- intervention(m, a~x, function(x) (x>0)*1)
sim(mm, 10)
#>    y a          x
#> 1  0 0 -0.2232628
#> 2  0 0 -0.8529459
#> 3  1 1  1.0928103
#> 4  0 1  0.1812314
#> 5  0 0 -0.6480695
#> 6  1 1  0.3027412
#> 7  0 0 -0.9283982
#> 8  0 0 -0.4408714
#> 9  0 0 -0.3724818
#> 10 1 0 -0.8176981
```
