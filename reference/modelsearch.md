# Model searching

Performs Wald or score tests

## Usage

``` r
modelsearch(x, k = 1, dir = "forward", type = "all", ...)
```

## Arguments

- x:

  `lvmfit`-object

- k:

  Number of parameters to test simultaneously. For `equivalence` the
  number of additional associations to be added instead of `rel`.

- dir:

  Direction to do model search. "forward" := add associations/arrows to
  model/graph (score tests), "backward" := remove associations/arrows
  from model/graph (wald test)

- type:

  If equal to 'correlation' only consider score tests for covariance
  parameters. If equal to 'regression' go through direct effects only
  (default 'all' is to do both)

- ...:

  Additional arguments to be passed to the low level functions

## Value

Matrix of test-statistics and p-values

## See also

[`compare`](http://kkholst.github.io/lava/reference/compare.md),
[`equivalence`](http://kkholst.github.io/lava/reference/equivalence.md)

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm();
regression(m) <- c(y1,y2,y3) ~ eta; latent(m) <- ~eta
regression(m) <- eta ~ x
m0 <- m; regression(m0) <- y2 ~ x
dd <- sim(m0,100)[,manifest(m0)]
e <- estimate(m,dd);
modelsearch(e,messages=0)
#>  Score: S P(S>s) Index  holm BH    
#>  0.0734   0.7864 y3~~x  1    0.7864
#>  0.0734   0.7864 y3~x   1    0.7864
#>  0.0734   0.7864 x~y3   1    0.7864
#>  0.0734   0.7864 y1~~y2 1    0.7864
#>  0.0734   0.7864 y1~y2  1    0.7864
#>  0.0734   0.7864 y2~y1  1    0.7864
#>  0.1991   0.6554 y1~~x  1    0.7864
#>  0.1991   0.6554 y1~x   1    0.7864
#>  0.1991   0.6554 x~y1   1    0.7864
#>  0.1991   0.6554 y2~~y3 1    0.7864
#>  0.1991   0.6554 y2~y3  1    0.7864
#>  0.1991   0.6554 y3~y2  1    0.7864
#>  0.675    0.4113 y1~~y3 1    0.7864
#>  0.675    0.4113 y1~y3  1    0.7864
#>  0.675    0.4113 y3~y1  1    0.7864
#>  0.675    0.4113 y2~~x  1    0.7864
#>  0.675    0.4113 y2~x   1    0.7864
#>  0.675    0.4113 x~y2   1    0.7864
modelsearch(e,messages=0,type="cor")
#>  Score: S P(S>s) Index  holm BH    
#>  0.0734   0.7864 y3~~x  1    0.7864
#>  0.0734   0.7864 y1~~y2 1    0.7864
#>  0.1991   0.6554 y1~~x  1    0.7864
#>  0.1991   0.6554 y2~~y3 1    0.7864
#>  0.675    0.4113 y1~~y3 1    0.7864
#>  0.675    0.4113 y2~~x  1    0.7864
```
