# Extract pathways in model graph

Extract all possible paths from one variable to another connected
component in a latent variable model. In an estimated model the effect
size is decomposed into direct, indirect and total effects including
approximate standard errors.

## Usage

``` r
# S3 method for class 'lvm'
path (object, to = NULL, from, all=FALSE, ...)
# S3 method for class 'lvmfit'
effects (object, to, from, ...)
```

## Arguments

- object:

  Model object (`lvm`)

- ...:

  Additional arguments to be passed to the low level functions

- to:

  Outcome variable (string). Alternatively a formula specifying response
  and predictor in which case the argument `from` is ignored.

- from:

  Response variable (string), not necessarily directly affected by `to`.

- all:

  If TRUE all simple paths (in undirected graph) is returned on/off.

## Value

If `object` is of class `lvmfit` a list with the following elements is
returned

- idx:

  A list where each element defines a possible pathway via a integer
  vector indicating the index of the visited nodes.

- V :

  A List of covariance matrices for each path.

- coef :

  A list of parameters estimates for each path

- path :

  A list where each element defines a possible pathway via a character
  vector naming the visited nodes in order.

- edges :

  Description of 'comp2'

If `object` is of class `lvm` only the `path` element will be returned.

The `effects` method returns an object of class `effects`.

## Note

For a `lvmfit`-object the parameters estimates and their corresponding
covariance matrix are also returned. The `effects`-function additionally
calculates the total and indirect effects with approximate standard
errors

## See also

`children`, `parents`

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm(c(y1,y2,y3)~eta)
regression(m) <- y2~x1
latent(m) <- ~eta
regression(m) <- eta~x1+x2
d <- sim(m,500)
e <- estimate(m,d)

path(Model(e),y2~x1)
#> [[1]]
#> [1] "x1" "y2"
#> 
#> [[2]]
#> [1] "x1"  "eta" "y2" 
#> 
parents(Model(e), ~y2)
#> [1] "eta" "x1" 
children(Model(e), ~x2)
#> [1] "eta"
children(Model(e), ~x2+eta)
#> [1] "eta" "y1"  "y2"  "y3" 
effects(e,y2~x1)
#>           Estimate Std.Err z value   Pr(>|z|)
#> Total       2.0041 0.06400   31.31 3.111e-215
#> Direct      0.9491 0.07118   13.33  1.460e-40
#> Indirect    1.0549 0.07285   14.48  1.599e-47
#> y2~eta~x1   1.0549 0.07285   14.48  1.599e-47
#> 
#>                      Estimate   2.5%  97.5%
#> Mediation proportion   0.5264 0.4633 0.5895
## All simple paths (undirected)
path(m,y1~x1,all=TRUE)
#> [[1]]
#> [1] "x1"  "y2"  "eta" "y1" 
#> 
#> [[2]]
#> [1] "x1"  "eta" "y1" 
#> 
```
