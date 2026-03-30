# Remove variables from (model) object.

Generic method for removing elements of object

## Usage

``` r
rmvar(x, ...) <- value
```

## Arguments

- x:

  Model object

- ...:

  additional arguments to lower level functions

- value:

  Vector of variables or formula specifying which nodes to remove

## See also

`cancel`

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm()
addvar(m) <- ~y1+y2+x
covariance(m) <- y1~y2
regression(m) <- c(y1,y2) ~ x
### Cancel the covariance between the residuals of y1 and y2
cancel(m) <- y1~y2
### Remove y2 from the model
rmvar(m) <- ~y2
```
