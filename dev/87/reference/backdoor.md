# Backdoor criterion

Check backdoor criterion of a lvm object

## Usage

``` r
backdoor(object, f, cond, ..., return.graph = FALSE)
```

## Arguments

- object:

  lvm object

- f:

  formula. Conditioning, z, set can be given as y~x\|z

- cond:

  Vector of variables to conditon on

- ...:

  Additional arguments to lower level functions

- return.graph:

  Return moral ancestral graph with z and effects from x removed

## Examples

``` r
m <- lvm(y~c2,c2~c1,x~c1,m1~x,y~m1, v1~c3, x~c3,v1~y,
         x~z1, z2~z1, z2~z3, y~z3+z2+g1+g2+g3)
ll <- backdoor(m, y~x)
backdoor(m, y~x|c1+z1+g1)
#> [1] TRUE
```
