# Newton-Raphson method

Newton-Raphson method

## Usage

``` r
NR(
  start,
  objective = NULL,
  gradient = NULL,
  hessian = NULL,
  control,
  args = NULL,
  ...
)
```

## Arguments

- start:

  Starting value

- objective:

  Optional objective function (used for selecting step length)

- gradient:

  gradient

- hessian:

  hessian (if NULL a numerical derivative is used)

- control:

  optimization arguments (see details)

- args:

  Optional list of arguments parsed to objective, gradient and hessian

- ...:

  additional arguments parsed to lower level functions

## Details

`control` should be a list with one or more of the following components:

- trace integer for which output is printed each 'trace'th iteration

- iter.max number of iterations

- stepsize: Step size (default 1)

- nstepsize: Increase stepsize every nstepsize iteration (from stepsize
  to 1)

- tol: Convergence criterion (gradient)

- epsilon: threshold used in pseudo-inverse

- backtrack: In each iteration reduce stepsize unless solution is
  improved according to criterion (gradient, armijo, curvature, wolfe)

## Examples

``` r
# Objective function with gradient and hessian as attributes
f <- function(z) {
   x <- z[1]; y <- z[2]
   val <- x^2 + x*y^2 + x + y
   structure(val, gradient=c(2*x+y^2+1, 2*y*x+1),
             hessian=rbind(c(2,2*y),c(2*y,2*x)))
}
NR(c(0,0),f)
#> $par
#> [1] -0.7324166  0.6825751
#> 
#> $iterations
#> [1] 12
#> 
#> $method
#> [1] "NR"
#> 
#> $gradient
#> [1] 2.451187e-07 7.301897e-07
#> 
#> $iH
#>            [,1]       [,2]
#> [1,] -0.3054540 -0.2849143
#> [2,] -0.2849143  0.4172596
#> attr(,"det")
#> [1] 1.937717
#> attr(,"pseudo")
#> [1] FALSE
#> attr(,"minSV")
#> [1] -2.473621
#> 

# Parsing arguments to the function and
g <- function(x,y) (x*y+1)^2
NR(0, gradient=g, args=list(y=2), control=list(trace=1,tol=1e-20))
#> 
#> Iter=0   ;   
#>      p= 0 
#> [1] "Numerical Hessian"
#> Iter=1   ;
#>  D= 1 
#>  p= -0.25 
#> [1] "Numerical Hessian"
#> Iter=2   ;
#>  D= 0.25 
#>  p= -0.375 
#> [1] "Numerical Hessian"
#> Iter=3   ;
#>  D= 0.0625 
#>  p= -0.4375 
#> [1] "Numerical Hessian"
#> Iter=4   ;
#>  D= 0.01562 
#>  p= -0.4688 
#> [1] "Numerical Hessian"
#> Iter=5   ;
#>  D= 0.003906 
#>  p= -0.4844 
#> [1] "Numerical Hessian"
#> Iter=6   ;
#>  D= 0.0009766 
#>  p= -0.4922 
#> [1] "Numerical Hessian"
#> Iter=7   ;
#>  D= 0.0002441 
#>  p= -0.4961 
#> [1] "Numerical Hessian"
#> Iter=8   ;
#>  D= 6.104e-05 
#>  p= -0.498 
#> [1] "Numerical Hessian"
#> Iter=9   ;
#>  D= 1.526e-05 
#>  p= -0.499 
#> [1] "Numerical Hessian"
#> Iter=10  ;
#>  D= 3.815e-06 
#>  p= -0.4995 
#> [1] "Numerical Hessian"
#> Iter=11  ;
#>  D= 9.537e-07 
#>  p= -0.4998 
#> [1] "Numerical Hessian"
#> Iter=12  ;
#>  D= 2.384e-07 
#>  p= -0.4999 
#> [1] "Numerical Hessian"
#> Iter=13  ;
#>  D= 5.96e-08 
#>  p= -0.4999 
#> [1] "Numerical Hessian"
#> Iter=14  ;
#>  D= 1.49e-08 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=15  ;
#>  D= 3.725e-09 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=16  ;
#>  D= 9.313e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=17  ;
#>  D= 2.328e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=18  ;
#>  D= 5.821e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=19  ;
#>  D= 1.455e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=20  ;
#>  D= 3.638e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=21  ;
#>  D= 9.095e-13 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=22  ;
#>  D= 2.274e-13 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=23  ;
#>  D= 5.684e-14 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=24  ;
#>  D= 1.421e-14 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=25  ;
#>  D= 3.553e-15 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=26  ;
#>  D= 8.882e-16 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=27  ;
#>  D= 2.22e-16 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=28  ;
#>  D= 5.551e-17 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=29  ;
#>  D= 1.388e-17 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=30  ;
#>  D= 3.469e-18 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=31  ;
#>  D= 8.674e-19 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=32  ;
#>  D= 2.168e-19 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=33  ;
#>  D= 5.421e-20 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=34  ;
#>  D= 1.355e-20 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=35  ;
#>  D= 3.388e-21 
#>  p= -0.5 
#> $par
#> [1] -0.5
#> 
#> $iterations
#> [1] 35
#> 
#> $method
#> [1] "NR"
#> 
#> $gradient
#> [1] 3.388132e-21
#> 
#> $iH
#>             [,1]
#> [1,] -4294965233
#> attr(,"det")
#>               [,1]
#> [1,] -2.328308e-10
#> 

```
