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
#>  D= 0.06254 
#>  p= -0.4375 
#> [1] "Numerical Hessian"
#> Iter=4   ;
#>  D= 0.01565 
#>  p= -0.4687 
#> [1] "Numerical Hessian"
#> Iter=5   ;
#>  D= 0.003918 
#>  p= -0.4843 
#> [1] "Numerical Hessian"
#> Iter=6   ;
#>  D= 0.0009826 
#>  p= -0.4921 
#> [1] "Numerical Hessian"
#> Iter=7   ;
#>  D= 0.0002472 
#>  p= -0.496 
#> [1] "Numerical Hessian"
#> Iter=8   ;
#>  D= 6.259e-05 
#>  p= -0.498 
#> [1] "Numerical Hessian"
#> Iter=9   ;
#>  D= 1.604e-05 
#>  p= -0.499 
#> [1] "Numerical Hessian"
#> Iter=10  ;
#>  D= 4.208e-06 
#>  p= -0.4995 
#> [1] "Numerical Hessian"
#> Iter=11  ;
#>  D= 1.152e-06 
#>  p= -0.4997 
#> [1] "Numerical Hessian"
#> Iter=12  ;
#>  D= 3.392e-07 
#>  p= -0.4998 
#> [1] "Numerical Hessian"
#> Iter=13  ;
#>  D= 1.115e-07 
#>  p= -0.4999 
#> [1] "Numerical Hessian"
#> Iter=14  ;
#>  D= 4.219e-08 
#>  p= -0.4999 
#> [1] "Numerical Hessian"
#> Iter=15  ;
#>  D= 1.859e-08 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=16  ;
#>  D= 9.411e-09 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=17  ;
#>  D= 5.347e-09 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=18  ;
#>  D= 3.327e-09 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=19  ;
#>  D= 2.221e-09 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=20  ;
#>  D= 1.567e-09 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=21  ;
#>  D= 1.154e-09 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=22  ;
#>  D= 8.799e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=23  ;
#>  D= 6.901e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=24  ;
#>  D= 5.54e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=25  ;
#>  D= 4.535e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=26  ;
#>  D= 3.774e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=27  ;
#>  D= 3.185e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=28  ;
#>  D= 2.721e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=29  ;
#>  D= 2.349e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=30  ;
#>  D= 2.047e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=31  ;
#>  D= 1.799e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=32  ;
#>  D= 1.593e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=33  ;
#>  D= 1.419e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=34  ;
#>  D= 1.272e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=35  ;
#>  D= 1.146e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=36  ;
#>  D= 1.038e-10 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=37  ;
#>  D= 9.444e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=38  ;
#>  D= 8.626e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=39  ;
#>  D= 7.909e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=40  ;
#>  D= 7.276e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=41  ;
#>  D= 6.716e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=42  ;
#>  D= 6.217e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=43  ;
#>  D= 5.77e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=44  ;
#>  D= 5.37e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=45  ;
#>  D= 5.01e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=46  ;
#>  D= 4.684e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=47  ;
#>  D= 4.389e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=48  ;
#>  D= 4.12e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=49  ;
#>  D= 3.876e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=50  ;
#>  D= 3.652e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=51  ;
#>  D= 3.447e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=52  ;
#>  D= 3.258e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=53  ;
#>  D= 3.085e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=54  ;
#>  D= 2.925e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=55  ;
#>  D= 2.776e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=56  ;
#>  D= 2.639e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=57  ;
#>  D= 2.512e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=58  ;
#>  D= 2.393e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=59  ;
#>  D= 2.283e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=60  ;
#>  D= 2.18e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=61  ;
#>  D= 2.084e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=62  ;
#>  D= 1.994e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=63  ;
#>  D= 1.91e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=64  ;
#>  D= 1.83e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=65  ;
#>  D= 1.756e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=66  ;
#>  D= 1.686e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=67  ;
#>  D= 1.62e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=68  ;
#>  D= 1.558e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=69  ;
#>  D= 1.5e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=70  ;
#>  D= 1.444e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=71  ;
#>  D= 1.392e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=72  ;
#>  D= 1.342e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=73  ;
#>  D= 1.295e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=74  ;
#>  D= 1.251e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=75  ;
#>  D= 1.208e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=76  ;
#>  D= 1.168e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=77  ;
#>  D= 1.13e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=78  ;
#>  D= 1.093e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=79  ;
#>  D= 1.059e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=80  ;
#>  D= 1.026e-11 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=81  ;
#>  D= 9.939e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=82  ;
#>  D= 9.638e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=83  ;
#>  D= 9.35e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=84  ;
#>  D= 9.075e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=85  ;
#>  D= 8.811e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=86  ;
#>  D= 8.559e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=87  ;
#>  D= 8.317e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=88  ;
#>  D= 8.086e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=89  ;
#>  D= 7.864e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=90  ;
#>  D= 7.651e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=91  ;
#>  D= 7.446e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=92  ;
#>  D= 7.25e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=93  ;
#>  D= 7.061e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=94  ;
#>  D= 6.879e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=95  ;
#>  D= 6.705e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=96  ;
#>  D= 6.537e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=97  ;
#>  D= 6.375e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=98  ;
#>  D= 6.219e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=99  ;
#>  D= 6.068e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=100 ;
#>  D= 5.923e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=101 ;
#>  D= 5.783e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=102 ;
#>  D= 5.648e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=103 ;
#>  D= 5.518e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=104 ;
#>  D= 5.392e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=105 ;
#>  D= 5.27e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=106 ;
#>  D= 5.153e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=107 ;
#>  D= 5.039e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=108 ;
#>  D= 4.929e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=109 ;
#>  D= 4.822e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=110 ;
#>  D= 4.719e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=111 ;
#>  D= 4.62e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=112 ;
#>  D= 4.523e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=113 ;
#>  D= 4.429e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=114 ;
#>  D= 4.338e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=115 ;
#>  D= 4.25e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=116 ;
#>  D= 4.165e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=117 ;
#>  D= 4.082e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=118 ;
#>  D= 4.002e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=119 ;
#>  D= 3.923e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=120 ;
#>  D= 3.848e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=121 ;
#>  D= 3.774e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=122 ;
#>  D= 3.702e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=123 ;
#>  D= 3.633e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=124 ;
#>  D= 3.565e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=125 ;
#>  D= 3.499e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=126 ;
#>  D= 3.435e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=127 ;
#>  D= 3.373e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=128 ;
#>  D= 3.313e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=129 ;
#>  D= 3.254e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=130 ;
#>  D= 3.196e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=131 ;
#>  D= 3.14e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=132 ;
#>  D= 3.086e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=133 ;
#>  D= 3.033e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=134 ;
#>  D= 2.981e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=135 ;
#>  D= 2.931e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=136 ;
#>  D= 2.882e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=137 ;
#>  D= 2.834e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=138 ;
#>  D= 2.787e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=139 ;
#>  D= 2.742e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=140 ;
#>  D= 2.697e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=141 ;
#>  D= 2.654e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=142 ;
#>  D= 2.611e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=143 ;
#>  D= 2.57e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=144 ;
#>  D= 2.53e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=145 ;
#>  D= 2.49e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=146 ;
#>  D= 2.452e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=147 ;
#>  D= 2.414e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=148 ;
#>  D= 2.377e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=149 ;
#>  D= 2.341e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=150 ;
#>  D= 2.306e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=151 ;
#>  D= 2.272e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=152 ;
#>  D= 2.238e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=153 ;
#>  D= 2.205e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=154 ;
#>  D= 2.173e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=155 ;
#>  D= 2.142e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=156 ;
#>  D= 2.111e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=157 ;
#>  D= 2.081e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=158 ;
#>  D= 2.051e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=159 ;
#>  D= 2.022e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=160 ;
#>  D= 1.994e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=161 ;
#>  D= 1.966e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=162 ;
#>  D= 1.939e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=163 ;
#>  D= 1.913e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=164 ;
#>  D= 1.887e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=165 ;
#>  D= 1.861e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=166 ;
#>  D= 1.836e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=167 ;
#>  D= 1.812e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=168 ;
#>  D= 1.788e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=169 ;
#>  D= 1.764e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=170 ;
#>  D= 1.741e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=171 ;
#>  D= 1.719e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=172 ;
#>  D= 1.697e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=173 ;
#>  D= 1.675e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=174 ;
#>  D= 1.653e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=175 ;
#>  D= 1.633e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=176 ;
#>  D= 1.612e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=177 ;
#>  D= 1.592e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=178 ;
#>  D= 1.572e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=179 ;
#>  D= 1.553e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=180 ;
#>  D= 1.534e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=181 ;
#>  D= 1.515e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=182 ;
#>  D= 1.497e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=183 ;
#>  D= 1.479e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=184 ;
#>  D= 1.461e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=185 ;
#>  D= 1.443e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=186 ;
#>  D= 1.426e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=187 ;
#>  D= 1.41e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=188 ;
#>  D= 1.393e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=189 ;
#>  D= 1.377e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=190 ;
#>  D= 1.361e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=191 ;
#>  D= 1.345e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=192 ;
#>  D= 1.33e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=193 ;
#>  D= 1.315e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=194 ;
#>  D= 1.3e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=195 ;
#>  D= 1.285e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=196 ;
#>  D= 1.271e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=197 ;
#>  D= 1.257e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=198 ;
#>  D= 1.243e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=199 ;
#>  D= 1.229e-12 
#>  p= -0.5 
#> [1] "Numerical Hessian"
#> Iter=200 ;
#>  D= 1.216e-12 
#>  p= -0.5 
#> $par
#> [1] -0.4999995
#> 
#> $iterations
#> [1] 200
#> 
#> $method
#> [1] "NR"
#> 
#> $gradient
#> [1] 1.215822e-12
#> 
#> $iH
#>           [,1]
#> [1,] -2472.735
#> attr(,"det")
#>               [,1]
#> [1,] -0.0004044106
#> 

```
