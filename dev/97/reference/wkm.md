# Weighted K-means

Weighted K-means via Lloyd's algorithm

## Usage

``` r
wkm(
  x,
  mu,
  data,
  weights = rep(1, NROW(x)),
  iter.max = 20,
  n.start = 5,
  init = "kmpp",
  ...
)
```

## Arguments

- x:

  Data (or formula)

- mu:

  Initial centers (or number centers chosen randomly among x)

- data:

  optional data frmae

- weights:

  Optional weights

- iter.max:

  Max number of iterations

- n.start:

  Number of restarts

- init:

  method to create initial centres (default kmeans++)

- ...:

  Additional arguments to lower level functions

## Author

Klaus K. Holst

## Examples

``` r
## Two well-separated Gaussian blobs in 2-D
set.seed(1)
x <- rbind(matrix(rnorm(100, mean = -3), ncol = 2),
           matrix(rnorm(100, mean =  3), ncol = 2))
res <- wkm(x, mu = 2)
table(res$cluster)
#> 
#>  1  2 
#> 50 50 
res$center
#> class: 1
#> [1] 2.847515 3.076869
#> ------------------------------------------------------------ 
#> class: 2
#> [1] -2.899552 -2.882674

## Supply explicit initial centers (as a list)
res2 <- wkm(x, mu = list(c(-3, -3), c(3, 3)))

## Weighted clustering: up-weight the second blob
w <- c(rep(1, 50), rep(10, 50))
res3 <- wkm(x, mu = 2, weights = w)

## Formula interface on a data.frame
wkm(~ Sepal.Length + Sepal.Width, data = iris, mu = 3)
#> $cluster
#>   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
#>   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3 
#>  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40 
#>   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3 
#>  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60 
#>   3   3   3   3   3   3   3   3   3   3   2   2   2   1   2   1   2   1   2   1 
#>  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79  80 
#>   1   1   1   1   1   2   1   1   1   1   1   1   1   1   2   2   2   2   1   1 
#>  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99 100 
#>   1   1   1   1   1   1   2   1   1   1   1   1   1   1   1   1   1   1   1   1 
#> 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 
#>   2   1   2   2   2   2   1   2   2   2   2   2   2   1   1   2   2   2   2   1 
#> 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 
#>   2   1   2   1   2   2   1   1   2   2   2   2   2   1   1   2   2   2   1   2 
#> 141 142 143 144 145 146 147 148 149 150 
#>   2   2   1   2   2   2   1   2   2   1 
#> 
#> $center
#> class: 1
#>  (Intercept) Sepal.Length  Sepal.Width 
#>     1.000000     5.773585     2.692453 
#> ------------------------------------------------------------ 
#> class: 2
#>  (Intercept) Sepal.Length  Sepal.Width 
#>     1.000000     6.812766     3.074468 
#> ------------------------------------------------------------ 
#> class: 3
#>  (Intercept) Sepal.Length  Sepal.Width 
#>        1.000        5.006        3.428 
#> 
#> $ssw
#> [1] 11.3000 12.6217 13.1290
#> 
```
