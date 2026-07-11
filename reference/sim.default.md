# Monte Carlo simulation

Applies a function repeatedly for a specified number of replications or
over a list/data.frame with plot and summary methods for summarizing the
Monte Carlo experiment. Can be parallelized via the future package (use
the
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
function).

## Usage

``` r
# Default S3 method
sim(
  x = NULL,
  R = 100,
  f = NULL,
  colnames = NULL,
  seed = NULL,
  args = list(),
  iter = FALSE,
  mc.cores,
  progressr.message = NULL,
  estimate.index = 1:2,
  ...
)
```

## Arguments

- x:

  function or 'sim' object

- R:

  Number of replications or data.frame with parameters

- f:

  Optional function (i.e., if x is a matrix)

- colnames:

  Optional column names

- seed:

  (optional) Seed (needed with cl=TRUE)

- args:

  (optional) list of named arguments passed to (mc)mapply

- iter:

  If TRUE the iteration number is passed as first argument to (mc)mapply

- mc.cores:

  Optional number of cores. Will use parallel::mcmapply instead of
  future

- progressr.message:

  Optional message for the progressr progress-bar

- estimate.index:

  If return object inherits from `estimate` then only these column
  indices are extracted (estimate, se, lower, upper, p-val)

- ...:

  Additional arguments to
  [`future.apply::future_mapply()`](https://future.apply.futureverse.org/reference/future_mapply.html)

## Details

To parallelize the calculation use the
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
function (e.g., `future::plan(multisession())` to distribute the
calculations over the `R` replications on all available cores). The
output is controlled via the progressr package (e.g.,
`progressr::handlers(global=TRUE)` to enable progress information).

## See also

[`summary.sim()`](https://kkholst.github.io/lava/reference/summary.sim.md)
[`plot.sim()`](https://kkholst.github.io/lava/reference/plot.sim.md)
[`sim.lvm()`](https://kkholst.github.io/lava/reference/sim.lvm.md)

## Examples

``` r
m <- lvm(y~x+e)
distribution(m,~y) <- 0
distribution(m,~x) <- uniform.lvm(a=-1.1,b=1.1)
transform(m,e~x) <- function(x) (1*x^4)*rnorm(length(x),sd=1)

onerun <- function(iter=NULL,...,n=2e3,b0=1,idx=2) {
    d <- sim(m,n,p=c("y~x"=b0))
    l <- lm(y~x,d)
    res <- c(coef(summary(l))[idx,1:2],
             confint(l)[idx,],
             coef(estimate(l), mat=TRUE)[idx,2:4])
    names(res) <- c("Estimate","Model.se","Model.lo","Model.hi",
                    "Sandwich.se","Sandwich.lo","Sandwich.hi")
    res
}
val <- sim(onerun,R=10,b0=1)
val
#>    Estimate Model.se Model.lo Model.hi Sandwich.se Sandwich.lo Sandwich.hi
#> 1  0.988453 0.033212 0.923320 1.053586 0.046984    0.896366    1.080541   
#> 2  1.013649 0.009959 0.994118 1.033180 0.014053    0.986105    1.041193   
#> 3  0.999245 0.001676 0.995959 1.002531 0.002348    0.994643    1.003847   
#> 4  0.995402 0.009091 0.977573 1.013230 0.012772    0.970370    1.020434   
#> 5  1.005428 0.003141 0.999267 1.011588 0.004419    0.996766    1.014089   
#> 6  1.025801 0.013000 1.000306 1.051296 0.018269    0.989995    1.061607   
#> 7  1.003724 0.005945 0.992066 1.015383 0.008351    0.987356    1.020093   
#> 8  1.003361 0.012299 0.979241 1.027481 0.017464    0.969131    1.037591   
#> 9  1.020572 0.011301 0.998409 1.042735 0.015674    0.989851    1.051293   
#> 10 1.012841 0.020191 0.973244 1.052439 0.028024    0.957916    1.067767   
#> 
#>      Estimate  Model.se Model.lo Model.hi Sandwich.se Sandwich.lo Sandwich.hi
#> Mean 1.006848 0.0119814 0.983350 1.030345    0.016836    0.973850    1.039845
#> SD   0.011455 0.0091467 0.023273 0.019088    0.012902    0.029998    0.025328

val <- sim(val,R=40,b0=1) ## append results
summary(val,estimate=c(1,1),confint=c(3,4,6,7),true=c(1,1))
#> 50 replications                  Time: 0.984s
#> 
#>           Estimate Estimate.1
#> Mean     1.0051843  1.0051843
#> SD       0.0139658  0.0139658
#> Coverage 0.8600000  0.9800000
#> Length   0.0399752  0.0562109
#>                              
#> Min      0.9637795  0.9637795
#> 2.5%     0.9810376  0.9810376
#> 50%      1.0030691  1.0030691
#> 97.5%    1.0340708  1.0340708
#> Max      1.0352231  1.0352231
#>                              
#> Missing  0.0000000  0.0000000
#>                              
#> True     1.0000000  1.0000000
#> Bias     0.0051843  0.0051843
#> RMSE     0.0148970  0.0148970
#> 

summary(val,estimate=c(1,1),se=c(2,5),names=c("Model","Sandwich"))
#> 50 replications                  Time: 0.984s
#> 
#>            Model Sandwich
#> Mean    1.005184 1.005184
#> SD      0.013966 0.013966
#> SE      0.010192 0.014340
#> SE/SD   0.729766 1.026779
#>                          
#> Min     0.963780 0.963780
#> 2.5%    0.981038 0.981038
#> 50%     1.003069 1.003069
#> 97.5%   1.034071 1.034071
#> Max     1.035223 1.035223
#>                          
#> Missing 0.000000 0.000000
#> 
summary(val,estimate=c(1,1),se=c(2,5),true=c(1,1),
        names=c("Model","Sandwich"),confint=TRUE)
#> 50 replications                  Time: 0.984s
#> 
#>              Model  Sandwich
#> Mean     1.0051843 1.0051843
#> SD       0.0139658 0.0139658
#> SE       0.0101918 0.0143398
#> SE/SD    0.7297662 1.0267787
#> Coverage 0.8600000 0.9800000
#> Length   0.0399510 0.0562109
#>                             
#> Min      0.9637795 0.9637795
#> 2.5%     0.9810376 0.9810376
#> 50%      1.0030691 1.0030691
#> 97.5%    1.0340708 1.0340708
#> Max      1.0352231 1.0352231
#>                             
#> Missing  0.0000000 0.0000000
#>                             
#> True     1.0000000 1.0000000
#> Bias     0.0051843 0.0051843
#> RMSE     0.0148970 0.0148970
#> 

if (interactive()) {
    plot(val,estimate=1,c(2,5),true=1,
         names=c("Model","Sandwich"),polygon=FALSE)
    plot(val,estimate=c(1,1),se=c(2,5),main=NULL,
         true=c(1,1),names=c("Model","Sandwich"),
         line.lwd=1,col=c("gray20","gray60"),
         rug=FALSE)
    plot(val,estimate=c(1,1),se=c(2,5),true=c(1,1),
         names=c("Model","Sandwich"))
}

f <- function(a=1, b=1) {
  rep(a*b, 5)
}
R <- Expand(a=1:3, b=1:3)
sim(f, R)
#>   [,1] [,2] [,3] [,4] [,5]
#> 1 1    1    1    1    1   
#> 2 2    2    2    2    2   
#> 3 3    3    3    3    3   
#> 4 2    2    2    2    2   
#> 5 4    4    4    4    4   
#> 6 6    6    6    6    6   
#> 7 3    3    3    3    3   
#> 8 6    6    6    6    6   
#> 9 9    9    9    9    9   
#> 
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> Mean 4.0000 4.0000 4.0000 4.0000 4.0000
#> SD   2.5495 2.5495 2.5495 2.5495 2.5495
sim(function(a,b) f(a,b), 3, args=c(a=5,b=5))
#>   [,1] [,2] [,3] [,4] [,5]
#> 1 25   25   25   25   25  
#> 2 25   25   25   25   25  
#> 3 25   25   25   25   25  
#> 
#>      [,1] [,2] [,3] [,4] [,5]
#> Mean   25   25   25   25   25
#> SD      0    0    0    0    0
sim(function(iter=1,a=5,b=5) iter*f(a,b), iter=TRUE, R=5)
#>   [,1] [,2] [,3] [,4] [,5]
#> 1  25   25   25   25   25 
#> 2  50   50   50   50   50 
#> 3  75   75   75   75   75 
#> 4 100  100  100  100  100 
#> 5 125  125  125  125  125 
#> 
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> Mean 75.000 75.000 75.000 75.000 75.000
#> SD   39.528 39.528 39.528 39.528 39.528

## Returning estimate objects with extra per-iteration information
onerun2 <- function(...) {
  x <- rnorm(50)
  y <- 0.5*x + rnorm(50)
  e <- estimate(lm(y ~ x))
  c(e, converged = 1, niter = sample(5:20, 1), drop.ic = TRUE)
  # works identically with summary.estimate objects
  # c(summary(e), converged = 1, niter = sample(5:20, 1), drop.ic = TRUE)
}
val2 <- sim(onerun2, R = 10)
val2
#>    (Intercept) x        (Intercept).Std.Err x.Std.Err converged niter   
#> 1   0.10011     0.48755  0.13145             0.09815   1.00000   6.00000
#> 2   0.20247     0.83953  0.15051             0.20878   1.00000  19.00000
#> 3   0.11318     0.48145  0.12845             0.11748   1.00000  13.00000
#> 4  -0.24943     0.56391  0.15076             0.16088   1.00000   5.00000
#> 5  -0.23902     0.56304  0.14119             0.13798   1.00000  16.00000
#> 6   0.24264     0.60289  0.13804             0.14833   1.00000  13.00000
#> 7   0.06325     0.46903  0.11287             0.13361   1.00000  10.00000
#> 8  -0.18888     0.47657  0.12326             0.09038   1.00000  10.00000
#> 9  -0.08492     0.33604  0.13019             0.14964   1.00000  12.00000
#> 10  0.13528     0.54151  0.12871             0.10612   1.00000  13.00000
#> 
#>       (Intercept)       x
#> Mean    0.0094681 0.53615
#> SD      0.1844915 0.12973
#> SE      0.1335432 0.13513
#> SE/SD   0.7238451 1.04166
# extra columns are stored but excluded from default summary statistics
idx <- c("(Intercept)", "x", "niter")
summary(val2, estimate=idx, true=c(0,0.5,NA))
#> 10 replications                  Time: 0.201s
#> 
#>          (Intercept)        x   niter
#> Mean       0.0094681 0.536154 11.7000
#> SD         0.1844915 0.129729  4.2177
#> SE         0.1335432 0.135133        
#> SE/SD      0.7238451 1.041659        
#> Coverage   1.0000000 1.000000        
#> Length     0.5234799 0.529713        
#>                                      
#> Min       -0.2494253 0.336044  5.0000
#> 2.5%      -0.2470833 0.365967  5.2250
#> 50%        0.0816777 0.514531 12.5000
#> 97.5%      0.2336029 0.786287 18.3250
#> Max        0.2426428 0.839530 19.0000
#>                                      
#> Missing    0.0000000 0.000000  0.0000
#>                                      
#> True       0.0000000 0.500000        
#> Bias       0.0094681 0.036154        
#> RMSE       0.1847343 0.134673        
#> 
```
