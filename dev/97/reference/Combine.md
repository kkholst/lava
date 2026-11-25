# Report estimates across different models

Report estimates across different models

## Usage

``` r
Combine(x, ...)
```

## Arguments

- x:

  list of model objects

- ...:

  additional arguments to lower-level functions

## Author

Klaus K. Holst

## Examples

``` r
data(serotonin)
m1 <- lm(cau ~ age*gene1 + age*gene2,data=serotonin)
m2 <- lm(cau ~ age + gene1,data=serotonin)
m3 <- lm(cau ~ age*gene2,data=serotonin)

Combine(list(A=m1,B=m2,C=m3),fun=function(x)
     c("_____"="",R2=" "%++%format(summary(x)$r.squared,digits=2)))
#>             A                            B                          
#> (Intercept)  0.88  [0.85;0.91]   p<0.001  0.89  [0.86;0.91]  p<0.001
#> age          0.02  [-0.01;0.05]  p=0.108 -0.01  [-0.02;0.01] p=0.414
#> gene1       -0.02  [-0.06;0.01]  p=0.228 -0.02  [-0.06;0.01] p=0.219
#> gene2        0     [-0.03;0.04]  p=0.93                             
#> age:gene1   -0.04  [-0.07;-0.01] p=0.019                            
#> age:gene2   -0.02  [-0.05;0.01]  p=0.218                            
#> _____                                                               
#> R2           0.039                        0.0083                    
#>             C                          
#> (Intercept)  0.88  [0.85;0.9]   p<0.001
#> age          0.01  [-0.02;0.03] p=0.549
#> gene1                                  
#> gene2        0     [-0.03;0.04] p=0.99 
#> age:gene1                              
#> age:gene2   -0.03  [-0.06;0.01] p=0.126
#> _____                                  
#> R2           0.012                     
```
