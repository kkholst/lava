# Calculate diagnostic tests for 2x2 table

Calculate prevalence, sensitivity, specificity, and positive and
negative predictive values

## Usage

``` r
diagtest(
  table,
  positive = 2,
  exact = FALSE,
  p0 = NA,
  confint = c("logit", "arcsin", "pseudoscore", "exact"),
  ...
)
```

## Arguments

- table:

  Table or (matrix/data.frame with two columns)

- positive:

  Switch reference

- exact:

  If TRUE exact binomial proportions CI/test will be used

- p0:

  Optional null hypothesis (test prevalenc, sensitivity, ...)

- confint:

  Type of confidence limits

- ...:

  Additional arguments to lower level functions

## Details

Table should be in the format with outcome in columns and test in rows.
Data.frame should be with test in the first column and outcome in the
second column.

## Author

Klaus Holst

## Examples

``` r
M <- as.table(matrix(c(42,12,
                       35,28),ncol=2,byrow=TRUE,
                     dimnames=list(rater=c("no","yes"),gold=c("no","yes"))))
diagtest(M,exact=TRUE)
#> Call: diagtest(table = M, exact = TRUE)
#> Confidence limits: exact
#> ────────────────────────────────────────────────────────────────────────────────
#>      gold
#> rater no yes        no   yes
#>   no  42  12     0.359 0.103
#>   yes 35  28     0.299 0.239
#> 
#> Positive outcome: 'yes'
#> ────────────────────────────────────────────────────────────────────────────────
#>                         Estimate    2.5%   97.5% P-value
#> Prevalence               0.34188 0.25670 0.43528        
#> Test                     0.53846 0.44389 0.63104        
#> Sensitivity              0.70000 0.53468 0.83437        
#> Specificity              0.54545 0.42790 0.65940        
#> PositivePredictiveValue  0.44444 0.31917 0.57511        
#> NegativePredictiveValue  0.77778 0.64400 0.87956        
#> Accuracy                 0.59829 0.50363 0.68785        
#> Homogeneity              0.74468 0.59650 0.86055  0.0011
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#> Prevalence:              Prob( outcome+ )
#> Test:                    Prob( test+ )
#> Sensitivity (True positive rate):    Prob( test+ | outcome+ )
#> Specificity (True negative rate):    Prob( test- | outcome- )
#> Positive predictive value (Precision):   Prob( outcome+ | test+ )
#> Negative predictive value:       Prob( outcome- | test- )
#> Accuracy:                Prob( correct classification )
#> Homogeneity/Symmetry:            Prob( outcome+, test- | discordant ), H0: p=0.5 
#> 
```
