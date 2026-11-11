# Apply a Function to a Data Frame Split by Factors

Apply a Function to a Data Frame Split by Factors

## Usage

``` r
By(x, INDICES, FUN, COLUMNS, array = FALSE, ...)
```

## Arguments

- x:

  Data frame

- INDICES:

  Indices (vector or list of indices, vector of column names, or formula
  of column names)

- FUN:

  A function to be applied to data frame subsets of 'data'.

- COLUMNS:

  (Optional) subset of columns of x to work on

- array:

  if TRUE an array/matrix is always returned

- ...:

  Additional arguments to lower-level functions

## Details

Simple wrapper of the 'by' function

## Author

Klaus K. Holst

## Examples

``` r
By(datasets::CO2,~Treatment+Type,colMeans,~conc)
#>             Type
#> Treatment    Quebec Mississippi
#>   nonchilled    435         435
#>   chilled       435         435
By(datasets::CO2,~Treatment+Type,colMeans,~conc+uptake)
#> Treatment: nonchilled
#> Type: Quebec
#>      conc    uptake 
#> 435.00000  35.33333 
#> ------------------------------------------------------------ 
#> Treatment: chilled
#> Type: Quebec
#>      conc    uptake 
#> 435.00000  31.75238 
#> ------------------------------------------------------------ 
#> Treatment: nonchilled
#> Type: Mississippi
#>      conc    uptake 
#> 435.00000  25.95238 
#> ------------------------------------------------------------ 
#> Treatment: chilled
#> Type: Mississippi
#>      conc    uptake 
#> 435.00000  15.81429 
```
