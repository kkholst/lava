# Expand long-format data

Expand long-format data to all subject and visit combinations

## Usage

``` r
complete_visits(data, id, time, response = NULL, levels = NULL)
```

## Arguments

- data:

  (data.frame) data in long format

- id:

  (character) name of the subject variable

- time:

  (character) name of the visit variable. In the returned data this
  variable is a factor with levels `levels`.

- response:

  (character) optional name of the response variable which is set to
  `NA` in added rows

- levels:

  (vector) visit levels. Defaults to the levels of the visit variable if
  it is a factor, and otherwise its sorted unique values.

## Value

data.frame with `length(levels)` rows per subject

## Details

Returns a long-format data.frame with exactly one row per subject and
visit (`levels`), sorted by subject (in order of first appearance) and
visit. Rows that are not present in `data` are added by copying all
columns from the previous existing row of the subject (last observation
carried forward), or from the next existing row if there is no previous
row. The response of added rows is set to `NA`.

## Examples

``` r
d <- data.frame(id = c(1, 1, 2), visit = c(1, 2, 2), y = 1:3, x = 4:6)
complete_visits(d, "id", "visit", "y")
#>   id visit  y x
#> 1  1     1  1 4
#> 2  1     2  2 5
#> 3  2     1 NA 6
#> 4  2     2  3 6
```
