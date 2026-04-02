# Extract graph

Extract or replace graph object

## Usage

``` r
Graph(x, ...)

Graph(x, ...) <- value
```

## Arguments

- x:

  Model object

- ...:

  Additional arguments to be passed to the low level functions

- value:

  New `graphNEL` object

## See also

[`Model`](https://kkholst.github.io/lava/reference/Model.md)

## Author

Klaus K. Holst

## Examples

``` r
m <- lvm(y~x)
Graph(m)
#> NULL
```
