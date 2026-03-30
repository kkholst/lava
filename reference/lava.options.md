# Set global options for `lava`

Extract and set global parameters of `lava`. In particular optimization
parameters for the `estimate` function.

## Usage

``` r
lava.options(...)
```

## Arguments

- ...:

  Arguments

## Value

`list` of parameters

## Details

- `param`: 'relative' (factor loading and variance of one endogenous
  variables in each measurement model are fixed to one), 'absolute'
  (mean and variance of latent variables are set to 0 and 1,
  respectively), 'hybrid' (intercept of latent variables is fixed to 0,
  and factor loading of at least one endogenous variable in each
  measurement model is fixed to 1), 'none' (no constraints are added)

- `layout`: One of 'dot','fdp','circo','twopi','neato','osage'

- `messages`: Set to 0 to disable various output messages

- ...

see `control` parameter of the `estimate` function.

## Author

Klaus K. Holst

## Examples

``` r
if (FALSE) { # \dontrun{
lava.options(iter.max=100,messages=0)
} # }
```
