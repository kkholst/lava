# Define constant risk difference or relative risk association for binary exposure

Set up model as defined in Richardson, Robins and Wang (2017).

## Usage

``` r
binomial.rd(
  x,
  response,
  exposure,
  target.model,
  nuisance.model,
  exposure.model = binomial.lvm(),
  ...
)
```

## Arguments

- x:

  model

- response:

  response variable (character or formula)

- exposure:

  exposure variable (character or formula)

- target.model:

  variable defining the linear predictor for the target model

- nuisance.model:

  variable defining the linear predictor for the nuisance model

- exposure.model:

  model for exposure (default binomial logit link)

- ...:

  additional arguments to lower level functions
