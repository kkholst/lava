# fplot

Faster plot via RGL

## Usage

``` r
fplot(
  x,
  y,
  z = NULL,
  xlab,
  ylab,
  ...,
  z.col = topo.colors(64),
  data = parent.frame(),
  add = FALSE,
  aspect = c(1, 1),
  zoom = 0.8
)
```

## Arguments

- x:

  X variable

- y:

  Y variable

- z:

  Z variable (optional)

- xlab:

  x-axis label

- ylab:

  y-axis label

- ...:

  additional arggument to lower-level plot functions

- z.col:

  color (use argument alpha to set transparency)

- data:

  data.frame

- add:

  if TRUE use current active device

- aspect:

  aspect ratio

- zoom:

  zoom level

## Examples

``` r
if (interactive()) {
data(iris)
fplot(Sepal.Length ~ Petal.Length+Species, data=iris, size=2, type="s")
}
```
