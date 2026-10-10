# Convert pdf to raster format

Convert PDF file to print quality png (default 300 dpi)

## Usage

``` r
pdfconvert(
  files,
  dpi = 300,
  resolution = 1024,
  gs,
  gsopt,
  resize,
  format = "png",
  ...
)
```

## Arguments

- files:

  Vector of (pdf-)filenames to process

- dpi:

  DPI

- resolution:

  Resolution of raster image file

- gs:

  Optional ghostscript command

- gsopt:

  Optional ghostscript arguments

- resize:

  Optional resize arguments (mogrify)

- format:

  Raster format (e.g. png, jpg, tif, ...)

- ...:

  Additional arguments

## Details

Access to ghostscript program 'gs' is needed

## See also

`dev.copy2pdf`, `printdev`

## Author

Klaus K. Holst
