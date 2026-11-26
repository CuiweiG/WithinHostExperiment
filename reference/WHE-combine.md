# Combine WithinHostExperiment objects by rows or columns

Methods for the
[`combineRows`](https://rdrr.io/pkg/S4Vectors/man/RectangularData-class.html)
and
[`combineCols`](https://rdrr.io/pkg/S4Vectors/man/RectangularData-class.html)
generics. `combineRows` combines objects vertically (adding variant
sites); all objects must have the same samples (columns). Samples are
matched by their column names: an object that holds them in another
order is put in the first object's order, and one that holds other
samples is an error naming them. The `colData` and the `metadata` are
the first object's. `combineCols` combines objects horizontally (adding
samples); all objects must have the same variant sites (rows). Sites are
matched by sequence name, position, `ref` and `alt`: an object that
holds them in another order is put in the first object's order, and one
that holds other sites is an error counting them. Objects read from
separate files usually hold different sites; read the files together
with [`readWithinHost`](readWithinHost.md) instead, which takes the
union. QC logs and filters from all objects are merged.

## Usage

``` r
# S4 method for class 'WithinHostExperiment'
combineRows(x, ...)

# S4 method for class 'WithinHostExperiment'
combineCols(x, ..., use.names = TRUE)
```

## Arguments

- x:

  A [`WithinHostExperiment`](WithinHostExperiment.md).

- ...:

  Further [`WithinHostExperiment`](WithinHostExperiment.md) objects to
  combine.

- use.names:

  Accepted for compatibility with the generic. Sites are matched by
  sequence name, position and alleles, not by row name.

## Value

A [`WithinHostExperiment`](WithinHostExperiment.md).

## Examples

``` r
data(example_whe)
whe1 <- example_whe[1:10, ]
whe2 <- example_whe[11:nrow(example_whe), ]
merged <- combineRows(whe1, whe2)
nrow(merged) == nrow(example_whe)
#> [1] TRUE

whe_d <- example_whe[, 1]
whe_r <- example_whe[, 2]
merged <- combineCols(whe_d, whe_r)
ncol(merged)
#> [1] 2
```
