# Plot a within-host site frequency spectrum

Creates a publication-quality bar plot of a
[`WithinHostSFS`](WithinHostSFS-class.md) object.

## Usage

``` r
plotSFS(sfs, normalize = FALSE)
```

## Arguments

- sfs:

  A [`WithinHostSFS`](WithinHostSFS-class.md) object, or a named list of
  `WithinHostSFS` objects for overlay.

- normalize:

  Logical. If `TRUE`, plot proportions instead of counts. Default:
  `FALSE`.

## Value

A ggplot object.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
    vcf <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    sfs <- buildSFS(whe, genomeLength = 1000L)
    plotSFS(sfs)
}
```
