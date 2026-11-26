# Plot iSNV frequency spectrum

Creates a publication-quality histogram or density plot of alternative
allele frequencies, optionally faceted by sample or metadata group. Only
calls are drawn: a recorded frequency of 0 is a site assessed without
the allele and is left out. A site at frequency 1 is drawn unless QC has
removed it.

## Usage

``` r
plotFrequencySpectrum(
  whe,
  type = c("histogram", "density"),
  facetBy = NULL,
  usePassedOnly = TRUE,
  threshold = NULL,
  freqRange = NULL
)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md).

- type:

  Character. `"histogram"` (default) or `"density"`.

- facetBy:

  Character (optional). Column name in `colData` to facet by.

- usePassedOnly:

  Logical. Only show QC-passed variants (default: TRUE).

- threshold:

  Numeric scalar (optional). Variant calling threshold to display as
  reference line. Default `NULL` (no line shown). Set to e.g. 0.03 for a
  3% threshold.

- freqRange:

  Numeric vector of length 2 (optional). Frequency range for the x-axis
  as proportions, e.g. `c(0, 0.5)`. Default `NULL` uses the data range.

## Value

A ggplot object.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
    vcf <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "d"),
        caller = "ivar")
    plotFrequencySpectrum(whe, threshold = 0.03)
}
```
