# Plot donor-recipient frequency scatter

Creates a publication-quality scatter plot comparing variant frequencies
between a donor and recipient, highlighting shared and unique variants
with distinct shapes and colours.

## Usage

``` r
plotPairScatter(whe, pairId, threshold = 0.03)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md) with transmission
  pairs.

- pairId:

  Character. Which pair to plot.

- threshold:

  Numeric. Variant calling threshold shown as reference lines (default:
  0.03).

## Value

A ggplot object.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
    vcf1 <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
        package = "WithinHostExperiment")
    whe <- readWithinHost(c(vcf1, vcf2),
        colData = S4Vectors::DataFrame(
            sample_id = c("donor", "recipient"),
            role = c("donor", "recipient"),
            pair_id = c("p1", "p1")),
        caller = "ivar")
    plotPairScatter(whe, pairId = "p1")
}
```
