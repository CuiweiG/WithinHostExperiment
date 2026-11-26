# Generate a QC summary report

Produces per-sample summary statistics from a
[`WithinHostExperiment`](WithinHostExperiment.md).

## Usage

``` r
qcSummary(whe)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md).

## Value

A `DataFrame` with one row per sample and columns:

- sample_id:

  Sample identifier.

- total_sites:

  The sample's calls: sites where its frequency is above 0. A recorded 0
  is a site assessed without the allele and is not counted.

- passed_sites, flagged_sites:

  The calls whose `qcPass` is `TRUE` and `FALSE`.

- pass_rate:

  `passed_sites / total_sites`, `NA` without calls.

- median_depth, median_freq:

  The median read depth and the median frequency of the calls.

## Examples

``` r
vcf <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(vcf,
    colData = S4Vectors::DataFrame(sample_id = "donor"),
    caller = "ivar")
whe <- flagISNV(whe, ISNVFilter(minDepth = 500L))
qcSummary(whe)
#> DataFrame with 1 row and 7 columns
#>     sample_id total_sites passed_sites flagged_sites pass_rate median_depth
#>   <character>   <integer>    <integer>     <integer> <numeric>    <numeric>
#> 1       donor          15           14             1  0.933333         1600
#>   median_freq
#>     <numeric>
#> 1        0.18
```
