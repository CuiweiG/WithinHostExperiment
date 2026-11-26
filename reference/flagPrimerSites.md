# Flag iSNVs in primer binding regions

For amplicon-based sequencing data, marks variants that fall within
primer binding sites in the `qcPass` assay.

## Usage

``` r
flagPrimerSites(whe, primers, action = c("flag", "remove"))
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md).

- primers:

  `GRanges`. Primer binding site coordinates. Can be loaded from a BED
  file using
  [`rtracklayer::import.bed()`](https://rdrr.io/pkg/rtracklayer/man/BEDFile-class.html).
  Their sequence names must match those of the sites: when the two share
  none, the function stops and prints both, rather than report that no
  site lies in a primer.

- action:

  Character. `"flag"` (default) sets qcPass to FALSE for affected sites,
  preserving all data for audit. `"remove"` is provided as an explicit
  convenience for users who require destructive filtering after
  inspection.

## Value

A [`WithinHostExperiment`](WithinHostExperiment.md).

## Examples

``` r
library(GenomicRanges)
vcf <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(vcf,
    colData = S4Vectors::DataFrame(sample_id = "donor"),
    caller = "ivar")
primers <- GRanges("test_segment", IRanges::IRanges(
    start = c(40, 880), width = c(25, 30)))
whe <- flagPrimerSites(whe, primers)
qcLog(whe)
#> DataFrame with 1 row and 6 columns
#>            step        parameter       value n_flagged  n_passed
#>     <character>      <character> <character> <integer> <integer>
#> 1 primer_filter n_primer_regions           2         2        13
#>             timestamp
#>           <character>
#> 1 2026-11-26T15:37:52
```
