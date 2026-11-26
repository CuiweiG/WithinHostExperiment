# Pairwise caller concordance matrix

For a list of [`WithinHostExperiment`](WithinHostExperiment.md) objects
(one per variant caller), computes pairwise concordance statistics
including shared sites, caller-unique sites, and Jaccard index.

## Usage

``` r
calcCallerConcordance(whe_list)
```

## Arguments

- whe_list:

  A named list of [`WithinHostExperiment`](WithinHostExperiment.md)
  objects, one per variant caller.

## Value

A `DataFrame` with columns:

- caller1:

  First caller name.

- caller2:

  Second caller name.

- n_shared:

  Number of sites detected by both callers.

- n_only1:

  Sites unique to caller1.

- n_only2:

  Sites unique to caller2.

- jaccard:

  Jaccard similarity index (shared / union).

A caller detects a site when some sample has a frequency above 0 there;
a recorded 0 is a site assessed without the allele.

## References

Roder AE et al. (2023). Optimized quantification of intra-host viral
diversity in SARS-CoV-2 and influenza virus sequence data. *mBio*
14:e01046-23.
[doi:10.1128/mbio.01046-23](https://doi.org/10.1128/mbio.01046-23)

## Examples

``` r
gr1 <- GenomicRanges::GRanges("seg1",
    IRanges::IRanges(c(100, 200, 300), width = 1))
S4Vectors::mcols(gr1)$ref <- c("A", "C", "G")
S4Vectors::mcols(gr1)$alt <- c("T", "G", "A")
whe1 <- WithinHostExperiment(
    assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
    rowRanges = gr1,
    colData = S4Vectors::DataFrame(sample_id = "S1"))
gr2 <- GenomicRanges::GRanges("seg1",
    IRanges::IRanges(c(100, 200, 400), width = 1))
S4Vectors::mcols(gr2)$ref <- c("A", "C", "T")
S4Vectors::mcols(gr2)$alt <- c("T", "G", "C")
whe2 <- WithinHostExperiment(
    assays = list(altFreq = matrix(c(0.12, 0.22, 0.15), ncol = 1)),
    rowRanges = gr2,
    colData = S4Vectors::DataFrame(sample_id = "S1"))
calcCallerConcordance(list(ivar = whe1, lofreq = whe2))
#> DataFrame with 1 row and 6 columns
#>       caller1     caller2  n_shared   n_only1   n_only2   jaccard
#>   <character> <character> <integer> <integer> <integer> <numeric>
#> 1        ivar      lofreq         2         1         1       0.5
```
