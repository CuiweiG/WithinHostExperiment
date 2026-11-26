# Multi-caller consensus iSNV detection

Takes a named list of [`WithinHostExperiment`](WithinHostExperiment.md)
objects (one per variant caller) and identifies variant sites detected
by multiple callers. Sites are matched by genomic position (seqnames +
start) and alternative allele. Only sites detected by at least
`min_callers` with a frequency spread no greater than `freq_tolerance`
are retained.

## Usage

``` r
consensusISNV(whe_list, min_callers = 2L, freq_tolerance = 0.05)
```

## Arguments

- whe_list:

  A named list of [`WithinHostExperiment`](WithinHostExperiment.md)
  objects, one per variant caller (e.g.,
  `list(ivar = whe1, lofreq = whe2)`).

- min_callers:

  Integer scalar. Minimum number of callers that must detect a site
  (default: 2L).

- freq_tolerance:

  Numeric scalar. Maximum allowed range (max - min) of allele
  frequencies across callers for a site to pass (default: 0.05); a range
  equal to it, to within rounding error, passes.

## Value

A [`WithinHostExperiment`](WithinHostExperiment.md) containing only
consensus sites, with a column for every sample of the inputs.
`mcols(rowRanges)` includes:

- n_callers:

  Number of callers detecting the site; a caller that lists the same
  site and allele twice counts once.

- caller_names:

  Comma-separated caller names.

- freq_spread:

  Range of the allele frequencies of the calls, across callers.

Where a site passes in more than one sample, these three columns
describe one of those samples, not all of them.

The `altFreq` assay contains the mean frequency of the calls for each
consensus site and `NA` where the site did not pass in that sample, so a
sample with no consensus site is a column of `NA`. Every call enters the
mean and the spread, so a caller that lists a site twice is weighted
twice there, though it counts once in `n_callers`. A sample held by
fewer than `min_callers` inputs cannot reach consensus; it is kept as a
column of `NA`, with a warning naming it, since that column would
otherwise read as a sample without iSNVs. The constructor adds a
`qcPass` assay that is `TRUE` throughout. Read counts, depths and the QC
history are not carried over, so run [`flagISNV`](flagISNV.md) on the
individual callers' objects rather than on the consensus.

The [`seqinfo()`](https://rdrr.io/pkg/Seqinfo/man/seqinfo.html) of the
result merges those of the inputs; a sequence whose length, circularity
or genome differs between inputs is an error. Each sample's `colData`
row comes from the first input that holds the sample, so host, role and
transmission-pair columns carry over.

## Details

When technical replicates are unavailable, using a multi-caller
consensus is recommended for reducing false positive iSNV calls (Roder
et al. 2023).

## References

Roder AE et al. (2023). Optimized quantification of intra-host viral
diversity in SARS-CoV-2 and influenza virus sequence data. *mBio*
14:e01046-23.
[doi:10.1128/mbio.01046-23](https://doi.org/10.1128/mbio.01046-23)

## Examples

``` r
gr <- GenomicRanges::GRanges("seg1",
    IRanges::IRanges(c(100, 200, 300), width = 1))
S4Vectors::mcols(gr)$ref <- c("A", "C", "G")
S4Vectors::mcols(gr)$alt <- c("T", "G", "A")
whe1 <- WithinHostExperiment(
    assays = list(altFreq = matrix(c(0.10, 0.20, 0.30), ncol = 1)),
    rowRanges = gr, colData = S4Vectors::DataFrame(sample_id = "S1"))
gr2 <- GenomicRanges::GRanges("seg1",
    IRanges::IRanges(c(100, 200, 400), width = 1))
S4Vectors::mcols(gr2)$ref <- c("A", "C", "T")
S4Vectors::mcols(gr2)$alt <- c("T", "G", "C")
whe2 <- WithinHostExperiment(
    assays = list(altFreq = matrix(c(0.12, 0.18, 0.15), ncol = 1)),
    rowRanges = gr2, colData = S4Vectors::DataFrame(sample_id = "S1"))
consensus <- consensusISNV(list(caller1 = whe1, caller2 = whe2))
consensus
#> WithinHostExperiment with 2 variant sites and 1 samples
#>   assays(2): altFreq, qcPass
#>   qcLog: 0 filtering step(s)
#>   transmission pairs: 0 
```
