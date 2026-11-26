# Detect emerging variants across time points

Identifies variants whose allele frequency crosses a threshold between
consecutive time points within each host. A variant is considered
emerging if it was below `freqThreshold` at one time point and increased
by at least `minIncrease` at the next. Replicate samples taken at one
time point are combined first: the variant's frequency there is the mean
of the replicates' calls, so replicates never form an interval of their
own and the result does not depend on the order of the columns.
Intervals join the consecutive time points at which the variant has a
call in the host; a time point without a call for it is passed over,
since no call is not a frequency of 0, so a variant first called at a
later time point is not reported.

## Usage

``` r
detectEmergingVariants(
  whe,
  hostCol = "host_id",
  timeCol = "timepoint",
  freqThreshold = 0.05,
  minIncrease = 0.05,
  usePassedOnly = TRUE
)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md) with longitudinal
  samples.

- hostCol:

  Character scalar. Column name in `colData` for host identifier
  (default: `"host_id"`). Samples whose host is `NA` are left out, with
  a warning.

- timeCol:

  Character scalar. Column name in `colData` for time point (default:
  `"timepoint"`), read as in [`trackFrequency`](trackFrequency.md): it
  is put in time order, and samples whose time point is `NA` are left
  out, with a warning.

- freqThreshold:

  Numeric scalar. Frequency threshold for emergence (default: 0.05).

- minIncrease:

  Numeric scalar. Minimum frequency increase between consecutive time
  points (default: 0.05).

- usePassedOnly:

  Logical. If TRUE (default), an interval is reported only when the
  later time point has a call that passed QC, and its frequency is the
  mean of those calls. The earlier time point is the mean of its passing
  calls too, or of all its calls when none passed:
  [`flagISNV`](flagISNV.md) fails a frequency below the filter's
  `minFreq`, and such a frequency is the low start that emergence is
  measured from. `qc_from` says which. Intervals are formed on every
  call before QC is applied, so a time point that failed QC, for example
  above `maxFreq`, is never passed over. With `FALSE` every call is
  used.

## Value

A `DataFrame` with columns:

- host_id:

  Host identifier.

- variant_key:

  Variant identifier (chrom:pos:alt).

- time_from:

  Earlier time point.

- time_to:

  Later time point.

- freq_from:

  Frequency at the earlier time point.

- freq_to:

  Frequency at the later time point.

- delta_freq:

  Change in frequency.

- qc_from:

  Logical. Whether a call at the earlier time point passed QC (`TRUE`
  without a `qcPass` assay). A `FALSE` start failed QC, which a
  frequency below `minFreq` does by design; one that failed for low
  depth may be set aside.

One row is one interval, between two distinct time points. Rows are
ordered by host, then by variant, each in order of first appearance in
the object, then by time point.

## Examples

``` r
gr <- GenomicRanges::GRanges("seg1",
    IRanges::IRanges(c(100, 200), width = 1))
S4Vectors::mcols(gr)$ref <- c("A", "C")
S4Vectors::mcols(gr)$alt <- c("T", "G")
freq <- matrix(c(0.02, 0.10, 0.12, 0.15), nrow = 2, ncol = 2)
whe <- WithinHostExperiment(
    assays = list(altFreq = freq),
    rowRanges = gr,
    colData = S4Vectors::DataFrame(
        sample_id = c("S1_t1", "S1_t2"),
        host_id = c("H1", "H1"),
        timepoint = c(1, 2)))
detectEmergingVariants(whe)
#> DataFrame with 1 row and 8 columns
#>       host_id variant_key time_from   time_to freq_from   freq_to delta_freq
#>   <character> <character> <numeric> <numeric> <numeric> <numeric>  <numeric>
#> 1          H1  seg1:100:T         1         2      0.02      0.12        0.1
#>     qc_from
#>   <logical>
#> 1      TRUE
```
