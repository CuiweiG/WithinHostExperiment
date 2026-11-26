# Track variant frequency trajectories over time

Extracts per-variant allele frequency trajectories across longitudinal
time points within each host. Returns a long-format `DataFrame` suitable
for downstream analysis or plotting.

## Usage

``` r
trackFrequency(
  whe,
  hostCol = "host_id",
  timeCol = "timepoint",
  usePassedOnly = TRUE
)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md) with longitudinal
  samples from the same host.

- hostCol:

  Character scalar. Column name in `colData` containing the host/patient
  identifier (default: `"host_id"`).

- timeCol:

  Character scalar. Column name in `colData` containing the time point
  (default: `"timepoint"`). It must hold numbers, dates, date-times,
  time differences or a factor whose levels are in time order; text is
  read as numbers when every value is one, and is otherwise an error, as
  in [`flagTemporalInconsistency`](flagTemporalInconsistency.md).
  Samples whose time point is `NA` are left out, with a warning.

- usePassedOnly:

  Logical. If TRUE (default), only use variants where `qcPass == TRUE`.
  With `FALSE` every call is returned and `qc_passed` marks those that
  failed.

## Value

A `DataFrame` with columns:

- host_id:

  Host identifier.

- variant_key:

  Unique variant identifier (chrom:pos:alt).

- timepoint:

  Time point value.

- frequency:

  Alternative allele frequency at that time point.

- qc_passed:

  Logical. Whether the call passed QC at that time point; an `NA`
  `qcPass` counts as not passed (`TRUE` if no QC assay present).

## Examples

``` r
gr <- GenomicRanges::GRanges("seg1",
    IRanges::IRanges(c(100, 200), width = 1))
S4Vectors::mcols(gr)$ref <- c("A", "C")
S4Vectors::mcols(gr)$alt <- c("T", "G")
freq <- matrix(c(0.05, 0.10, 0.15, 0.20), nrow = 2, ncol = 2)
whe <- WithinHostExperiment(
    assays = list(altFreq = freq),
    rowRanges = gr,
    colData = S4Vectors::DataFrame(
        sample_id = c("S1_t1", "S1_t2"),
        host_id = c("H1", "H1"),
        timepoint = c(1, 2)))
trackFrequency(whe)
#> DataFrame with 4 rows and 5 columns
#>       host_id variant_key timepoint frequency qc_passed
#>   <character> <character> <numeric> <numeric> <logical>
#> 1          H1  seg1:100:T         1      0.05      TRUE
#> 2          H1  seg1:200:G         1      0.10      TRUE
#> 3          H1  seg1:100:T         2      0.15      TRUE
#> 4          H1  seg1:200:G         2      0.20      TRUE
```
