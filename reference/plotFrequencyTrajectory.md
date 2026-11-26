# Plot allele frequency trajectories over time

Creates a ggplot2 spaghetti plot of allele frequency over time for each
variant within a specified host. Each line represents one variant,
coloured with the Okabe-Ito palette (Wong 2011); with more than eight
variants the colours repeat.

## Usage

``` r
plotFrequencyTrajectory(
  whe,
  hostId,
  hostCol = "host_id",
  timeCol = "timepoint",
  usePassedOnly = TRUE
)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md) with longitudinal
  samples.

- hostId:

  Character scalar. Host identifier to plot.

- hostCol:

  Character scalar. Column name in `colData` for host identifier
  (default: `"host_id"`).

- timeCol:

  Character scalar. Column name in `colData` for time point (default:
  `"timepoint"`), read as in [`trackFrequency`](trackFrequency.md).
  Numbers given as text are drawn in numeric order.

- usePassedOnly:

  Logical. If TRUE (default), only use QC-passed variants.

## Value

A ggplot object.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
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
    plotFrequencyTrajectory(whe, hostId = "H1")
}
```
