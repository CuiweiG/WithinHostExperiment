# Benchmark variant callers against truth or each other

If a truth `WithinHostExperiment` is provided, computes true positives,
false positives, false negatives, precision, recall, and F1 for each
caller. If no truth is provided, computes pairwise agreement using
[`calcCallerConcordance`](calcCallerConcordance.md).

## Usage

``` r
benchmarkCallers(whe_list, truth = NULL)
```

## Arguments

- whe_list:

  A named list of [`WithinHostExperiment`](WithinHostExperiment.md)
  objects, one per variant caller.

- truth:

  An optional [`WithinHostExperiment`](WithinHostExperiment.md)
  containing the true variant set. Default `NULL`.

## Value

A `DataFrame`. When `truth` is provided:

- caller:

  Caller name.

- n_total:

  Total sites called by this caller.

- TP:

  True positives.

- FP:

  False positives.

- FN:

  False negatives.

- precision:

  TP / (TP + FP).

- recall:

  TP / (TP + FN).

- F1:

  Harmonic mean of precision and recall.

When `truth` is `NULL`, the result of
[`calcCallerConcordance`](calcCallerConcordance.md) is returned. A site
counts as called, by a caller or in the truth set, when some sample has
a frequency above 0 there; a recorded 0 is a site assessed without the
allele and is not a call.

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
benchmarkCallers(list(caller1 = whe1), truth = whe1)
#> DataFrame with 1 row and 8 columns
#>        caller   n_total        TP        FP        FN precision    recall
#>   <character> <integer> <integer> <integer> <integer> <numeric> <numeric>
#> 1     caller1         3         3         0         0         1         1
#>          F1
#>   <numeric>
#> 1         1
```
