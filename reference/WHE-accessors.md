# Access the QC audit log

Returns the quality control log recording all filtering steps applied to
a [`WithinHostExperiment`](WithinHostExperiment.md).

Returns a subset of the
[`WithinHostExperiment`](WithinHostExperiment.md) containing only
variant sites where at least one sample passed QC.

Identifies donor-recipient pairs based on `role` and `pair_id` columns
in `colData`.

## Usage

``` r
qcLog(x)

qcFilters(x)

passedISNV(x)

transmissionPairs(x)

qcLog(x) <- value

qcFilters(x) <- value

# S4 method for class 'WithinHostExperiment'
qcLog(x)

# S4 method for class 'WithinHostExperiment'
qcFilters(x)

# S4 method for class 'WithinHostExperiment'
qcLog(x) <- value

# S4 method for class 'WithinHostExperiment'
qcFilters(x) <- value

# S4 method for class 'WithinHostExperiment'
passedISNV(x)

# S4 method for class 'WithinHostExperiment'
transmissionPairs(x)
```

## Arguments

- x:

  A [`WithinHostExperiment`](WithinHostExperiment.md) object.

- value:

  A list replacement value.

## Value

A `DataFrame` with columns: step, parameter, value, n_flagged, n_passed,
timestamp.

A list of [`ISNVFilter`](ISNVFilter.md) objects.

A [`WithinHostExperiment`](WithinHostExperiment.md) (subset of rows).

A `DataFrame` with columns: donor (sample_id), recipient (sample_id),
pair_id.

The modified [`WithinHostExperiment`](WithinHostExperiment.md).

The modified [`WithinHostExperiment`](WithinHostExperiment.md).

## Examples

``` r
gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(1:3, width = 1))
whe <- WithinHostExperiment(
    assays = list(altFreq = matrix(runif(3), ncol = 1)),
    rowRanges = gr,
    colData = S4Vectors::DataFrame(sample_id = "S1"))
qcLog(whe)
#> DataFrame with 0 rows and 6 columns
gr <- GenomicRanges::GRanges("seg", IRanges::IRanges(1:3, width = 1))
whe <- WithinHostExperiment(
    assays = list(altFreq = matrix(runif(3), ncol = 1)),
    rowRanges = gr,
    colData = S4Vectors::DataFrame(sample_id = "S1"))
passedISNV(whe)
#> WithinHostExperiment with 3 variant sites and 1 samples
#>   assays(2): altFreq, qcPass
#>   qcLog: 0 filtering step(s)
#>   transmission pairs: 0 
```
