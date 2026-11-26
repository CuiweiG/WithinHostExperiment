# Compare diversity between sample groups

Compare diversity between sample groups

## Usage

``` r
compareDiversity(whe, groupBy, index = "shannon")
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md).

- groupBy:

  Character. Column name in `colData` to group by.

- index:

  Character. Which index to compare (default: "shannon"). `"pi"` and
  `"watterson"` need the genome length recorded in
  `seqlengths(rowRanges(whe))`; see [`calcDiversity`](calcDiversity.md).

## Value

A `DataFrame` with test results.

## Examples

``` r
vcf <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(vcf,
    colData = S4Vectors::DataFrame(sample_id = "donor"),
    caller = "ivar")
## Single-group, returns diversity table:
compareDiversity(whe, groupBy = "sample_id", index = "shannon")
#> Need at least 2 groups for comparison; returning diversity table.
#> DataFrame with 1 row and 4 columns
#>     sample_id n_sites_assessed detection_limit   shannon
#>   <character>        <integer>       <numeric> <numeric>
#> 1       donor               15              NA   7.00462
```
