# Create a WithinHostExperiment from components

Constructs a `WithinHostExperiment` from assay matrices, genomic ranges,
and sample metadata.

## Usage

``` r
WithinHostExperiment(assays, rowRanges, colData, metadata = list())
```

## Arguments

- assays:

  Named list of matrices. Must include at least `altFreq`. If `qcPass`
  is not provided, it is automatically created as an all-TRUE matrix.

- rowRanges:

  `GRanges`. Genomic coordinates for each variant site.

- colData:

  `DataFrame` or `data.frame`. Sample metadata. Must contain a
  `sample_id` column.

- metadata:

  Named list (optional). Additional metadata such as `reference_genome`,
  `genome_length`, `calling_info`.

## Value

A `WithinHostExperiment` object.

## Examples

``` r
library(GenomicRanges)
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following object is masked from ‘package:utils’:
#> 
#>     data
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, scale, sequence, table,
#>     tapply, transform, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
gr <- GRanges("seg1", IRanges::IRanges(c(100, 200, 300), width = 1))
whe <- WithinHostExperiment(
    assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
    rowRanges = gr,
    colData = S4Vectors::DataFrame(sample_id = "S1"))
whe
#> WithinHostExperiment with 3 variant sites and 1 samples
#>   assays(2): altFreq, qcPass
#>   qcLog: 0 filtering step(s)
#>   transmission pairs: 0 
```
