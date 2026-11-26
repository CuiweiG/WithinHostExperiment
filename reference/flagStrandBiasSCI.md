# Flag strand-biased iSNVs using the SCI metric

Applies Strand Coverage Imbalance filtering to a
[`WithinHostExperiment`](WithinHostExperiment.md). Sites whose absolute
SCI score exceeds `maxSCI` are flagged in the `qcPass` assay.

## Usage

``` r
flagStrandBiasSCI(whe, maxSCI = 1.5)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md) object. Must
  contain the assays `fwdAltCount`, `revAltCount`, `fwdTotalCount`, and
  `revTotalCount`.

- maxSCI:

  Numeric scalar. Maximum allowable absolute SCI score (default: 1.5).
  Sites with `|SCI| > maxSCI` are flagged.

## Value

A [`WithinHostExperiment`](WithinHostExperiment.md) with updated
`qcPass` and `qcLog`.

## Details

Unlike the simple fold-difference strand bias in
[`flagISNV`](flagISNV.md), this function uses the SCI score, which
normalises by the global strand coverage ratio (see
[`sciStrandBias`](sciStrandBias.md)). This is preferred for
amplicon-based sequencing where forward and reverse strand coverage is
inherently unbalanced.

The global strand ratio is computed for each sample over the sites that
sample called, that is, the sites where its four count assays are not
`NA`. Sites the sample did not call keep their `qcPass` value.

## References

Mostefai F et al. (2024). Refining SARS-CoV-2 intra-host variation by
leveraging large-scale sequencing data. *NAR Genomics and
Bioinformatics* 6:lqae145.
[doi:10.1093/nargab/lqae145](https://doi.org/10.1093/nargab/lqae145)

## See also

[`sciStrandBias`](sciStrandBias.md) for the underlying computation,
[`flagISNV`](flagISNV.md) for standard QC.

## Examples

``` r
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
library(GenomicRanges)
gr <- GRanges("seg", IRanges::IRanges(start = 1:3, width = 1))
assays <- list(
    altFreq      = matrix(c(0.1, 0.05, 0.2), ncol = 1),
    totalDepth   = matrix(c(1000L, 1000L, 1000L), ncol = 1),
    qcPass       = matrix(TRUE, nrow = 3, ncol = 1),
    fwdAltCount  = matrix(c(50L, 2L, 30L), ncol = 1),
    revAltCount  = matrix(c(45L, 40L, 28L), ncol = 1),
    fwdTotalCount = matrix(c(500L, 500L, 500L), ncol = 1),
    revTotalCount = matrix(c(480L, 480L, 480L), ncol = 1)
)
whe <- WithinHostExperiment(
    rowRanges = gr, assays = assays,
    colData = S4Vectors::DataFrame(sample_id = "s1"))
whe <- flagStrandBiasSCI(whe, maxSCI = 1.5)
qcLog(whe)
#> DataFrame with 1 row and 6 columns
#>                     step   parameter       value n_flagged  n_passed
#>              <character> <character> <character> <integer> <integer>
#> 1 strand_bias_sci_filter      maxSCI         1.5         1         2
#>             timestamp
#>           <character>
#> 1 2026-11-26T14:39:41
```
