# Compare two within-host SFS objects

Computes summary statistics and a chi-squared test for homogeneity
between two [`WithinHostSFS`](WithinHostSFS-class.md) objects (e.g., two
hosts or two timepoints). A spectrum built before QC can be compared
with one built after it only when every filter applied used
`maxFreq = 1`; otherwise their caps differ and the function stops.

## Usage

``` r
compareSFS(sfs1, sfs2)
```

## Arguments

- sfs1:

  A [`WithinHostSFS`](WithinHostSFS-class.md) object.

- sfs2:

  A [`WithinHostSFS`](WithinHostSFS-class.md) object with the same
  breaks, folding, polarisation and `maxFreq` as `sfs1`.

## Value

A named list with components:

- chisq_stat:

  Chi-squared statistic, `NA` when no test is run.

- chisq_p:

  Chi-squared p-value, `NA` when no test is run.

- df:

  Degrees of freedom, `NA` when no test is run.

- min_expected:

  Smallest expected count under homogeneity, or `NA` when either
  spectrum is empty or fewer than two bins remain, and no test is run.

- chisq_approx_ok:

  `TRUE` when every expected count is at least 5, the condition under
  which [`chisq.test`](https://rdrr.io/r/stats/chisq.test.html) does not
  warn; `NA` alongside the other `NA` diagnostics when no test is run.

- n1:

  Total variants in sfs1.

- n2:

  Total variants in sfs2.

- proportions1:

  Proportion in each bin, sfs1.

- proportions2:

  Proportion in each bin, sfs2.

## Details

The test compares raw variant counts in matching bins, so the two
spectra must share their breaks, their folding, their polarisation and
their QC frequency cap (the `maxFreq` slot), and neither may be
bias-corrected: the counts of a corrected spectrum are rescaled and no
longer multinomial. The function stops, naming the difference, when any
of these does not hold. Bins where both spectra are empty are dropped,
and when only two bins remain
[`chisq.test`](https://rdrr.io/r/stats/chisq.test.html) applies Yates'
continuity correction. No test is run, and the statistic, p-value and
diagnostics are all `NA`, when either spectrum is empty or fewer than
two bins hold variants. The chi-squared approximation is unreliable when
an expected count falls below 5; instead of a warning, that condition is
reported in `chisq_approx_ok`, and the p-value should then be read with
caution.

## Examples

``` r
data(example_whe)
whe1 <- example_whe
sfs1 <- buildSFS(whe1, genomeLength = 1000L)
sfs2 <- buildSFS(whe1, sampleIdx = 2L, genomeLength = 1000L)
compareSFS(sfs1, sfs2)
#> $chisq_stat
#> [1] 7.638889
#> 
#> $chisq_p
#> [1] 0.9072783
#> 
#> $df
#> [1] 14
#> 
#> $min_expected
#> [1] 0.4
#> 
#> $chisq_approx_ok
#> [1] FALSE
#> 
#> $n1
#> [1] 15
#> 
#> $n2
#> [1] 10
#> 
#> $proportions1
#>  [1] 0.06666667 0.06666667 0.13333333 0.13333333 0.00000000 0.06666667
#>  [7] 0.06666667 0.00000000 0.06666667 0.00000000 0.06666667 0.00000000
#> [13] 0.06666667 0.06666667 0.06666667 0.00000000 0.06666667 0.06666667
#> [19] 0.00000000 0.00000000
#> 
#> $proportions2
#>  [1] 0.1 0.1 0.1 0.0 0.0 0.1 0.0 0.0 0.1 0.0 0.1 0.1 0.0 0.1 0.0 0.1 0.1 0.0 0.0
#> [20] 0.0
#> 
```
