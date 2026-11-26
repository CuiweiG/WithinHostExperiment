# Multiple testing correction for p-values

Thin wrapper around [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html)
for convenient multiple testing correction of p-value vectors.

## Usage

``` r
multitestCorrection(pvalues, method = "BH")
```

## Arguments

- pvalues:

  Numeric vector of raw p-values.

- method:

  Character scalar. Correction method passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) (default: `"BH"`
  for Benjamini-Hochberg).

## Value

Numeric vector of adjusted p-values (same length as `pvalues`).

## Examples

``` r
raw_p <- c(0.001, 0.04, 0.03, 0.5)
multitestCorrection(raw_p)
#> [1] 0.00400000 0.05333333 0.05333333 0.50000000
multitestCorrection(raw_p, method = "bonferroni")
#> [1] 0.004 0.160 0.120 1.000
```
