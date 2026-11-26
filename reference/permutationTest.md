# Permutation test for two-group comparison

Non-parametric permutation test comparing a summary statistic between
two groups. The observed test statistic is compared to a null
distribution obtained by randomly permuting group labels.

## Usage

``` r
permutationTest(values, groups, n_perm = 999L, statistic = "mean_diff")
```

## Arguments

- values:

  Numeric vector. The measured values (e.g.\\ a diversity index).

- groups:

  Factor or character vector (same length as `values`). Must contain
  exactly two unique levels.

- n_perm:

  Integer scalar. Number of permutations (default: 999). Permutation
  uses the random number generator, so set a seed for a reproducible
  p-value.

- statistic:

  Character scalar. The test statistic to compute: `"mean_diff"`
  (default) for difference in group means.

## Value

A list with components:

- observed:

  Numeric scalar. The observed test statistic.

- p_value:

  Numeric scalar. Two-sided permutation p-value.

- n_perm:

  Integer. Number of permutations performed.

## Examples

``` r
set.seed(42)
vals <- c(rnorm(10, 2), rnorm(10, 3))
grps <- rep(c("A", "B"), each = 10)
permutationTest(vals, grps, n_perm = 499)
#> $observed
#> [1] -0.2892465
#> 
#> $p_value
#> [1] 0.672
#> 
#> $n_perm
#> [1] 499
#> 
```
