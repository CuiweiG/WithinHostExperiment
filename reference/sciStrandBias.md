# Strand Coverage Imbalance (SCI) score

Computes a strand bias score that normalises the per-site
forward/reverse alternative-allele ratio by the global strand coverage
ratio, making it suitable for amplicon data where strand coverage is
inherently unequal. Mostefai et al. (2024) showed that strand bias
filters which ignore this imbalance can discard genuine variants; their
own metric is a per-strand binomial likelihood.

## Usage

``` r
sciStrandBias(fwd_alt, rev_alt, fwd_total, rev_total)
```

## Arguments

- fwd_alt:

  Numeric vector. Forward-strand alternative allele read counts. All
  four count vectors must be numeric and hold non-negative, finite
  values or `NA`; `NA` marks a site with no call (see Details).

- rev_alt:

  Numeric vector. Reverse-strand alternative allele read counts.

- fwd_total:

  Numeric vector. Forward-strand total read counts.

- rev_total:

  Numeric vector. Reverse-strand total read counts.

## Value

A numeric vector of SCI scores (log2 scale). Values near zero indicate
no strand bias; large absolute values indicate strand-specific
enrichment or depletion.

## Details

The SCI is defined as: \$\$\mathrm{SCI} = \log_2 \\\left(
\frac{f\_{\mathrm{alt}} / r\_{\mathrm{alt}}} {F\_{\mathrm{total}} /
R\_{\mathrm{total}}} \right)\$\$

where \\F\_{\mathrm{total}}\\ and \\R\_{\mathrm{total}}\\ are the sums
of forward and reverse total counts across all sites. A pseudocount of 1
is added to strand counts to avoid division by zero.

A site where any of the four counts is `NA`, as in a sample that made no
call there, scores `NA` and is left out of the sums, so the global ratio
is taken over the sites with all four counts.

## References

Mostefai F et al. (2024). Refining SARS-CoV-2 intra-host variation by
leveraging large-scale sequencing data. *NAR Genomics and
Bioinformatics* 6:lqae145.
[doi:10.1093/nargab/lqae145](https://doi.org/10.1093/nargab/lqae145)

## Examples

``` r
sci <- sciStrandBias(
    fwd_alt   = c(50, 2, 30),
    rev_alt   = c(45, 40, 28),
    fwd_total = c(500, 500, 500),
    rev_total = c(480, 480, 480)
)
sci
#> [1]  0.09008968 -3.83136321  0.03744161
```
