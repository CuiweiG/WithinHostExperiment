# Simpson diversity index for iSNV frequencies

Calculates Simpson diversity summed over biallelic iSNV sites: \$\$D =
\sum_i \left\[1 - p_i^2 - (1 - p_i)^2\right\] = \sum_i 2 p_i (1 -
p_i)\$\$ where \\p_i\\ is the alternative allele frequency at site
\\i\\. Each site contributes the probability that two alleles drawn at
random at that site differ (at most 0.5), so the index is the expected
number of differing sites. Higher values indicate greater diversity.

## Usage

``` r
simpsonISNV(freq)
```

## Arguments

- freq:

  Numeric vector. Alternative allele frequencies between 0 and 1;
  anything else is an error. Values of 0 and 1, which carry no
  diversity, and NA are excluded.

## Value

Numeric scalar. Simpson diversity index.

## Examples

``` r
simpsonISNV(c(0.05, 0.15, 0.30, 0.45))
#> [1] 1.265
simpsonISNV(0.5) # maximum per-site diversity = 0.5
#> [1] 0.5
```
