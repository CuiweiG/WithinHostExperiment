# Chao1 richness estimator for iSNV data

Estimates total species (variant) richness using the Chao1 estimator
(Chao 1984): \$\$\hat{S}\_{Chao1} = S\_{obs} + \frac{f_1^2}{2 f_2}\$\$
where \\S\_{obs}\\ is the observed number of variants, \\f_1\\ is the
number of singletons (variants with `altCount == 1`), and \\f_2\\ is the
number of doubletons (variants with `altCount == 2`). When \\f_2 = 0\\,
the bias-corrected form is used: \\S\_{obs} + f_1 (f_1 - 1) / 2\\.

## Usage

``` r
chao1ISNV(counts, detected)
```

## Arguments

- counts:

  Integer vector. Alternative allele read counts (`altCount`) for each
  variant.

- detected:

  Integer scalar. Number of observed variants (\\S\_{obs}\\).

## Value

Numeric scalar. Chao1 richness estimate.

## References

Chao A (1984). Nonparametric estimation of the number of classes in a
population. *Scand J Statist* 11:265-270.

## Examples

``` r
chao1ISNV(c(1, 1, 2, 5, 10), detected = 5)
#> [1] 7
```
