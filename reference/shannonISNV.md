# Shannon entropy from iSNV frequencies

Calculates per-site biallelic Shannon entropy, summed across all iSNV
sites: \$\$H = -\sum_i \[p_i \log(p_i) + (1-p_i) \log(1-p_i)\]\$\$ where
\\p_i\\ is the alternative allele frequency at site \\i\\. This
formulation treats each iSNV as an independent biallelic locus,
following McCrone & Lauring (2016) and Popa et al. (2020). Higher values
indicate greater within-host diversity.

## Usage

``` r
shannonISNV(freq)
```

## Arguments

- freq:

  Numeric vector. Alternative allele frequencies between 0 and 1;
  anything else is an error. Values of 0 and 1, which carry no
  diversity, and NA are excluded.

## Value

Numeric scalar. Shannon entropy.

## References

Lauring AS (2020). Within-Host Viral Diversity: A Window into Viral
Evolution. *Annu Rev Virol* 7:63-81.
[doi:10.1146/annurev-virology-010320-061642](https://doi.org/10.1146/annurev-virology-010320-061642)

## Examples

``` r
shannonISNV(c(0.05, 0.15, 0.30, 0.45))
#> [1] 1.920227
shannonISNV(0.5) # log(2)
#> [1] 0.6931472
```
