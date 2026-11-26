# Quick bottleneck size estimation (beta-binomial approximate MLE)

Estimates transmission bottleneck size using the approximate
beta-binomial method of Sobel Leonard et al. (2017), which does not
model sequencing depth explicitly. For the exact beta-binomial method
(which accounts for finite read depth) or other approaches
(presence-absence, KL divergence, Wright-Fisher), use the
ViralBottleneck package (Zheng et al. 2025, *Virus Evolution*
11:veaf071) via [`asViralBottleneckInput`](asViralBottleneckInput.md).

## Usage

``` r
quickBottleneck(whe, pairId, maxNb = 200L, threshold = 0.03)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md) with transmission
  pairs defined in `colData`.

- pairId:

  Character. Which pair to estimate.

- maxNb:

  Integer. Maximum bottleneck size to evaluate (default: 200); a single
  whole number of at least 1.

- threshold:

  Numeric. Variant calling threshold (default: 0.03), a single number
  above 0 and below 0.5.

## Value

A named list: Nb (MLE), ci (95% CI), loglik, n_variants, n_excluded,
pair_id. `n_excluded` is a named integer, as from
[`exactBottleneck`](exactBottleneck.md), of the donor variants at or
above `threshold` left out: `recipient_failed_qc`, sites where QC
removed a recipient call at or above `threshold`, which is not an
absence (see [`asViralBottleneckInput`](asViralBottleneckInput.md));
`unknown_depth`, always 0 here, since this method does not use depths;
and `donor_near_fixation`, donor variants above `1 - threshold`.

## Details

As in the published method, the donor variants used have frequencies
from `threshold` to `1 - threshold`; those above are left out, with a
message, since a donor variant at frequency 1 that the recipient lacks
has probability 0 under every `Nb`. The function stops with an error
when the log-likelihood is not finite at any `Nb` of the grid, where it
would otherwise return `Nb = 1`.

## References

Sobel Leonard A et al. (2017). Transmission Bottleneck Size Estimation
from Pathogen Deep-Sequencing Data, with an Application to Human
Influenza A Virus. *J Virol* 91:e00171-17.
[doi:10.1128/JVI.00171-17](https://doi.org/10.1128/JVI.00171-17)

## Examples

``` r
vcf1 <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
vcf2 <- system.file("extdata", "test_recipient.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(c(vcf1, vcf2),
    colData = S4Vectors::DataFrame(
        sample_id = c("donor", "recipient"),
        role = c("donor", "recipient"),
        pair_id = c("p1", "p1")),
    caller = "ivar")
quickBottleneck(whe, pairId = "p1", maxNb = 50L)
#> $Nb
#> [1] 3
#> 
#> $ci
#> [1] 1 7
#> 
#> $loglik
#> [1] -5.102543
#> 
#> $n_variants
#> [1] 15
#> 
#> $n_excluded
#> recipient_failed_qc       unknown_depth donor_near_fixation 
#>                   0                   0                   0 
#> 
#> $pair_id
#> [1] "p1"
#> 
```
