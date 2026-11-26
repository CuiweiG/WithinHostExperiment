# Exact bottleneck size estimation

Estimates transmission bottleneck size using one of three methods:

- `exact_bb`:

  Beta-binomial likelihood over the number of transmitted variant
  virions, which models finite recipient read depth, in the spirit of
  Sobel Leonard et al. (2017). Requires `altCount` and `totalDepth`
  assays. Two deviations from the published method are deliberate: when
  all or none of the transmitted virions carry the variant, the
  recipient read counts are modelled with a fixed sequencing error
  probability of `threshold / 2`, and the recipient detection threshold
  is not integrated out. The likelihood needs the recipient's read depth
  at each donor variant, so it cannot use a site where that depth is
  unknown: such a site would add the same constant at every `Nb`. This
  includes every site that the recipient's VCF file does not list, so
  with objects read from VCF files a donor variant the recipient lacks
  usually enters `exact_bb` not at all, and loss in transmission is then
  scored only by the other estimators. Such sites are left out, with a
  message, and counted in `n_excluded`, not in `n_variants`; a bootstrap
  interval (`nboot > 0`) resamples only the sites used. Where the depth
  is known but the recipient's alternative count is not, as in an iVar
  table without `ALT_DP`, the count is `round(altFreq * totalDepth)`;
  where the recipient has no frequency there either, it has no call, and
  the count is 0.

- `wright_fisher`:

  Approximate score in which one generation of binomial drift is
  combined with a beta kernel of concentration `max(Nb, 2)` around the
  drifted frequency. It mixes a density for detected variants with
  probabilities for lost variants, so it is a heuristic rather than a
  likelihood; prefer `exact_bb` or `presence_absence` for inference.
  Like `presence_absence`, it uses the donor variants with frequencies
  from `threshold` to `1 - threshold`.

- `presence_absence`:

  Presence-absence likelihood using \\P(\mathrm{detect} \| Nb, \nu) =
  1 - (1-\nu)^{Nb}\\. As in the published method, it uses the donor
  variants with frequencies from `threshold` to `1 - threshold`, and
  leaves out, with a message, those above: a donor variant at frequency
  1 that the recipient lacks would have probability 0 under every `Nb`.

## Usage

``` r
exactBottleneck(
  whe,
  pairId,
  maxNb = 500L,
  method = c("exact_bb", "wright_fisher", "presence_absence"),
  nboot = 0L,
  threshold = 0.03
)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md) with transmission
  pairs defined in `colData`.

- pairId:

  Character. Which pair to estimate.

- maxNb:

  Integer. Maximum bottleneck size to evaluate (default: 500); a single
  whole number of at least 1.

- method:

  Character. One of `"exact_bb"`, `"wright_fisher"`, or
  `"presence_absence"`.

- nboot:

  Integer. Number of bootstrap replicates for confidence intervals
  (default: 0, uses a likelihood-ratio interval). Bootstrapping
  resamples variant sites and uses the random number generator, so set a
  seed for reproducible intervals.

- threshold:

  Numeric. Variant calling threshold (default: 0.03), a single number
  above 0 and below 0.5: at 0 every recipient frequency counts as
  present, and from 0.5 no donor variant lies between `threshold` and
  `1 - threshold`.

## Value

A named list with components:

- `Nb`:

  Numeric. Maximum likelihood estimate of bottleneck size over the grid
  `1:maxNb`. A warning is issued when the estimate reaches `maxNb`,
  because the value is then censored by the grid.

- `ci`:

  Numeric vector of length 2. 95% confidence interval (likelihood-ratio
  or bootstrap percentile), also bounded by `maxNb`.

- `loglik`:

  Numeric. Log-likelihood at the MLE; for `"wright_fisher"`, the
  heuristic score at its maximiser.

- `method`:

  Character. Method used. For `"wright_fisher"` the `Nb` and `ci` are
  the maximiser and interval of a heuristic score, and a message says so
  when the method is run.

- `n_variants`:

  Integer. Number of donor variants used.

- `n_excluded`:

  Named integer. Donor variants at or above `threshold` that were left
  out, by reason: `recipient_failed_qc`, sites where QC removed a
  recipient call at or above `threshold`, which is not an absence (see
  [`asViralBottleneckInput`](asViralBottleneckInput.md)), left out by
  every method; `unknown_depth`, sites of unknown recipient depth, which
  only `exact_bb` leaves out; and `donor_near_fixation`, donor variants
  above `1 - threshold`, which `presence_absence` and `wright_fisher`
  leave out, as [`quickBottleneck`](quickBottleneck.md) does.

- `pair_id`:

  Character. Pair identifier.

## Details

Every method stops with an error when its log-likelihood is not finite
at any `Nb` of the grid, where it would otherwise return `Nb = 1`.

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
exactBottleneck(whe, pairId = "p1", maxNb = 50L,
    method = "presence_absence")
#> $Nb
#> [1] 2
#> 
#> $ci
#> [1] 1 4
#> 
#> $loglik
#> [1] -5.315958
#> 
#> $method
#> [1] "presence_absence"
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
