# Export data for ViralBottleneck package

Converts a [`WithinHostExperiment`](WithinHostExperiment.md) to the
input format expected by the ViralBottleneck package.

## Usage

``` r
asViralBottleneckInput(
  whe,
  pairId = NULL,
  threshold = 0.03,
  usePassedOnly = TRUE
)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md).

- pairId:

  Character (optional). Specific pair.

- threshold:

  Numeric. Variant calling threshold (default: 0.03).

- usePassedOnly:

  Logical (default: TRUE).

## Value

A `data.frame` with columns: pos, donor_freq, recipient_freq – a
simplified format for beta-binomial MLE. For the full ViralBottleneck
input format (per-base read counts), prepare data externally using the
read count assays (`altCount`, `refCount`).

## Details

One row is written for each donor variant at or above `threshold`. Where
the recipient has no call, its frequency is written as 0, which an
estimator reads as a variant lost in transmission. A recipient call at
or above `threshold` that QC removed is not such an absence: the call
may have been removed for lying above `maxFreq`, as when the variant
rose towards fixation in the recipient. With `usePassedOnly = TRUE`, a
donor variant at such a site is therefore left out, and a message gives
the number of sites. A removed call below `threshold` is an absence by
the estimators' own definition, and is written as 0.
(`n_donor_only_recipient_failed_qc` in
[`calcSharedVariants`](calcSharedVariants.md) counts every removed call
above 0, whatever its frequency.)
[`quickBottleneck`](quickBottleneck.md) and
[`exactBottleneck`](exactBottleneck.md) use the same donor variants.

For comprehensive bottleneck analysis with multiple methods
(presence-absence, KL divergence, binomial, beta-binomial
approximate/exact, and Wright-Fisher), use the ViralBottleneck package
(Zheng et al. 2025, *Virus Evolution* 11:veaf071). This function
prepares a simplified input suitable for the beta-binomial methods.

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
asViralBottleneckInput(whe, pairId = "p1")
#>    pos donor_freq recipient_freq
#> 1   45       0.35           0.28
#> 2  112       0.08           0.00
#> 3  198       0.22           0.00
#> 4  267       0.45           0.39
#> 5  334       0.12           0.00
#> 6  401       0.06           0.00
#> 7  456       0.18           0.00
#> 8  523       0.41           0.35
#> 9  589       0.09           0.00
#> 10 634       0.33           0.41
#> 11 701       0.15           0.00
#> 12 768       0.27           0.00
#> 13 812       0.04           0.00
#> 14 889       0.38           0.30
#> 15 945       0.11           0.00
```
