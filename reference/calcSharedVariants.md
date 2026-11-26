# Calculate shared variant statistics

Counts, for one transmission pair, the variant sites present in both
hosts, in the donor alone and in the recipient alone, using the calls
that passed QC.

## Usage

``` r
calcSharedVariants(whe, pairId, minDepth = NULL)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md).

- pairId:

  Character. Which pair to analyse.

- minDepth:

  Integer. The recipient depth at or above which an absent donor variant
  counts as absent at adequate depth. Defaults to the strictest
  `minDepth` of the filters applied with [`flagISNV`](flagISNV.md), or
  to 100 when none was applied.

## Value

A named list: `n_donor_only` (the four groups below added together),
`n_donor_only_adequate_depth`, `n_donor_only_low_depth`,
`n_donor_only_unknown_depth`, `n_donor_only_recipient_failed_qc`,
`n_recipient_only`, `n_shared`, `n_total`, `sharing_rate` (shared sites
over all sites) and `min_depth`, the threshold used.

## Details

A host carries a variant at a site when its frequency there is known and
above 0, as in [`exportPairFrequencies`](exportPairFrequencies.md); a
recorded frequency of 0 is a site that was assessed and found without
the variant.

A donor variant absent from the recipient means one thing when the
recipient was sequenced deeply enough at that site to have shown it, and
another when it was not. The donor-only count is therefore split by the
recipient's read depth at the site, from the `totalDepth` assay:
`n_donor_only_adequate_depth` where that depth is at least `minDepth`,
`n_donor_only_low_depth` where it is lower, and
`n_donor_only_unknown_depth` where it is missing. Depth is missing at
every site the recipient's caller did not report when the object was
read from VCF files, which list only the sites with a call, so with such
input most donor-only sites fall in the last group; a depth matrix from
a coverage file, or a table listing every assessed site with its depth
and a frequency of 0 where the recipient had no alternative reads
([`readWithinHostTable`](readWithinHostTable.md) with
`format = "generic"` and `aaf = 0`), is needed to tell the first two
apart. A fourth group, `n_donor_only_recipient_failed_qc`, holds the
sites where the recipient did make a call, at a frequency above 0, but
QC removed it, for example a call above `maxFreq` in a recipient where
the variant rose towards fixation: such a site is not absent, so it is
kept out of the three depth groups. Only a donor-only site at adequate
depth is evidence of loss in transmission.

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
calcSharedVariants(whe, pairId = "p1")
#> $n_donor_only
#> [1] 10
#> 
#> $n_donor_only_adequate_depth
#> [1] 0
#> 
#> $n_donor_only_low_depth
#> [1] 0
#> 
#> $n_donor_only_unknown_depth
#> [1] 10
#> 
#> $n_donor_only_recipient_failed_qc
#> [1] 0
#> 
#> $n_recipient_only
#> [1] 5
#> 
#> $n_shared
#> [1] 5
#> 
#> $n_total
#> [1] 20
#> 
#> $sharing_rate
#> [1] 0.25
#> 
#> $min_depth
#> [1] 100
#> 
```
