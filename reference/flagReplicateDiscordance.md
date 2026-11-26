# Flag iSNVs with replicate discordance

Identifies variant sites where technical replicates disagree and flags
them in the `qcPass` assay.

## Usage

``` r
flagReplicateDiscordance(whe, freqTolerance = 0.05, requireBothDetected = TRUE)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment-class.md) with
  `replicate_group` and `replicate_id` columns in `colData`.

- freqTolerance:

  Numeric scalar. Maximum allowable difference in `altFreq` between
  replicates (default: 0.05).

- requireBothDetected:

  Logical scalar. If `TRUE` (default), sites detected in only one
  replicate are also flagged. A `freqTolerance` of 1 or more with
  `requireBothDetected = FALSE` would flag nothing, since frequencies
  lie in \[0, 1\], and is an error, as it is for a
  [`ConcordanceRule`](ConcordanceRule-class.md).

## Value

A [`WithinHostExperiment`](WithinHostExperiment-class.md) with updated
`qcPass` and `qcLog`.

## Details

For each replicate group, the function compares `altFreq` between
replicates and flags sites where the difference exceeds `freqTolerance`;
a difference equal to it, to within rounding error, passes. Sites are
flagged in ALL samples of the replicate group.

Roder et al. (2023) reported that the choice of variant caller and the
use of replicate sequencing had the largest effects on single-nucleotide
variant discovery. An alternative approach is replicate intersection
(requiring detection in both replicates); use
`requireBothDetected = TRUE` with `freqTolerance = Inf` to achieve this.

## References

Roder AE et al. (2023). Optimized quantification of intra-host viral
diversity in SARS-CoV-2 and influenza virus sequence data. *mBio*
14:e01046-23.
[doi:10.1128/mbio.01046-23](https://doi.org/10.1128/mbio.01046-23)

## See also

[`flagISNV`](flagISNV.md) for standard QC filtering.

## Examples

``` r
rep1 <- system.file("extdata", "test_donor_rep1.vcf",
    package = "WithinHostExperiment")
rep2 <- system.file("extdata", "test_donor_rep2.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(c(rep1, rep2),
    colData = S4Vectors::DataFrame(
        sample_id = c("rep1", "rep2"),
        replicate_group = c("donor_reps", "donor_reps"),
        replicate_id = c("1", "2")),
    caller = "ivar")
whe <- flagReplicateDiscordance(whe, freqTolerance = 0.05)
qcLog(whe)
#> DataFrame with 2 rows and 6 columns
#>                    step           parameter       value n_flagged  n_passed
#>             <character>         <character> <character> <integer> <integer>
#> 1 replicate_discordance       freqTolerance        0.05         6        26
#> 2 replicate_discordance requireBothDetected        TRUE         0        26
#>             timestamp
#>           <character>
#> 1 2026-11-26T12:26:11
#> 2 2026-11-26T12:26:11
```
