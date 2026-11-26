# Flag iSNV sites based on quality criteria

Generic for applying quality control filters to within-host variant
data. The default method operates on
[`WithinHostExperiment`](WithinHostExperiment.md) objects. Other
packages can define methods for alternative containers.

Applies quality control filters to a
[`WithinHostExperiment`](WithinHostExperiment.md). Instead of removing
variants, this function sets the corresponding entries in the `qcPass`
assay to `FALSE`. Each filtering step is recorded in the `qcLog` for
full reproducibility.

## Usage

``` r
flagISNV(x, filter = ISNVFilter(), ...)

# S4 method for class 'WithinHostExperiment'
flagISNV(x, filter = ISNVFilter(), ...)
```

## Arguments

- x:

  A [`WithinHostExperiment`](WithinHostExperiment.md) object.

- filter:

  An [`ISNVFilter`](ISNVFilter.md) object (default:
  [`ISNVFilter()`](ISNVFilter.md)).

- ...:

  Must be empty. Any further argument is an error, so that a misspelled
  argument name cannot be silently ignored.

## Value

The filtered object (same class as input).

A [`WithinHostExperiment`](WithinHostExperiment.md) with updated
`qcPass` assay and `qcLog`.

## Details

Filtering steps are applied in order:

1.  Minimum read depth (`minDepth`)

2.  Minimum alternative allele frequency (`minFreq`)

3.  Maximum alternative allele frequency (`maxFreq`)

4.  Minimum alternative allele read count (`minAltReads`)

5.  Strand bias (`maxStrandBias`) – only if `fwdAltCount` and
    `revAltCount` assays exist

6.  Replicate concordance (`replicateConc`) – only if the filter holds a
    [`ConcordanceRule`](ConcordanceRule.md)

The filter is checked with
[`validObject()`](https://rdrr.io/r/methods/validObject.html) first, so
a threshold or a rule changed with `@<-` into one its class refuses is
an error rather than a step applied.

Base and mapping quality are not filtered here: the object holds no
per-site quality, so those thresholds belong in the variant caller.

The replicate step compares every pair of samples that share a value of
the `replicate_group` column of `colData` and flags a site that fails
the rule in every sample of the group. It compares the calls as made,
before the thresholds above, so it flags what
[`flagReplicateDiscordance`](flagReplicateDiscordance.md) flags with the
same settings, and it records the same two `qcLog` rows. A rule with no
`replicate_group` column, or with no group of two or more samples, is an
error rather than a step that checks nothing. Samples outside any such
group are named in a warning and are filtered by the other thresholds
only.

## See also

[`ISNVFilter`](ISNVFilter.md) for filter configuration,
[`ConcordanceRule`](ConcordanceRule.md) for the replicate rule,
[`passedISNV`](WHE-accessors.md) for extracting passed variants,
[`qcSummary`](qcSummary.md) for QC report,
[`flagReplicateDiscordance`](flagReplicateDiscordance.md) for replicate
QC.

## Examples

``` r
vcf <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(vcf,
    colData = S4Vectors::DataFrame(sample_id = "donor"),
    caller = "ivar")
whe <- flagISNV(whe, ISNVFilter(minDepth = 500L))
qcLog(whe)
#> DataFrame with 4 rows and 6 columns
#>               step   parameter       value n_flagged  n_passed
#>        <character> <character> <character> <integer> <integer>
#> 1     depth_filter    minDepth         500         1        14
#> 2  min_freq_filter     minFreq        0.03         0        14
#> 3  max_freq_filter     maxFreq         0.5         0        14
#> 4 alt_reads_filter minAltReads          10         0        14
#>             timestamp
#>           <character>
#> 1 2026-11-26T08:58:34
#> 2 2026-11-26T08:58:34
#> 3 2026-11-26T08:58:34
#> 4 2026-11-26T08:58:34
passedISNV(whe)
#> WithinHostExperiment with 14 variant sites and 1 samples
#>   assays(5): altFreq, totalDepth, refCount, altCount, qcPass
#>   qcLog: 4 filtering step(s)
#>   transmission pairs: 0 
```
