# ISNVFilter: Quality control filter for iSNV analysis

An S4 class holding parameters for filtering intra-host single
nucleotide variants. Used by [`flagISNV`](flagISNV.md).

## Usage

``` r
# S4 method for class 'ISNVFilter'
show(object)
```

## Arguments

- object:

  An `ISNVFilter` object.

## Value

An object of class `ISNVFilter`. Objects are created with the
constructor [`ISNVFilter`](ISNVFilter.md); the `show` method prints the
filter thresholds and returns `NULL` invisibly.

## Details

Default values are based on commonly used thresholds in the within-host
pathogen diversity literature:

- `minDepth = 100`:

  Minimum coverage for variant detection. This is a permissive floor:
  Grubaugh et al. (2019, Genome Biology) measured iSNVs above 3%
  accurately from at least 1,000 RNA copies sequenced to at least 400x.

- `minFreq = 0.03`:

  3% minor allele frequency threshold, the lower limit of accurate iSNV
  measurement in Grubaugh et al. (2019); Popa et al. (2020, Science
  Translational Medicine) used 2%.

- `maxFreq = 0.50`:

  Variants above 50% are considered consensus-level changes, not
  intra-host variants.

- `minAltReads = 10`:

  Minimum supporting reads to reduce false positives from sequencing
  errors.

- `maxStrandBias = 10`:

  Maximum fold-difference between forward and reverse strand alt allele
  counts. Excessive strand bias suggests systematic error.

## Slots

- `minDepth`:

  Integer scalar. Minimum read depth at iSNV site (default: 100).

- `minFreq`:

  Numeric scalar. Minimum alternative allele frequency (default: 0.03).

- `maxFreq`:

  Numeric scalar. Maximum alternative allele frequency; variants above
  this are excluded as near-fixed (default: 0.50).

- `minAltReads`:

  Integer scalar. Minimum reads supporting the alternative allele
  (default: 10).

- `maxStrandBias`:

  Numeric scalar. Maximum strand bias as fold-difference between strands
  (default: 10.0).

- `replicateConc`:

  A [`ConcordanceRule`](ConcordanceRule-class.md) saying how technical
  replicates must agree, or `NULL` (default) for no replicate check.
  [`flagISNV`](flagISNV.md) applies it after the other thresholds.

## References

Roder AE et al. (2023). Optimized quantification of intra-host viral
diversity in SARS-CoV-2 and influenza virus sequence data. *mBio*
14:e01046-23.
[doi:10.1128/mbio.01046-23](https://doi.org/10.1128/mbio.01046-23)

McCrone JT, Lauring AS (2016). Measurements of intrahost viral diversity
are extremely sensitive to systematic errors in variant calling. *J
Virol* 90:6884-6895.
[doi:10.1128/JVI.00667-16](https://doi.org/10.1128/JVI.00667-16)

Grubaugh ND et al. (2019). An amplicon-based sequencing framework for
accurately measuring intrahost virus diversity using PrimalSeq and iVar.
*Genome Biology* 20:8.
[doi:10.1186/s13059-018-1618-7](https://doi.org/10.1186/s13059-018-1618-7)

Mostefai F et al. (2024). Refining SARS-CoV-2 intra-host variation by
leveraging large-scale sequencing data. *NAR Genomics and
Bioinformatics* 6:lqae145.
[doi:10.1093/nargab/lqae145](https://doi.org/10.1093/nargab/lqae145)

## See also

[`ISNVFilter`](ISNVFilter.md) for the constructor,
[`flagISNV`](flagISNV.md) for applying filters.

## Examples

``` r
# Default filter
ISNVFilter()
#> ISNVFilter object
#>   minDepth:      100 
#>   freq range:   [0.03, 0.5]
#>   minAltReads:   10 
#>   maxStrandBias: 10 
#>   replicateConc: none

# Stringent filter
ISNVFilter(minDepth = 500L, minFreq = 0.05)
#> ISNVFilter object
#>   minDepth:      500 
#>   freq range:   [0.05, 0.5]
#>   minAltReads:   10 
#>   maxStrandBias: 10 
#>   replicateConc: none
```
