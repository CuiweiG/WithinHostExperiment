# Create an ISNVFilter object

Constructs a quality control filter configuration for use with
[`flagISNV`](flagISNV.md).

## Usage

``` r
ISNVFilter(
  minDepth = 100L,
  minFreq = 0.03,
  maxFreq = 0.5,
  minAltReads = 10L,
  maxStrandBias = 10,
  replicateConc = NULL
)
```

## Arguments

- minDepth:

  Integer. Minimum read depth (default: 100).

- minFreq:

  Numeric. Minimum alt allele frequency (default: 0.03).

- maxFreq:

  Numeric. Maximum alt allele frequency (default: 0.50).

- minAltReads:

  Integer. Minimum alt allele reads (default: 10).

- maxStrandBias:

  Numeric. Max strand bias fold-diff (default: 10).

- replicateConc:

  A [`ConcordanceRule`](ConcordanceRule.md) saying how technical
  replicates must agree, which [`flagISNV`](flagISNV.md) applies after
  the other thresholds, or `NULL` (default) for none.

## Value

An `ISNVFilter` object.

## Examples

``` r
ISNVFilter()
#> ISNVFilter object
#>   minDepth:      100 
#>   freq range:   [0.03, 0.5]
#>   minAltReads:   10 
#>   maxStrandBias: 10 
#>   replicateConc: none
ISNVFilter(minDepth = 500L, minFreq = 0.05)
#> ISNVFilter object
#>   minDepth:      500 
#>   freq range:   [0.05, 0.5]
#>   minAltReads:   10 
#>   maxStrandBias: 10 
#>   replicateConc: none
ISNVFilter(replicateConc = ConcordanceRule(freqTolerance = 0.05))
#> ISNVFilter object
#>   minDepth:      100 
#>   freq range:   [0.03, 0.5]
#>   minAltReads:   10 
#>   maxStrandBias: 10 
#>   replicateConc: freqTolerance 0.05, requireBothDetected TRUE
```
