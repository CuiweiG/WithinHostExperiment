# Create a ConcordanceRule object

Says how the technical replicates of a sample must agree, for use as the
`replicateConc` of an [`ISNVFilter`](ISNVFilter.md).

## Usage

``` r
ConcordanceRule(freqTolerance = 0.05, requireBothDetected = TRUE)
```

## Arguments

- freqTolerance:

  Numeric. The largest difference in allele frequency allowed between
  two replicates that both made the call (default: 0.05); `Inf` compares
  no frequencies.

- requireBothDetected:

  Logical. Whether a site called in one replicate of a pair but not in
  the other fails (default: `TRUE`).

## Value

A [`ConcordanceRule`](ConcordanceRule-class.md) object.

## Details

Neither argument is converted from another type: a factor, a string or a
logical `freqTolerance` is an error, and so is a `requireBothDetected`
that is not `TRUE` or `FALSE`, as in
[`flagReplicateDiscordance`](flagReplicateDiscordance.md). An integer
tolerance is stored as a double.

## See also

[`flagISNV`](flagISNV.md), which applies the rule, and
[`flagReplicateDiscordance`](flagReplicateDiscordance.md).

## Examples

``` r
rule <- ConcordanceRule(freqTolerance = 0.05)
rule
#> ConcordanceRule object
#>   freqTolerance:       0.05 
#>   requireBothDetected: TRUE 
ISNVFilter(replicateConc = rule)
#> ISNVFilter object
#>   minDepth:      100 
#>   freq range:   [0.03, 0.5]
#>   minAltReads:   10 
#>   maxStrandBias: 10 
#>   replicateConc: freqTolerance 0.05, requireBothDetected TRUE
```
