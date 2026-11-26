# ConcordanceRule: agreement required between technical replicates

An S4 class saying how the technical replicates of a sample must agree
for their calls to pass. It is held in the `replicateConc` slot of an
[`ISNVFilter`](ISNVFilter-class.md), and [`flagISNV`](flagISNV.md)
applies it after the per-sample thresholds.

## Usage

``` r
# S4 method for class 'ConcordanceRule'
show(object)
```

## Arguments

- object:

  A `ConcordanceRule` object.

## Value

An object of class `ConcordanceRule`. Objects are created with the
constructor [`ConcordanceRule`](ConcordanceRule.md); the `show` method
prints the two settings and returns `NULL` invisibly.

## Details

Replicates are the samples that share a value of the `replicate_group`
column of `colData`. Every pair of replicates in a group is compared,
and a site that fails for any pair is flagged in every sample of the
group. A call is a non-missing `altFreq`, taken before any threshold is
applied, so the rule flags what
[`flagReplicateDiscordance`](flagReplicateDiscordance.md) flags with the
same two settings. Replicate intersection, which keeps only the calls
made in both replicates, is `ConcordanceRule(freqTolerance = Inf)`.
Frequencies lie in \[0, 1\], so no two differ by more than 1: a rule
with `freqTolerance` of 1 or more, `Inf` included, and
`requireBothDetected = FALSE` would flag nothing and is refused.

## Slots

- `freqTolerance`:

  Numeric scalar. The largest difference in `altFreq` allowed between
  two replicates that both made the call (default: 0.05); a difference
  equal to it, to within rounding error, passes. `Inf` compares no
  frequencies.

- `requireBothDetected`:

  Logical scalar. If `TRUE` (default), a site called in one replicate of
  a pair but not in the other fails.

## See also

[`ConcordanceRule`](ConcordanceRule.md) for the constructor,
[`ISNVFilter`](ISNVFilter.md) and [`flagISNV`](flagISNV.md) for applying
it, [`flagReplicateDiscordance`](flagReplicateDiscordance.md) for the
same check on its own.

## Examples

``` r
# Frequencies within 5 points, and called in both replicates
ConcordanceRule()
#> ConcordanceRule object
#>   freqTolerance:       0.05 
#>   requireBothDetected: TRUE 

# Replicate intersection only
ConcordanceRule(freqTolerance = Inf)
#> ConcordanceRule object
#>   freqTolerance:       Inf 
#>   requireBothDetected: TRUE 
```
