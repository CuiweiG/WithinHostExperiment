# Example WithinHostExperiment dataset

A small synthetic [`WithinHostExperiment`](WithinHostExperiment.md) with
3 samples (donor, recipient, independent) containing iSNV data at
positions 33 to 945 of a synthetic genome segment. The donor and
recipient form a transmission pair (pair_id = "pair_1").

## Usage

``` r
data(example_whe)
```

## Format

A [`WithinHostExperiment`](WithinHostExperiment.md) with 28 variant
sites and 3 samples. Assays: altFreq, totalDepth, refCount, altCount,
qcPass. The length of the segment, 1000 bases, is recorded in
[`seqlengths()`](https://rdrr.io/pkg/Seqinfo/man/seqinfo.html) from
`test_ref.fa`.

## Source

Synthetic. Built by `inst/scripts/save_example_data.R` from
`test_donor.vcf`, `test_recipient.vcf` and `test_independent.vcf` in
`inst/extdata`, read with [`readWithinHost`](readWithinHost.md).

## Examples

``` r
data(example_whe)
example_whe
#> WithinHostExperiment with 28 variant sites and 3 samples
#>   assays(5): altFreq, totalDepth, refCount, altCount, qcPass
#>   qcLog: 0 filtering step(s)
#>   transmission pairs: 1 
transmissionPairs(example_whe)
#> DataFrame with 1 row and 3 columns
#>         donor   recipient     pair_id
#>   <character> <character> <character>
#> 1       donor   recipient      pair_1
```
