# Compute within-host diversity indices

Generic for computing diversity metrics from within-host variant data.
Methods should return a `DataFrame` with one row per sample.

Calculates diversity metrics per sample from a
[`WithinHostExperiment`](WithinHostExperiment.md). Implements standard
population genetics indices adapted for deep sequencing iSNV data,
following Lauring (2020).

## Usage

``` r
calcDiversity(x, ...)

# S4 method for class 'WithinHostExperiment'
calcDiversity(
  x,
  indices = c("shannon", "pi", "watterson", "richness", "simpson", "chao1"),
  usePassedOnly = TRUE,
  ...
)
```

## Arguments

- x:

  A [`WithinHostExperiment`](WithinHostExperiment.md).

- ...:

  Not used; supplying anything is an error.

- indices:

  Character vector. Which indices to compute: `"shannon"`, `"pi"`,
  `"watterson"`, `"richness"`, `"simpson"`, `"chao1"` (default: all
  six).

- usePassedOnly:

  Logical. If TRUE (default), only use variants where `qcPass == TRUE`.

## Value

A `DataFrame` of diversity metrics.

A `DataFrame` with one row per sample and columns `sample_id`,
`n_sites_assessed`, `detection_limit` and one column per requested
index.

## Details

**Genome length.** `"pi"` and `"watterson"` are per site and divide by
the genome length, which is read from `seqlengths(rowRanges(x))` as the
summed length of every sequence recorded there, including segments of a
segmented genome on which no site lies. It is a property of the
reference the calls were made against, not a number to type:
[`readWithinHost`](readWithinHost.md) records it from the VCF header's
`##contig` lines or from a `DNAStringSet` passed as `genome`, and
[`annotateCodonChange`](annotateCodonChange.md) from its reference. When
it is not recorded the function stops and says how to record it, for
example `seqlengths(x) <- c(seg = 29903L)`.

For `"watterson"`, the sample's median read depth is passed as the
sample size \\n\\, which for the reasons given in
[`tajimaD`](tajimaD.md) makes \\\theta_W\\ small;
[`calcNeutralityTests`](calcNeutralityTests.md) caps \\n\\ instead, so
the two need not agree.

**Reading a zero.** A diversity of zero means something different when
the sample had no call at all, when every call failed QC, and when calls
were present but below the frequency QC allows. Two columns keep these
apart: `n_sites_assessed`, the number of sites at which the sample has a
frequency, before QC, and `detection_limit`, the lowest frequency a call
could have kept through the filters applied (the strictest `minFreq`
among them), `NA` when no filter was applied or `usePassedOnly = FALSE`.
`richness` counts the sites that entered the indices after QC: the
segregating sites, whose frequency lies strictly between 0 and 1. A
frequency of 0, or of 1, a fixed difference from the reference such as
iVar reports, is not segregating within the host, so it enters none of
the indices; the site count of `"watterson"` and the singletons and
doubletons of `"chao1"` come from the same sites. A site counts as
assessed only when the sample has a call there: a VCF lists only the
sites a caller reported, so a site sequenced without a call cannot be
told apart from one not sequenced at all, and neither is counted.

## References

Lauring AS (2020). Within-Host Viral Diversity: A Window into Viral
Evolution. *Annu Rev Virol* 7:63-81.
[doi:10.1146/annurev-virology-010320-061642](https://doi.org/10.1146/annurev-virology-010320-061642)

## Examples

``` r
vcf <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
## the genome length comes from the reference, recorded on import
ref <- Biostrings::readDNAStringSet(system.file("extdata",
    "test_ref.fa", package = "WithinHostExperiment"))
whe <- readWithinHost(vcf,
    colData = S4Vectors::DataFrame(sample_id = "donor"),
    caller = "ivar", genome = ref)
calcDiversity(whe)
#> DataFrame with 1 row and 9 columns
#>     sample_id n_sites_assessed detection_limit  richness   shannon         pi
#>   <character>        <integer>       <numeric> <integer> <numeric>  <numeric>
#> 1       donor               15              NA        15   7.00462 0.00454524
#>    watterson   simpson     chao1
#>    <numeric> <numeric> <numeric>
#> 1 0.00188569    4.5424        15
```
