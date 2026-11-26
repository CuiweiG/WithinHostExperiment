# Build a site frequency spectrum from a WithinHostExperiment

Constructs a [`WithinHostSFS`](WithinHostSFS-class.md) object from the
allele frequencies in a
[`WithinHostExperiment`](WithinHostExperiment.md). One SFS is built per
sample. Only segregating sites enter it, those whose frequency lies
strictly between 0 and 1, as in
[`calcNeutralityTests`](calcNeutralityTests.md): a frequency of 1, a
fixed difference from the reference, would otherwise fill the top bin of
an unfolded spectrum.

## Usage

``` r
buildSFS(
  whe,
  sampleIdx = 1L,
  nBins = 20L,
  fold = TRUE,
  usePassedOnly = TRUE,
  threshold = 0.03,
  genomeLength = NULL
)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md).

- sampleIdx:

  Integer or character. The one sample (column) to use, by position or
  by a column name that occurs once. Default: 1.

- nBins:

  Integer. Number of frequency bins, a whole number of at least 1
  (default: 20).

- fold:

  Logical. If `TRUE`, fold the SFS to minor allele frequency (0 to 0.5),
  giving `polarization = "none"`. If `FALSE`, bin the alternative allele
  frequency from 0 to 1, treating the reference allele as ancestral
  (`polarization = "reference"`). Default: `TRUE`.

- usePassedOnly:

  Logical. If `TRUE`, only use QC-passed variants. Default: `TRUE`. The
  spectrum's `maxFreq` slot then records the smallest `maxFreq` among
  the filters applied with [`flagISNV`](flagISNV.md), above which no
  call can appear; with the default [`ISNVFilter`](ISNVFilter.md) that
  is 0.5, so an unfolded spectrum is empty above 0.5. To keep the
  high-frequency classes, flag the unfiltered object with `maxFreq = 1`:
  [`flagISNV()`](flagISNV.md) adds to the flags already set, so it
  cannot restore calls that an earlier filter with a lower `maxFreq`
  removed.

- threshold:

  Numeric. Detection threshold (default: 0.03). The first bin starts at
  this value, not 0. It must be at least 0 and below the top of the
  spectrum: 0.5 when folded, 1 otherwise.

- genomeLength:

  Integer (optional). The genome length that the per-site statistics of
  [`neutralityFromSFS`](neutralityFromSFS.md) divide by. By default it
  is read from `seqlengths(rowRanges(whe))`, summed over every recorded
  sequence as in [`calcDiversity`](calcDiversity.md), and is `NA` when
  that is not recorded. The spectrum itself does not depend on it.

## Value

A [`WithinHostSFS`](WithinHostSFS-class.md) object.

## Examples

``` r
vcf <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(vcf,
    colData = S4Vectors::DataFrame(sample_id = "donor"),
    caller = "ivar")
sfs <- buildSFS(whe, genomeLength = 1000L)
sfs
#> WithinHostSFS with 15 variants in 20 bins
#>   sample: donor 
#>   folded: TRUE 
#>   polarization: none 
#>   threshold: 0.03 
#>   corrected: FALSE 
#>   genome: 1000 bp
#>   median depth: 1600 x
```
