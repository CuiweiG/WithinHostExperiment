# Build SFS for all samples in a WithinHostExperiment

Convenience wrapper that calls [`buildSFS`](buildSFS.md) for each sample
and returns a named list of [`WithinHostSFS`](WithinHostSFS-class.md)
objects.

## Usage

``` r
buildSFSList(
  whe,
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

A named list of [`WithinHostSFS`](WithinHostSFS-class.md) objects.

## Examples

``` r
data(example_whe)
sfs_list <- buildSFSList(example_whe, genomeLength = 1000L)
length(sfs_list)
#> [1] 3
```
