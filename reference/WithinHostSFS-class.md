# WithinHostSFS: Site Frequency Spectrum for within-host variants

An S4 class representing the site frequency spectrum (SFS) of
within-host variants from deep sequencing data. The SFS bins variant
allele frequencies into discrete categories and stores metadata needed
for population-genetic inference.

## Usage

``` r
# S4 method for class 'WithinHostSFS'
show(object)
```

## Arguments

- object:

  A `WithinHostSFS` object.

## Value

An object of class `WithinHostSFS`, normally created by
[`buildSFS`](buildSFS.md); the `show` method prints a summary of the
spectrum and returns `NULL` invisibly.

## Slots

- `counts`:

  Integer vector. Number of iSNVs in each frequency bin. Length equals
  `length(breaks) - 1`.

- `breaks`:

  Numeric vector. Bin boundaries. Bins are right-closed, and the lowest
  break is included in the first bin. Default: 20 equal bins from the
  detection threshold to the maximum frequency.

- `folded`:

  Logical scalar. If `TRUE`, the SFS is folded (minor allele frequency,
  0 to 0.5); if `FALSE`, unfolded: the frequency of the alternative
  allele, 0 to 1.

- `polarization`:

  Character scalar: which allele is treated as ancestral. `"none"` for a
  folded spectrum, which carries no polarisation; `"reference"` for an
  unfolded spectrum built by [`buildSFS`](buildSFS.md), whose
  frequencies are those of the alternative allele, so the reference
  allele stands in for the ancestral one; `"outgroup"` for a spectrum
  whose frequencies are of alleles derived relative to an outgroup,
  which this package does not build but accepts. Statistics that need
  derived-allele frequencies, such as Fay and Wu's H, refuse `"none"`.

- `sampleId`:

  Character scalar. Sample identifier.

- `genomeLength`:

  Integer scalar. Reference genome length.

- `nSites`:

  Integer scalar. Number of variant sites counted in the spectrum, that
  is, with a frequency between the lowest and the highest break;
  [`buildSFS`](buildSFS.md) counts only segregating sites, whose
  frequency lies strictly between 0 and 1.

- `meanDepth`:

  Numeric scalar. Typical sequencing depth; [`buildSFS`](buildSFS.md)
  stores the median depth across sites.

- `threshold`:

  Numeric scalar. Detection threshold applied (e.g., 0.03). Sites below
  this are unobservable.

- `corrected`:

  Logical scalar. Whether ascertainment bias correction has been
  applied.

- `maxFreq`:

  Numeric scalar, above 0 and at most 1 (default 1). The highest
  alternative allele frequency a call could have kept through QC:
  [`buildSFS`](buildSFS.md) stores the smallest `maxFreq` among the
  filters applied with [`flagISNV`](flagISNV.md) when it uses QC-passed
  calls only, and 1 otherwise. Below 1, the spectrum holds no call above
  it, so the high-frequency classes are empty by construction;
  [`neutralityFromSFS`](neutralityFromSFS.md) refuses such an unfolded
  spectrum, and [`compareSFS`](compareSFS.md) compares only spectra with
  the same value.

## See also

[`buildSFS`](buildSFS.md) for constructing from a
`WithinHostExperiment`, [`neutralityFromSFS`](neutralityFromSFS.md) for
computing all neutrality tests from the SFS.

## Examples

``` r
data(example_whe)
## a folded spectrum of the first sample, in ten bins from 3% to 50%;
## the summary reports its polarisation, "none" for a folded spectrum,
## and the genome length, read from seqlengths(example_whe)
sfs <- buildSFS(example_whe, sampleIdx = 1L, nBins = 10L)
sfs
#> WithinHostSFS with 15 variants in 10 bins
#>   sample: donor 
#>   folded: TRUE 
#>   polarization: none 
#>   threshold: 0.03 
#>   corrected: FALSE 
#>   genome: 1000 bp
#>   median depth: 1600 x
```
