# Correct ascertainment bias in the within-host SFS

Deep sequencing imposes a detection threshold (e.g., 3%): sites with
true frequency below this threshold are unobservable. This inflates the
apparent fraction of moderate-frequency variants relative to rare
variants.

## Usage

``` r
correctSFSBias(sfs, method = c("truncation", "binomial"))
```

## Arguments

- sfs:

  A [`WithinHostSFS`](WithinHostSFS-class.md) object.

- method:

  Character. Correction method: `"truncation"` (default) rescales the
  whole spectrum for the unobservable part of the frequency range;
  `"binomial"` models the probability of detecting a variant at each
  frequency given the read depth.

## Value

A corrected [`WithinHostSFS`](WithinHostSFS-class.md) with
`corrected = TRUE`. Counts are no longer raw variant counts, and
`nSites` changes accordingly, which also changes the number of
segregating sites used by [`neutralityFromSFS`](neutralityFromSFS.md).
Use uncorrected spectra for [`compareSFS`](compareSFS.md).

## Details

This function applies a correction based on the expected number of
unobserved sites under a neutral beta distribution model: given the
observed SFS and depth-dependent detectability, it estimates the missing
mass in the lowest frequency bins.

The `"truncation"` method assumes a uniform density of variant
frequencies and multiplies every bin by \\f\_{max} / (f\_{max} -
\mathrm{threshold})\\, the reciprocal of the observable fraction of the
frequency range. It therefore rescales the spectrum without changing its
shape. The `"binomial"` method uses the stored depth to compute
\\P(\mathrm{detect} \| \nu, n) = P(\mathrm{Alt} \ge k\_{\min} \| n,
\nu)\\ where \\k\_{\min} = \lceil n \times \mathrm{threshold} \rceil\\,
and reweights each bin by the inverse of its detectability.
Detectability is floored at 0.01, so a bin is inflated at most one
hundredfold.

Either way the corrected counts stay whole numbers of sites, as
[`fusFs`](fusFs.md) needs a whole number of segregating sites. The
corrected total is rounded once, and its sites are shared out among the
bins by largest remainder, ties going to the lower bins; rounding each
bin on its own would leave every bin of a few sites unchanged. `nSites`
is that total.

## Examples

``` r
vcf <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(vcf,
    colData = S4Vectors::DataFrame(sample_id = "donor"),
    caller = "ivar")
sfs <- buildSFS(whe, genomeLength = 1000L)
sfs_corr <- correctSFSBias(sfs, method = "binomial")
sfs_corr
#> WithinHostSFS with 15 variants in 20 bins
#>   sample: donor 
#>   folded: TRUE 
#>   polarization: none 
#>   threshold: 0.03 
#>   corrected: TRUE 
#>   genome: 1000 bp
#>   median depth: 1600 x
```
