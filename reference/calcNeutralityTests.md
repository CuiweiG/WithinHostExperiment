# Compute neutrality test statistics

Generic for computing neutrality tests from within-host variant data.
Methods should return a `DataFrame` with one row per sample.

Calculates neutrality test statistics per sample from a
[`WithinHostExperiment`](WithinHostExperiment.md). Implements Tajima's D
(Tajima 1989) and Fu's Fs (Fu 1997), adapted for deep sequencing iSNV
data.

## Usage

``` r
calcNeutralityTests(x, ...)

# S4 method for class 'WithinHostExperiment'
calcNeutralityTests(
  x,
  indices = c("tajimaD", "fusFs"),
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

  Character vector. Which tests to compute: `"tajimaD"`, `"fusFs"`
  (default: both).

- usePassedOnly:

  Logical. If TRUE (default), only use variants where `qcPass == TRUE`.

## Value

A `DataFrame` of neutrality statistics.

A `DataFrame` with one row per sample and columns `sample_id`, `n_sites`
(the number of segregating sites each statistic rests on, that is, the
sites at which the sample's frequency after QC lies strictly between 0
and 1; a frequency of 0, or of 1, a fixed difference from the reference,
is not segregating within the host) and one column per requested test.
With \\n\\ the effective sample size, the sample's median depth capped
at 100, a value is `NA` in four cases (see [`tajimaD`](tajimaD.md) and
[`fusFs`](fusFs.md)): Tajima's D when `n_sites` is below three, or when
\\n\\ is below four; Fu's Fs when `n_sites` is 0, or when `n_sites + 1`
exceeds \\n\\, which is the case for every sample with 100 or more
segregating sites. For each case that occurs, one warning names how many
samples it affected.

## Details

No genome length is needed. Both statistics depend on the nucleotide
diversity summed over the genome, \\\hat{\pi} L\\, which is computed
directly from the frequencies, so a per-site value would only be
multiplied back by the length it was divided by.

Fu's Fs sets the number of haplotypes to `n_sites + 1`, the
infinite-sites maximum, as [`fusFs`](fusFs.md) explains. It therefore
carries no haplotype information, depends only on `n_sites`, the
effective sample size and \\\hat{\pi} L\\, and is pushed below 0 under
neutrality, increasingly so as \\\theta\\ grows (about -7 at \\\theta =
5\\ with \\n = 100\\, near 0 at \\\theta\\ of 1 or less); a negative
value is not evidence of expansion.

## References

Tajima F (1989). Statistical method for testing the neutral mutation
hypothesis by DNA polymorphism. *Genetics* 123:585-595.

Fu YX (1997). Statistical tests of neutrality of mutations against
population growth, hitchhiking and background selection. *Genetics*
147:915-925.

## Examples

``` r
vcf <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(vcf,
    colData = S4Vectors::DataFrame(sample_id = "donor"),
    caller = "ivar")
calcNeutralityTests(whe)
#> DataFrame with 1 row and 4 columns
#>     sample_id   n_sites   tajimaD     fusFs
#>   <character> <integer> <numeric> <numeric>
#> 1       donor        15   1.61544 -0.382119
```
