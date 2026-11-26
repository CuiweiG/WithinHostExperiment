# Compute all neutrality statistics from a WithinHostSFS

Computes Tajima's D, Fu's Fs, and Fay and Wu's H from a single
[`WithinHostSFS`](WithinHostSFS-class.md) object. This is the
recommended interface for neutrality testing: build the SFS once, then
derive all statistics from it.

## Usage

``` r
neutralityFromSFS(sfs, n = NULL)
```

## Arguments

- sfs:

  A [`WithinHostSFS`](WithinHostSFS-class.md) object.

- n:

  Integer. Effective sample size for neutrality tests. For deep
  sequencing data, do **not** use raw depth. Default:
  `min(100, meanDepth)`, where the spectrum's `meanDepth` slot holds the
  median depth across sites, as recommended in the
  [`tajimaD`](tajimaD.md) documentation.

## Value

A named list with components:

- S:

  Number of segregating sites.

- n:

  Effective sample size used.

- pi:

  Nucleotide diversity per site, Fay and Wu's \\\hat{\pi}\\ on the
  classes \\i\\; this is the value passed to [`tajimaD`](tajimaD.md)
  and, times the genome length, to [`fusFs`](fusFs.md).

- theta_W:

  Watterson's theta per site.

- theta_H:

  Fay and Wu's theta_H per site.

- tajimaD:

  Tajima's D.

- fusFs:

  Fu's Fs, with the number of haplotypes set to \\S + 1\\, the
  infinite-sites maximum. As [`fusFs`](fusFs.md) explains, it then
  carries no haplotype information, depends only on \\S\\, \\n\\ and
  \\\hat{\pi}\\, and is pushed below 0 under neutrality, increasingly so
  as \\\theta\\ grows, so a negative value is not evidence of expansion.

- fayWuH:

  Fay and Wu's H.

- sfs_corrected:

  Whether the input SFS was bias-corrected.

- polarization:

  The input's polarisation, `"reference"` or `"outgroup"`: which allele
  `theta_H` and `fayWuH` treat as ancestral.

## Details

**Fay and Wu's H** measures an excess of high-frequency derived
variants, which can indicate positive selection or genetic hitchhiking
(Fay and Wu 2000). It needs the frequency of the derived allele at each
site, so the function stops on a folded spectrum
(`polarization = "none"`): folding maps a derived allele at 0.9 and one
at 0.1 to the same bin. For Tajima's D and Fu's Fs alone, which need no
polarisation, use [`calcNeutralityTests`](calcNeutralityTests.md). An
unfolded spectrum built by [`buildSFS`](buildSFS.md) has
`polarization = "reference"`: its frequencies are those of the
alternative allele, so H here assumes that the reference allele is
ancestral at every site, which holds only as far as the reference is
close to the population the host was infected from. H is defined as:
\$\$H = \hat{\pi} - \hat{\theta}\_H\$\$ where \\\hat{\pi} =
\frac{2}{n(n-1)} \sum_i i (n - i) \xi_i\\ and \\\hat{\theta}\_H =
\frac{2}{n(n-1)} \sum_i i^2 \xi_i\\ are Fay and Wu's estimators on the
same classes, and \\\xi_i\\ is the count of sites with \\i\\ derived
alleles. Each bin stands for the class \\i\\ nearest its midpoint,
`round(mid * n)`, held within 1 to \\n - 1\\. On matching classes both
estimators have the same expectation under neutrality, so H is centred
on 0 when every class from 1 to \\n - 1\\ is observed and each bin holds
one class; the detection threshold, which removes the lowest classes,
and wider bins move it (with the default threshold of 0.03 and \\n =
100\\, leaving out classes 1 and 2 shifts it by about \\-0.04\theta\\),
and \\\hat{\pi}\\ carries the \\n/(n-1)\\ correction that
[`calcNeutralityTests`](calcNeutralityTests.md) applies. A negative H
indicates an excess of high-frequency variants.

Because H weights the high-frequency classes most, it is uninterpretable
when QC has removed them. The function therefore stops on a spectrum
whose `maxFreq` slot is below 1, as it is after
[`flagISNV`](flagISNV.md) with the default
[`ISNVFilter`](ISNVFilter.md), which drops every call above 0.5: H could
not then reveal an excess of high-frequency variants, whatever the data.
Flag the unfiltered object, before any filter with a lower `maxFreq`,
with `ISNVFilter(maxFreq = 1)`, then `buildSFS(fold = FALSE)`.

The per-site values need the spectrum's genome length, which
[`buildSFS`](buildSFS.md) takes from the object; a spectrum without one
is an error rather than a set of per-genome values labelled per site.

## References

Fay JC, Wu CI (2000). Hitchhiking under positive Darwinian selection.
*Genetics* 155:1405-1413.

## Examples

``` r
vcf <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(vcf,
    colData = S4Vectors::DataFrame(sample_id = "donor"),
    caller = "ivar")
sfs <- buildSFS(whe, genomeLength = 1000L, fold = FALSE)
neutralityFromSFS(sfs)
#> $S
#> [1] 15
#> 
#> $n
#> [1] 100
#> 
#> $pi
#> [1] 0.004580202
#> 
#> $theta_W
#> [1] 0.00289722
#> 
#> $theta_H
#> [1] 0.001924848
#> 
#> $tajimaD
#> [1] 1.60772
#> 
#> $fusFs
#> [1] -0.3911511
#> 
#> $fayWuH
#> [1] 0.002655354
#> 
#> $sfs_corrected
#> [1] FALSE
#> 
#> $polarization
#> [1] "reference"
#> 
```
