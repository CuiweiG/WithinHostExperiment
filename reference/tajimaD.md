# Tajima's D statistic

Computes Tajima's D (Tajima 1989) from summary statistics of a
within-host pathogen population. The statistic measures the difference
between nucleotide diversity (\\\hat{\pi}\\) and Watterson's estimator
(\\\hat{\theta}\_W = S / a_1\\), normalised by its expected standard
deviation under neutrality: \$\$D = \frac{\hat{\pi} - S/a_1}{\sqrt{e_1
S + e_2 S(S-1)}}\$\$ A significantly negative D indicates an excess of
rare variants (e.g., population expansion or purifying selection), while
a positive D suggests balancing selection or population contraction.

## Usage

``` r
tajimaD(S, n, pi_hat, genomeLength)
```

## Arguments

- S:

  Integer scalar. Number of segregating sites.

- n:

  Numeric scalar. Sample size. For deep sequencing data, do **not** use
  raw read depth directly – the classical Tajima's D assumes independent
  lineages, whereas sequencing reads oversample a much smaller viral
  population. Use an effective sample size (e.g., `n = 50` or `n = 100`)
  that reflects the number of independently sampled viral genomes, or
  cap the depth at a biologically meaningful value. See Details.

- pi_hat:

  Numeric scalar. Nucleotide diversity estimate (per-site).

- genomeLength:

  Integer scalar. Reference genome length in bp.

## Value

Numeric scalar. Tajima's D, or `NA` with a warning when `S` is below
three or `n` rounds to fewer than four.

## Details

**Fewer than three segregating sites:** with one segregating site D is a
monotone function of that site's frequency, and with two it is a
function of two frequencies; neither describes a frequency spectrum, and
at `S = 0` the statistic is 0/0. The function therefore returns `NA`,
with a warning, whenever `S < 3`. Report `S` alongside every D, as
[`calcNeutralityTests`](calcNeutralityTests.md) does in its `n_sites`
column. With two or three sequences the variance of D is zero (\\c_1 =
c_2 = 0\\), so `n` must round to at least four; below that D is likewise
`NA`.

**Deep sequencing caveat:** The classical Tajima (1989) formula was
designed for independent haplotype samples (e.g., Sanger-sequenced
clones or individual genomes). In deep sequencing, reads are not
independent lineages; they are technical replicates of a viral
population whose effective size is much smaller than the read depth.
Using raw read depth (e.g., n = 2000) inflates the harmonic number \\a_1
= \sum 1/i\\, making \\\theta_W\\ artificially small and biasing D
towards positive values.

We recommend capping `n` at a biologically meaningful effective sample
size. Lauring (2020 *Annu Rev Virol*) suggests n = 100 as a conservative
proxy. The high-level [`calcNeutralityTests`](calcNeutralityTests.md)
method applies this cap automatically (`n = min(depth, 100)`).

**Interpretation warning:** Deep-sequencing D values are *not* directly
comparable to classical population-genetic D computed from independently
sampled haplotypes. We recommend reporting D as a *descriptive
statistic* of the site frequency spectrum (e.g., "the majority of
samples show negative D, consistent with an excess of rare variants"),
rather than as a formal neutrality test with p-value interpretation. For
rigorous neutrality testing on deep sequencing data, consider
simulation-based approaches that account for the read-sampling process.

## References

Tajima F (1989). Statistical method for testing the neutral mutation
hypothesis by DNA polymorphism. *Genetics* 123:585-595.

Lauring AS (2020). Within-Host Viral Diversity: A Window into Viral
Evolution. *Annu Rev Virol* 7:63-81.
[doi:10.1146/annurev-virology-010320-061642](https://doi.org/10.1146/annurev-virology-010320-061642)

## Examples

``` r
tajimaD(S = 10, n = 1000, pi_hat = 0.001, genomeLength = 10000)
#> [1] 12.88484
```
