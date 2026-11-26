# Fu's Fs statistic

Computes Fu's Fs (Fu 1997), a test of neutrality based on the
probability of observing at least the observed number of haplotypes
given the estimate of theta from nucleotide diversity. Under the
infinite-sites model: \$\$F_s = \ln\\\left(\frac{S'}{1 - S'}\right)\$\$
where \\S'\\ is the probability of observing a number of haplotypes
greater than or equal to the observed value, given \\\theta =
\hat{\pi}\\.

## Usage

``` r
fusFs(S, n, theta_pi)
```

## Arguments

- S:

  Integer scalar. Number of segregating sites; the haplotype count is
  taken to be \\S + 1\\ (see Details).

- n:

  Numeric scalar. Sample size. As in [`tajimaD`](tajimaD.md), do **not**
  use raw read depth: use an effective sample size that reflects the
  number of independently sampled genomes.

- theta_pi:

  Numeric scalar. Theta estimated from nucleotide diversity (\\\hat{\pi}
  \times L\\).

## Value

Numeric scalar. Fu's Fs, or `NA` with a warning in three cases: at \\S =
0\\, where \\S' = 1\\ and Fs is infinite; at `theta_pi = 0` with \\S \>
0\\, where there is no diversity to condition on; and when \\S + 1 \>
n\\, since the statistic uses \\S + 1\\ as the number of haplotypes,
which cannot exceed the number of sequences sampled.

## Details

**The haplotype count is not observed.** Per-site calls cannot say how
many distinct haplotypes a sample holds, so the function sets that
number to \\S + 1\\. Under infinite sites without recombination this is
the largest number of haplotypes \\S\\ segregating sites can produce,
and it is usually far above the number present. The value therefore
carries no haplotype information: it depends only on \\S\\, \\n\\ and
`theta_pi`, the summaries Tajima's D already contrasts, and under
neutrality it is pushed below 0, increasingly so as \\\theta\\ grows.
With \\n = 100\\ and \\S\\ and `theta_pi` at their neutral expectations,
it is -1.46 at \\\theta = 2\\ (\\S = 10\\), -7.00 at \\\theta = 5\\ (\\S
= 26\\) and -25.0 at \\\theta = 10\\ (\\S = 52\\), where the expected
numbers of haplotypes are 8.4, 15.7 and 24.4 rather than 11, 27 and 53;
at \\\theta = 1\\ (\\S = 5\\) it is -0.37, and at very small \\\theta\\
it can be positive (+0.55 at \\\theta = 0.2\\, \\S = 1\\). A negative
value is therefore not evidence of population expansion or hitchhiking,
and no fixed threshold calibrates it.

The probability is computed exactly, for every `n`, from Ewens' sampling
formula, \\P(K = k) = \|s(n, k)\| \theta^k / \theta\_{(n)}\\, with the
unsigned Stirling numbers of the first kind built in log space. \\S'\\
is held below `1 - .Machine$double.eps`, so the value cannot exceed
about 36.

## References

Fu YX (1997). Statistical tests of neutrality of mutations against
population growth, hitchhiking and background selection. *Genetics*
147:915-925.

## Examples

``` r
fusFs(S = 10, n = 100, theta_pi = 5)
#> [1] 3.014689
```
