# Watterson's theta from iSNV count

Estimates Watterson's theta (Watterson 1975): \$\$\hat{\theta}\_W = S /
(a_n \cdot L)\$\$ where \\S\\ is the number of segregating sites, \\a_n
= \sum\_{i=1}^{n-1} 1/i\\ is the harmonic number, and \\L\\ is the
genome length, with \\n\\ the number of sampled viral genomes.

## Usage

``` r
wattersonISNV(nSites, genomeLength, meanDepth)
```

## Arguments

- nSites:

  Integer scalar. Number of segregating (iSNV) sites.

- genomeLength:

  Integer scalar. Genome length in bp.

- meanDepth:

  Numeric scalar. Sample size \\n\\: an effective number of sampled
  genomes, not raw read depth. See Details.

## Value

Numeric scalar. Watterson's theta estimate.

## Details

As [`tajimaD`](tajimaD.md) sets out, raw read depth should **not** be
passed as \\n\\ for deep sequencing data: reads are technical replicates
of a much smaller viral population, and a depth-sized \\n\\ inflates
\\a_n\\ and makes \\\theta_W\\ correspondingly small. Pass an effective
sample size, or a depth capped at one, as
[`calcNeutralityTests`](calcNeutralityTests.md) does.

## References

Watterson GA (1975). On the number of segregating sites in genetical
models without recombination. *Theor Popul Biol* 7:256-276.
[doi:10.1016/0040-5809(75)90020-9](https://doi.org/10.1016/0040-5809%2875%2990020-9)

## Examples

``` r
## n is an effective sample size, not the raw depth
wattersonISNV(15, genomeLength = 13588, meanDepth = 100)
#> [1] 0.000213219
```
