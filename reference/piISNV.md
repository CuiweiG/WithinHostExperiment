# Nucleotide diversity (pi) from iSNV frequencies

Estimates nucleotide diversity following Nei & Li (1979), adapted for
deep sequencing allele frequencies: \$\$\pi = \frac{n}{n-1} \cdot
\frac{2}{L} \sum_i p_i(1 - p_i)\$\$ where \\p_i\\ is the alternative
allele frequency at site \\i\\, \\L\\ is the genome length, and \\n\\ is
the number of reads sampled at a site, for which mean read depth stands
in. The \\n/(n-1)\\ correction removes the finite-sample bias of
estimating a frequency from \\n\\ draws (Nei 1987, Eq. 10.5). This is a
different use of "sample size" from the one in [`tajimaD`](tajimaD.md),
where reads must not stand in for independent lineages. When `meanDepth`
is not provided, the uncorrected estimator is returned (equivalent to
infinite sample size).

## Usage

``` r
piISNV(freq, genomeLength, meanDepth = NULL)
```

## Arguments

- freq:

  Numeric vector. Alternative allele frequencies, between 0 and 1;
  anything else is an error. NA values are dropped.

- genomeLength:

  Integer scalar. Reference genome length in bp.

- meanDepth:

  Numeric scalar (optional). Mean read depth for finite-sample
  correction. If `NULL`, no correction is applied.

## Value

Numeric scalar. Nucleotide diversity.

## References

Nei M, Li WH (1979). Mathematical model for studying genetic variation
in terms of restriction endonucleases. *Proc Natl Acad Sci USA*
76:5269-5273.
[doi:10.1073/pnas.76.10.5269](https://doi.org/10.1073/pnas.76.10.5269)

Nei M (1987). *Molecular Evolutionary Genetics.* Columbia University
Press, New York. Chapter 10.

## Examples

``` r
piISNV(c(0.1, 0.3, 0.5), genomeLength = 13588)
#> [1] 8.095378e-05
piISNV(c(0.1, 0.3, 0.5), genomeLength = 13588, meanDepth = 1000)
#> [1] 8.103482e-05
```
