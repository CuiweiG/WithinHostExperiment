# Estimate selection coefficient from frequency trajectory

Given a vector of allele frequencies observed over equally-spaced time
points, estimates the selection coefficient \\s\\ using the change in
logit-transformed frequency per unit time: \$\$s = \frac{\Delta
\mathrm{logit}(p)}{\Delta t}\$\$ where \\\mathrm{logit}(p) = \log(p /
(1 - p))\\.

## Usage

``` r
estimateSelectionCoefficient(freqs, times = NULL, generation_time = 1)
```

## Arguments

- freqs:

  Numeric vector. Allele frequencies at successive time points (values
  must be in (0, 1)).

- times:

  Numeric vector (optional). Time values for each observation. If
  `NULL`, integer time steps 1, 2, ... are used.

- generation_time:

  Numeric scalar. Time per generation (default: 1). Used to scale the
  selection coefficient.

## Value

Numeric vector of selection coefficient estimates, one per consecutive
time interval. Length is `length(freqs) - 1`.

## Details

This assumes a simple deterministic selection model. Each consecutive
pair of time points yields one estimate of \\s\\.

## References

Feder AF et al. (2016). More effective drugs lead to harder selective
sweeps in the evolution of drug resistance in HIV-1. *eLife* 5:e10670.
[doi:10.7554/eLife.10670](https://doi.org/10.7554/eLife.10670)

## Examples

``` r
freqs <- c(0.05, 0.10, 0.25, 0.50)
estimateSelectionCoefficient(freqs)
#> [1] 0.7472144 1.0986123 1.0986123
estimateSelectionCoefficient(freqs, times = c(0, 3, 7, 14))
#> [1] 0.2490715 0.2746531 0.1569446
```
