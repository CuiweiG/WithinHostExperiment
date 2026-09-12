<div align="center">

# WithinHostExperiment

*Bioconductor Infrastructure for Within-Host Pathogen Variant QC,
Diversity, and Transmission Bottleneck Workflows*

[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
[![BioC status](https://img.shields.io/badge/BioC-0.99.0-orange.svg)](https://bioconductor.org/)

</div>

---

## Why does this package exist?

Researchers studying within-host pathogen evolution face four
practical problems that no existing Bioconductor package solves:

1. **False-positive iSNVs inflate diversity estimates** -- without
   replicate-aware QC, a third of variant calls can be noise
   (Figure 1).
2. **Fragmented tooling** -- importing from iVar, annotating genes,
   filtering, computing diversity, and exporting for bottleneck
   estimation each require separate scripts (Figure 3).
3. **No longitudinal tracking infrastructure** -- following variant
   frequencies across timepoints requires custom code for every
   study (Figure 4).
4. **Hard-filtering destroys provenance** -- deleting variants that
   fail a threshold makes it impossible to re-analyse with different
   criteria. The flag-don't-delete design preserves all data while
   removing noise from downstream results (Figures 1c, 6).

`WithinHostExperiment` solves all four in a single
`RangedSummarizedExperiment`-based container with 50 exported
functions covering the complete workflow from import to publication.

---

## Data sources

All figures use real, publicly available data:

| Dataset | Reference | Role |
|---------|-----------|------|
| **McCrone *et al.* (2023)** *Nat. Commun.* 14:235 | 159 iSNVs, 83 samples, duplicate sequencing, 132 transmission pairs | QC, diversity, transmission, population genetics (Figs 1--3, 5--6) |
| **Farjo *et al.* (2022)** bioRxiv, Brooke Lab | 1 patient, **9 timepoints**, iVar output | Longitudinal within-host evolution (Fig 4) |
| **NCBI RefSeq** NC_045512.2 GFF3 | SARS-CoV-2 gene annotation | `annotateFromGFF()` demo (Figs 3b, 4d) |

All figures were generated locally in R 4.5.3 via
`inst/scripts/generate_readme_figures.R`. No data were simulated.

---

## Problem 1: Your iSNVs are noisier than you think

<div align="center">
<img src="man/figures/fig1_replicate_qc.png?r20260330155119" width="720" alt="Replicate-aware QC impact"/>
</div>

> **Figure 1 | Replicate-aware QC and post-hoc threshold exploration.**
> **(a)** Replicate-frequency scatter for 159 iSNV calls; concordant
> (|freq diff| <= 2%, *n* = 106, green) vs. discordant (*n* = 53,
> red); R^2 = 0.948.
> **(b)** Per-sample nucleotide diversity before (naive) and after
> replicate QC; black diamonds = medians; Wilcoxon *p* = 4.0 x 10^-8.
> **(c)** The flag-don't-delete advantage: median pi as a function of
> concordance threshold, computed post-hoc from a single flagged
> dataset -- no pipeline re-run needed.

<div align="center">
<img src="man/figures/fig2_frequency_spectrum.png?r20260330155119" width="480" alt="iSNV frequency spectrum"/>
</div>

> **Figure 2 | Where the noise lives.**
> Discordant calls (red) cluster below 10% frequency -- exactly where
> sequencing errors masquerade as genuine minority variants.

---

## Problem 2: One package, complete workflow

<div align="center">
<img src="man/figures/fig3_qc_overview.png?r20260330155119" width="700" alt="From data to biology"/>
</div>

> **Figure 3 | From raw data to biological discovery.**
> **(a)** Sequencing depth (median 2,033x).
> **(b)** Gene-level iSNV distribution annotated via
> `annotateFromGFF()` with NCBI GFF3.
> **(c)** Pair HH46: all 5 donor iSNVs lost at transmission,
> with gene and frequency annotation.
> **(d)** Variant sharing across 52 transmission pairs: **90%
> share zero iSNVs** -- tight bottleneck evidence.

---

## Problem 3: Longitudinal tracking needs infrastructure

<div align="center">
<img src="man/figures/fig4_diversity_landscape.png?r20260330155119" width="700" alt="Longitudinal evolution"/>
</div>

> **Figure 4 | Nine days of within-host evolution.**
> One SARS-CoV-2 patient sampled across 9 timepoints
> (Farjo *et al.* 2022).
> **(a)** Diversity arc: richness (bars) and pi (red line) peak at
> day 8 with 286 QC-passed iSNVs.
> **(b)** Frequency trajectories of the 10 most dynamic variants --
> each line is one iSNV tracked by `trackFrequency()`.
> **(c)** SFS evolution via `buildSFS()`: the frequency spectrum
> shape shifts across infection stages (day 1, peak, day 9).
> **(d)** Temporal QC via `flagTemporalInconsistency()`: **97% of
> iSNVs are transient** (detected at only one timepoint) -- a new
> QC dimension unique to longitudinal data.

---

## Evolutionary analysis from the same flagged data

<div align="center">
<img src="man/figures/fig5_evolutionary_analysis.png?r20260330155119" width="700" alt="Evolutionary analysis"/>
</div>

> **Figure 5 | Population genetics: selection or drift?**
> **(a)** Tajima's D (n capped at 100): median = -0.57, 65% negative
> -- consistent with purifying selection.
> **(b)** Per-gene dN/dS: ORF1b = 0.6 (purifying), ORF3a/ORF7a = 2.0
> (positive selection signal); Intergenic = undefined (no synonymous
> sites).
> **(c)** Per-sample pi ranked by naive estimate; connecting segments
> show QC-induced reduction.
> **(d)** SFS as first-class object: naive (orange), QC (green),
> and bias-corrected (blue) spectra via `buildSFS()` and
> `correctSFSBias()` -- the first structured within-host SFS
> implementation in Bioconductor.

---

## Problem 4: QC should remove noise, not signal

<div align="center">
<img src="man/figures/fig6_consensus_validation.png?r20260330155119" width="700" alt="QC signal preservation"/>
</div>

> **Figure 6 | Replicate QC removes noise without distorting signal.**
> **(a)** Tajima's D distributions before and after QC; Wilcoxon
> *p* = n.s., confirming QC preserves population-genetic inference.
> **(b)** Formal SFS comparison via `compareSFS()`: chi-squared
> test (p = 0.97) confirms QC preserves the frequency spectrum
> shape -- the first quantitative SFS-level QC validation.

---

## Installation

```r
# After Bioconductor acceptance:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("WithinHostExperiment")

# Development version from GitHub:
BiocManager::install("CuiweiG/WithinHostExperiment")
```

## Quick start

```r
library(WithinHostExperiment)

# 1. Import donor and recipient VCFs
vcf_d <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
vcf_r <- system.file("extdata", "test_recipient.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(c(vcf_d, vcf_r),
    colData = S4Vectors::DataFrame(
        sample_id = c("donor", "recipient"),
        role      = c("donor", "recipient"),
        pair_id   = c("pair_1", "pair_1")),
    caller = "ivar")

# 2. QC -- flag, don't delete
whe <- flagISNV(whe, ISNVFilter(minDepth = 200L, minFreq = 0.03))

# 3. Diversity (finite-sample corrected pi)
div <- calcDiversity(passedISNV(whe), genomeLength = 29903L)

# 4. Transmission bottleneck
nb <- quickBottleneck(whe, pairId = "pair_1")
```

---

## Key design principles

1. **Flag, don't delete.** QC marks variants in `qcPass` rather than
   removing rows, preserving all data for re-analysis under different
   thresholds (Cavallo *et al.* 2023).
2. **Replicate-aware.** Technical replicate concordance is a
   first-class QC criterion, not a post-hoc script.
3. **Interoperable, don't reinvent.** Standardised export to
   ViralBottleneck (Archaman *et al.* 2025) and VRanges-based
   Bioconductor workflows.
4. **Bioconductor-native.** Built on `RangedSummarizedExperiment`;
   subsetting, combining, and accessors follow Bioconductor
   conventions.

## Relationship to existing tools

| Tool | Role | Relationship |
|------|------|-------------|
| iVar, LoFreq, Freebayes | Variant calling | Upstream -- we import their output |
| deepSNV | Low-frequency variant detection | Upstream -- complementary |
| ViralBottleneck | Bottleneck estimation (6 methods) | Downstream -- we export to it; also self-contained via `exactBottleneck()` |
| QSutils | Quasispecies diversity | Parallel -- different granularity (haplotype vs iSNV) |
| CliqueSNV, ShoRAH | Haplotype phasing | Complementary -- their output can feed into WHE |

## Reproducibility

All README figures are generated from public data via a single script:

```r
# Clone the repository, then from the package root directory:
source("inst/scripts/generate_readme_figures.R")
```

Data files in `inst/scripts/real_data/` (GitHub only, not in installed package):

| File | Source | Description |
|------|--------|-------------|
| `all_variants_filtered.tsv` | McCrone *et al.* 2023 | 159 iSNV calls, 83 samples, both replicates |
| `AvgCoverage.all` | McCrone *et al.* 2023 | Per-replicate mean amplicon depth, 188 samples |
| `Transmission_pairs.csv` | McCrone *et al.* 2023 | Household transmission pair metadata |
| `farjo_longitudinal/` | Farjo *et al.* 2022 | 9-timepoint iVar TSVs, patient 432870 |
| `sars2_NC045512.gff3` | NCBI RefSeq | SARS-CoV-2 gene annotation (NC_045512.2) |

## Key references

- Farjo M *et al.* (2022) Within-host evolutionary dynamics and
  tissue compartmentalization during acute SARS-CoV-2 infection.
  *bioRxiv*. doi:10.1101/2022.06.21.497047.
- Cavallo I *et al.* (2023) Optimized quantification of intra-host
  viral diversity. *mSphere* 8:e00173-23.
- McCrone JT *et al.* (2023) Rapid transmission and tight bottlenecks
  constrain the evolution of highly transmissible SARS-CoV-2 variants.
  *Nat. Commun.* 14:235.
- Sobel Leonard A *et al.* (2017) Transmission bottleneck size
  estimation from pathogen deep-sequencing data. *J Virol*
  91:e00171-17.
- Farkas C *et al.* (2024) Refining SARS-CoV-2 intra-host variation
  calling. *NAR Genomics Bioinformatics* 6:lqae145.
- Archaman B *et al.* (2025) ViralBottleneck: an R package for
  estimating viral transmission bottlenecks. *Virus Evolution*
  11:veaf071.
- Tajima F (1989) Statistical method for testing the neutral mutation
  hypothesis by DNA polymorphism. *Genetics* 123:585-595.
- Fu YX (1997) Statistical tests of neutrality of mutations against
  population growth. *Genetics* 147:915-925.
- Nei M, Gojobori T (1986) Simple methods for estimating the numbers
  of synonymous and nonsynonymous substitutions. *Mol Biol Evol*
  3:418-426.

## Documentation

- [Vignette: Introduction to WithinHostExperiment](vignettes/WithinHostExperiment.Rmd)
- [GitHub repository](https://github.com/CuiweiG/WithinHostExperiment)

The continuous-integration configuration lives on the `ci` branch, so that the
default branch holds package code only, as the _Bioconductor_ package
submission instructions require.
