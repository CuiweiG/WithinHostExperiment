<div align="center">

# WithinHostExperiment

*Bioconductor Infrastructure for Within-Host Pathogen Variant QC,
Diversity, and Transmission Bottleneck Workflows*

[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
[![BioC status](https://img.shields.io/badge/BioC-0.99.0-orange.svg)](https://bioconductor.org/)

</div>

---

## Why does this package exist?

Researchers studying within-host pathogen evolution run into four
practical problems:

1. **Unreplicated iSNV calls include noise** -- in the duplicate
   sequencing of Bendall *et al.* (2023), {{fig1::a::discordant calls::%d}} of
   {{data::::iSNV calls::%d}} calls ({{fig1::a::percentage of calls discordant::%.0f}}%) differ by more
   than 2 percentage points between replicates (Figure 1).
2. **Fragmented tooling** -- importing from iVar, annotating genes,
   filtering, computing diversity and exporting for bottleneck
   estimation are often handled by separate scripts (Figure 3).
3. **Longitudinal data need their own checks** -- following variant
   frequencies across sampling days calls for tracking and consistency
   flags that single samples do not (Figure 4).
4. **Hard filtering loses provenance** -- deleting variants that fail
   a threshold prevents re-analysis under other criteria. Flagging
   instead keeps every call, while downstream summaries use only those
   that pass (Figures 1c, 6).

`WithinHostExperiment` brings these steps into one
`RangedSummarizedExperiment`-based container, from import through
replicate-aware QC, diversity and longitudinal summaries to export for
bottleneck estimation.

---

## Data sources

All figures use real, publicly available data:

| Source | Data | Role |
|---------|-----------|------|
| **Bendall *et al.* (2023)** *Nat. Commun.* 14:272 | {{data::::iSNV calls::%d}} iSNVs, {{data::::samples::%d}} samples, duplicate sequencing, {{data::::transmission pairs in metadata::%d}} transmission pairs | QC, diversity, transmission, population genetics (Figs 1--3, 5--6) |
| **Farjo *et al.* (2024)** *J. Virol.* 98:e01618-23 | Participant 432870: {{fig4::::saliva samples in the source data::%d}} daily saliva samples; the {{fig4::::samples analysed (mean coverage >= 1000x)::%d}} that the study analysed (mean coverage of at least 1000x; days {{fig4::::days of infection analysed::%s}}) are used | Longitudinal within-host evolution (Fig 4) |
| **NCBI RefSeq** NC_045512.2 GFF3 | SARS-CoV-2 gene annotation | `annotateFromGFF()` demo (Figs 3b, 5b) |

All figures were generated with {{session::::R version::%s}} by
`inst/scripts/generate_readme_figures.R`, which also writes every number
quoted below to `readme_figure_values.csv`; the captions are filled from
that file by `inst/scripts/fill_readme.R`. No data were simulated.

---

## Problem 1: Unreplicated calls include noise

<div align="center">
<img src="man/figures/fig1_replicate_qc.png" width="720" alt="Replicate-aware QC impact"/>
</div>

> **Figure 1 | Replicate-aware QC and post-hoc threshold exploration.**
> **(a)** Replicate frequencies for {{data::::iSNV calls::%d}} iSNV calls:
> concordant (|frequency difference| <= 2 percentage points,
> *n* = {{fig1::a::concordant calls (|freq diff| <= 0.02)::%d}}) and discordant
> (*n* = {{fig1::a::discordant calls::%d}}); R^2 = {{fig1::a::replicate R2::%.3f}}.
> **(b)** Per-sample nucleotide diversity before (naive) and after
> replicate QC; black diamonds mark the medians
> ({{fig1::b::median pi naive (x1e-4)::%.2f}} and {{fig1::b::median pi QC (x1e-4)::%.2f}} x 10^-4). QC only removes
> calls, so pi cannot rise; it fell in {{fig1::b::samples in which QC lowered pi::%d}} of
> {{fig1::b::samples::%d}} samples, and the median fall, across the samples whose naive
> estimate is above zero, is {{fig1::b::median per-sample fall in pi (%)::%.0f}}%.
> **(c)** Median pi as the concordance threshold varies, recomputed from
> the per-call frequency difference that QC keeps rather than by
> re-importing and re-filtering the data; the dashed lines mark 2% and 5%.

<div align="center">
<img src="man/figures/fig2_frequency_spectrum.png" width="480" alt="iSNV frequency spectrum"/>
</div>

> **Figure 2 | Where discordant calls fall.**
> Mean replicate frequencies of the {{data::::iSNV calls::%d}} calls, which fall at
> {{fig2::::distinct sites::%d}} distinct variants (position and substitution):
> {{fig2::::discordant calls below 10% mean frequency::%d}} of the
> {{fig1::a::discordant calls::%d}} discordant calls sit below 10%, the range in which
> sequencing and amplification errors are hardest to tell apart from
> genuine minority variants.

---

## Problem 2: One package, complete workflow

<div align="center">
<img src="man/figures/fig3_qc_overview.png" width="700" alt="From data to annotated transmission pairs"/>
</div>

> **Figure 3 | From raw data to annotated transmission pairs.**
> **(a)** Mean read depth, averaged over the two replicates, for
> {{fig3::a::samples with depth::%d}} samples (median {{fig3::a::median mean read depth::comma}}x).
> **(b)** Concordant iSNV calls by gene, annotated with
> `annotateFromGFF()` from the NCBI GFF3.
> **(c)** Pair {{fig3::c::example pair::%s}}: {{fig3::c::donor iSNVs not detected in the recipient::%d}} of
> {{fig3::c::donor iSNVs in the example pair::%d}} concordant donor iSNVs have no iSNV call in the
> recipient.
> **(d)** Of the {{fig3::d::pairs with at least one concordant donor iSNV::%d}} pairs whose donor carried at least one
> concordant iSNV ({{data::::transmission pairs in metadata::%d}} pairs in the metadata),
> {{fig3::d::pairs sharing no donor iSNV::%d}} ({{fig3::d::percentage sharing none::%.0f}}%) share none with the recipient.
> Recipient calls are limited to frequencies of 2--98%, so a donor
> variant that became fixed in the recipient also counts as not shared;
> the panel describes detection and is not a bottleneck estimate
> (Bendall *et al.* 2023 report those).

---

## Problem 3: Longitudinal tracking needs infrastructure

<div align="center">
<img src="man/figures/fig4_diversity_landscape.png" width="700" alt="Longitudinal within-host dynamics"/>
</div>

> **Figure 4 | Within-host dynamics across sampling days.**
> Daily saliva samples from one SARS-CoV-2 participant (Farjo *et al.*
> 2024), limited to the {{fig4::::samples analysed (mean coverage >= 1000x)::%d}} samples the study analysed (days
> {{fig4::::days of infection analysed::%s}}); a call needs at least 1000 reads and a frequency
> between 3% and 97%.
> **(a)** Richness (bars) and pi (red line): QC-passed iSNVs peak at
> {{fig4::a::peak QC-passed iSNVs::%d}} on day {{fig4::a::day of peak iSNV count::%d}}, and pi on day {{fig4::a::day of maximum pi::%d}}.
> **(b)** Frequency trajectories of the {{fig4::b::variants with trajectories shown::%d}} variants with the
> widest frequency range, tracked with `trackFrequency()`; the dotted
> line marks the 3% threshold.
> **(c)** Folded spectra from `buildSFS()` on days {{fig4::c::days shown::%s}} -- the
> first, peak and last sampled days, which merge when they coincide.
> **(d)** `flagTemporalInconsistency()`: {{fig4::d::transient sites (1 timepoint)::%d}} of {{fig4::d::variant sites classified::%d}} variant
> sites ({{fig4::d::percentage transient::%.0f}}%) are detected on a single day. A site can also be
> missed on a day when its depth is below 1000 reads, so transient sites
> are flagged for review rather than treated as errors.

---

## Evolutionary summaries from the same flagged data

<div align="center">
<img src="man/figures/fig5_evolutionary_analysis.png" width="700" alt="Population-genetic summaries"/>
</div>

> **Figure 5 | Population-genetic summaries of the flagged data.**
> **(a)** Tajima's D for the {{fig5::a::samples with >= 1 concordant iSNV::%d}} samples with at least one
> concordant iSNV and a finite D, using read depth capped at 100 as the sample size:
> median {{fig5::a::median Tajima D (QC, n capped at 100)::%.2f}}, {{fig5::a::percentage of samples with D < 0::%.0f}}% below zero. Calls below the 2% frequency
> threshold are absent, and negative D is expected both under purifying
> selection and after population growth, so these values describe the
> spectra rather than test either.
> **(b)** Nei-Gojobori ratio of nonsynonymous to synonymous iSNVs per
> site, pooling concordant calls across samples
> ({{fig5::b::pooled nonsynonymous iSNVs::%d}} nonsynonymous, {{fig5::b::pooled synonymous iSNVs::%d}} synonymous; pooled ratio {{fig5::b::pooled dN/dS (Jukes-Cantor)::%.2f}}).
> A further {{fig5::b::concordant iSNVs outside CDS or unannotated::%d}} concordant calls lie outside the annotated
> coding sequences, or carry no codon annotation, and enter neither count.
> Per-gene ratios rest on a handful of calls and are shown with their
> counts; a gene with no synonymous call gives no ratio and appears in
> grey with its nonsynonymous count. None of them is evidence of
> selection on any gene.
> **(c)** Per-sample pi ranked by the naive estimate; segments show the
> reduction from QC.
> **(d)** Naive, QC and bias-corrected spectra as `WithinHostSFS`
> objects from `buildSFS()` and `correctSFSBias()` (binomial correction).

---

## Problem 4: QC should remove noise, not signal

<div align="center">
<img src="man/figures/fig6_consensus_validation.png" width="700" alt="What replicate QC changes"/>
</div>

> **Figure 6 | What replicate QC changes.**
> **(a)** Tajima's D before and after QC in the {{fig6::a::samples with finite D before and after QC::%d}} samples with a
> finite value in both: median {{fig6::a::median Tajima D naive::%.2f}} before QC and {{fig6::a::median Tajima D QC::%.2f}} after.
> Across all of those samples the paired median change is
> {{fig6::a::median paired change in D (QC minus naive)::%.2f}} (95% CI {{fig6::a::95% CI lower (order statistics)::%.2f}} to {{fig6::a::95% CI upper (order statistics)::%.2f}}); QC altered D in
> {{fig6::a::samples whose D changed under QC::%d}} samples, and among those the median change is
> {{fig6::a::median change in D among samples QC altered::%.2f}} (95% CI {{fig6::a::95% CI lower, altered samples::%.2f}} to {{fig6::a::95% CI upper, altered samples::%.2f}}). Both intervals are
> distribution-free, from binomial order statistics.
> **(b)** Frequency spectra of the calls QC kept and the calls it removed
> ({{fig6::b::calls kept by QC in the spectrum::%d}} and {{fig6::b::calls removed by QC in the spectrum::%d}} of them fall inside the binned
> frequency range), compared with `compareSFS()`:
> chi-squared = {{fig6::b::compareSFS chi-squared (kept vs removed)::%.1f}}, df = {{fig6::b::compareSFS df::%d}}, asymptotic *p* {{fig6::b::compareSFS asymptotic p::p}} and
> Monte Carlo *p* from 10^5 tables {{fig6::b::Monte Carlo p (1e5 tables)::p}}; the smallest expected count is
> {{fig6::b::smallest expected count::%.1f}}. Calls pooled across samples are not independent draws, so the
> comparison is descriptive.

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
   thresholds.
2. **Replicate-aware.** Technical replicate concordance is a
   first-class QC criterion, not a post-hoc script.
3. **Interoperable, don't reinvent.** Standardised export to
   ViralBottleneck (Zheng *et al.* 2025) and VRanges-based
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

All README figures are generated from public data by one script, and the
captions are then filled from the values it records:

```r
# Clone the repository, then from the package root directory:
source("inst/scripts/generate_readme_figures.R")
source("inst/scripts/fill_readme.R")
```

Data files in `inst/scripts/real_data/` (GitHub only, not in installed package):

| File | Source | Description |
|------|--------|-------------|
| `all_variants_filtered.tsv` | Bendall *et al.* 2023 | {{data::::iSNV calls::%d}} iSNV calls, {{data::::samples::%d}} samples, both replicates |
| `AvgCoverage.all` | Bendall *et al.* 2023 | Per-replicate mean amplicon depth |
| `Transmission_pairs.csv` | Bendall *et al.* 2023 | Household transmission pair metadata |
| `farjo_longitudinal/` | Farjo *et al.* 2024 | Daily saliva iVar TSVs for participant 432870, with `samples_432870.csv` (dates, coverage and study inclusion, built by `inst/scripts/prepare_farjo_metadata.R` from github.com/BROOKELAB/SARS-CoV-2-within-host-evolution) |
| `sars2_NC045512.gff3` | NCBI RefSeq | SARS-CoV-2 gene annotation (NC_045512.2) |

## Key references

- Farjo M *et al.* (2024) Within-host evolutionary dynamics and
  tissue compartmentalization during acute SARS-CoV-2 infection.
  *J. Virol.* 98:e01618-23.
- Roder AE *et al.* (2023) Optimized quantification of intra-host
  viral diversity in SARS-CoV-2 and influenza virus sequence data.
  *mBio* 14:e01046-23.
- Bendall EE *et al.* (2023) Rapid transmission and tight bottlenecks
  constrain the evolution of highly transmissible SARS-CoV-2 variants.
  *Nat. Commun.* 14:272.
- Sobel Leonard A *et al.* (2017) Transmission bottleneck size
  estimation from pathogen deep-sequencing data. *J Virol*
  91:e00171-17.
- Mostefai F *et al.* (2024) Refining SARS-CoV-2 intra-host variation
  by leveraging large-scale sequencing data. *NAR Genom. Bioinform.*
  6:lqae145.
- Zheng B, Johnson PCD, Hughes J (2025) ViralBottleneck: an R package for
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

The continuous-integration configuration lives on the `ci` branch; the default
branch holds the package source only.
