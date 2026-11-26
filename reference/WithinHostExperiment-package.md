# WithinHostExperiment: Within-Host Pathogen Variant QC, Diversity, and Bottleneck

Provides a
[`RangedSummarizedExperiment-class`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)-based
container and workflow functions for analysing intra-host single
nucleotide variants (iSNV) from pathogen deep sequencing data.

## Details

The package fills an infrastructure gap between upstream variant callers
(iVar, LoFreq, Freebayes) and downstream transmission bottleneck
workflows. Its core contribution is replicate-aware,
uncertainty-preserving quality control – variants are flagged rather
than deleted, preserving full provenance for reproducible analysis.
Standardised export functions provide interoperability with tools such
as ViralBottleneck.

## Key function groups

- Import:

  [`readWithinHost`](readWithinHost.md) (from VCF),
  [`readWithinHostTable`](readWithinHostTable.md) (from TSV/CSV)

- Quality control:

  [`ISNVFilter`](ISNVFilter.md) (configure thresholds),
  [`flagISNV`](flagISNV.md) (apply standard QC),
  [`flagReplicateDiscordance`](flagReplicateDiscordance.md)
  (replicate-aware QC), [`flagPrimerSites`](flagPrimerSites.md)
  (amplicon primer masking), [`qcSummary`](qcSummary.md) (QC report),
  [`passedISNV`](WHE-accessors.md) (extract clean subset)

- Diversity:

  [`calcDiversity`](calcDiversity.md), [`shannonISNV`](shannonISNV.md),
  [`piISNV`](piISNV.md), [`wattersonISNV`](wattersonISNV.md),
  [`compareDiversity`](compareDiversity.md)

- Transmission analysis:

  [`transmissionPairs`](WHE-accessors.md),
  [`exportPairFrequencies`](exportPairFrequencies.md),
  [`asViralBottleneckInput`](asViralBottleneckInput.md),
  [`calcSharedVariants`](calcSharedVariants.md),
  [`quickBottleneck`](quickBottleneck.md)

- Neutrality tests:

  [`tajimaD`](tajimaD.md), [`fusFs`](fusFs.md),
  [`calcNeutralityTests`](calcNeutralityTests.md),
  [`simpsonISNV`](simpsonISNV.md), [`chao1ISNV`](chao1ISNV.md)

- Annotation:

  [`annotateGenes`](annotateGenes.md) (labels each site with the genes
  it overlaps, taking gene ranges as a `GRanges` or `TxDb`),
  [`annotateCodonChange`](annotateCodonChange.md) (the codon and amino
  acid consequence of each SNV, from
  [`VariantAnnotation::predictCoding()`](https://rdrr.io/pkg/VariantAnnotation/man/predictCoding-methods.html),
  taking CDS records as a `GRanges`, `GRangesList` or `TxDb` and the
  reference as a `DNAStringSet`)

- Selection analysis:

  [`pNpSWithinHost`](pNpSWithinHost.md) (within-host pN/pS),
  [`consensusSequences`](consensusSequences.md) (per-sample consensus
  sequences as a `DNAStringSet`, for dN/dS between hosts with MSA2dist),
  [`estimateSelectionCoefficient`](estimateSelectionCoefficient.md)

- Multi-caller consensus:

  [`consensusISNV`](consensusISNV.md),
  [`calcCallerConcordance`](calcCallerConcordance.md),
  [`benchmarkCallers`](benchmarkCallers.md)

- Longitudinal analysis:

  [`trackFrequency`](trackFrequency.md),
  [`detectEmergingVariants`](detectEmergingVariants.md),
  [`plotFrequencyTrajectory`](plotFrequencyTrajectory.md)

- Exact bottleneck estimation:

  [`exactBottleneck`](exactBottleneck.md)

- Strand bias:

  [`sciStrandBias`](sciStrandBias.md),
  [`flagStrandBiasSCI`](flagStrandBiasSCI.md)

- Temporal QC:

  [`flagTemporalInconsistency`](flagTemporalInconsistency.md)
  (longitudinal consistency filtering – flag transient iSNVs that appear
  at only one timepoint)

- Cross-contamination detection:

  [`detectCrossContamination`](detectCrossContamination.md) (identify
  index- hopping signatures between co-sequenced samples)

- Adaptive filtering:

  [`adaptiveFilter`](adaptiveFilter.md)

- Diversity modelling:

  [`fitDiversityModel`](fitDiversityModel.md),
  [`checkDiversityPlausibility`](checkDiversityPlausibility.md)

- Site frequency spectrum:

  [`WithinHostSFS-class`](WithinHostSFS-class.md) (S4 class),
  [`buildSFS`](buildSFS.md) and [`buildSFSList`](buildSFSList.md)
  (construct from WHE), [`correctSFSBias`](correctSFSBias.md)
  (ascertainment bias correction),
  [`neutralityFromSFS`](neutralityFromSFS.md) (Tajima's D + Fu's Fs +
  Fay & Wu's H from a single SFS), [`compareSFS`](compareSFS.md)
  (chi-squared SFS comparison), [`plotSFS`](plotSFS.md)
  (publication-quality SFS plot)

- Visualisation:

  [`plotFrequencySpectrum`](plotFrequencySpectrum.md),
  [`plotPairScatter`](plotPairScatter.md),
  [`plotQCDashboard`](plotQCDashboard.md)

## Scope and limitations

- **Single-site iSNV resolution**:

  This package operates at the level of individual single nucleotide
  variants (iSNVs). It does not perform haplotype phasing or
  quasispecies reconstruction, which require read-level (BAM) analysis
  that is outside the scope of `RangedSummarizedExperiment`-based
  containers. For haplotype-aware within-host analysis, consider
  CliqueSNV (Knyazev et al. 2021 *Nucleic Acids Res* 49:e102), ShoRAH
  (Zagordi et al. 2011 *BMC Bioinformatics* 12:119), or similar tools.
  Their output can be imported into `WithinHostExperiment` for
  downstream diversity and QC analysis.

- **Co-infection vs diversity**:

  Distinguishing genuine within-host evolution from co-infection
  (superinfection) requires phylogenetic analysis of full haplotypes.
  The multi-caller consensus and replicate-concordance modules reduce
  artefacts but do not resolve this distinction.

- **Selection coefficient estimation**:

  [`estimateSelectionCoefficient`](estimateSelectionCoefficient.md)
  assumes a simple deterministic model. For stochastic (drift-aware)
  estimates, consider dedicated population-genetics frameworks.

## Design philosophy

1.  **Flag, don't delete.** QC functions mark variants in the `qcPass`
    assay rather than removing rows, preserving uncertainty and enabling
    re-analysis with different thresholds.

2.  **Replicate-aware.** Technical replicate concordance is a
    first-class QC criterion, not a post-hoc check.

3.  **Interoperable, don't reinvent.** The package provides standardised
    export functions compatible with existing downstream tools rather
    than reimplementing their methods.

4.  **Bioconductor-native.** Built on `RangedSummarizedExperiment`,
    interoperable via `VRanges` coercion.

## References

Roder AE et al. (2023). Optimized quantification of intra-host viral
diversity in SARS-CoV-2 and influenza virus sequence data. *mBio*
14:e01046-23.
[doi:10.1128/mbio.01046-23](https://doi.org/10.1128/mbio.01046-23)

Bendall EE et al. (2023). Rapid transmission and tight bottlenecks
constrain the evolution of highly transmissible SARS-CoV-2 variants.
*Nat. Commun.* 14:272.
[doi:10.1038/s41467-023-36001-5](https://doi.org/10.1038/s41467-023-36001-5)

Sobel Leonard A et al. (2017). Transmission Bottleneck Size Estimation
from Pathogen Deep-Sequencing Data, with an Application to Human
Influenza A Virus. *J Virol* 91:e00171-17.
[doi:10.1128/JVI.00171-17](https://doi.org/10.1128/JVI.00171-17)

Mostefai F et al. (2024). Refining SARS-CoV-2 intra-host variation by
leveraging large-scale sequencing data. *NAR Genomics and
Bioinformatics* 6:lqae145.
[doi:10.1093/nargab/lqae145](https://doi.org/10.1093/nargab/lqae145)

Zheng B, Johnson PCD, Hughes J (2025). ViralBottleneck: an R package for
estimating viral transmission bottlenecks from deep sequencing data
using multiple methods. *Virus Evolution* 11:veaf071.
[doi:10.1093/ve/veaf071](https://doi.org/10.1093/ve/veaf071)

## See also

Useful links:

- <https://github.com/CuiweiG/WithinHostExperiment>

- Report bugs at
  <https://github.com/CuiweiG/WithinHostExperiment/issues>

## Author

**Maintainer**: Cuiwei Gao <48gaocuiwei@gmail.com>
([ORCID](https://orcid.org/0009-0002-4896-0601)) \[copyright holder\]

## Examples

``` r
# List exported functions
ls("package:WithinHostExperiment")
#>  [1] "ConcordanceRule"              "ISNVFilter"                  
#>  [3] "WithinHostExperiment"         "adaptiveFilter"              
#>  [5] "annotateCodonChange"          "annotateGenes"               
#>  [7] "asViralBottleneckInput"       "benchmarkCallers"            
#>  [9] "buildSFS"                     "buildSFSList"                
#> [11] "calcCallerConcordance"        "calcDiversity"               
#> [13] "calcNeutralityTests"          "calcSharedVariants"          
#> [15] "chao1ISNV"                    "checkDiversityPlausibility"  
#> [17] "combineCols"                  "combineRows"                 
#> [19] "compareDiversity"             "compareSFS"                  
#> [21] "consensusISNV"                "consensusSequences"          
#> [23] "correctSFSBias"               "detectCrossContamination"    
#> [25] "detectEmergingVariants"       "estimateSelectionCoefficient"
#> [27] "exactBottleneck"              "exportPairFrequencies"       
#> [29] "fitDiversityModel"            "flagISNV"                    
#> [31] "flagPrimerSites"              "flagReplicateDiscordance"    
#> [33] "flagStrandBiasSCI"            "flagTemporalInconsistency"   
#> [35] "fusFs"                        "multitestCorrection"         
#> [37] "neutralityFromSFS"            "pNpSWithinHost"              
#> [39] "passedISNV"                   "permutationTest"             
#> [41] "piISNV"                       "plotFrequencySpectrum"       
#> [43] "plotFrequencyTrajectory"      "plotPairScatter"             
#> [45] "plotQCDashboard"              "plotSFS"                     
#> [47] "qcFilters"                    "qcFilters<-"                 
#> [49] "qcLog"                        "qcLog<-"                     
#> [51] "qcSummary"                    "quickBottleneck"             
#> [53] "readWithinHost"               "readWithinHostTable"         
#> [55] "sciStrandBias"                "shannonISNV"                 
#> [57] "show"                         "simpsonISNV"                 
#> [59] "tajimaD"                      "trackFrequency"              
#> [61] "transmissionPairs"            "wattersonISNV"               
```
