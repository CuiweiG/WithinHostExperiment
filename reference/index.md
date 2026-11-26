# Package index

## Package

An overview of the package, its classes and its workflow.

- [`WithinHostExperiment-package`](WithinHostExperiment-package.md) :
  WithinHostExperiment: Within-Host Pathogen Variant QC, Diversity, and
  Bottleneck

## Classes

The container for within-host variant calls, the filter that QC applies
to it, the rule technical replicates must meet, and the site frequency
spectrum built from the calls.

- [`show(`*`<WithinHostExperiment>`*`)`](WithinHostExperiment-class.md)
  [`` `[`( ``*`<WithinHostExperiment>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](WithinHostExperiment-class.md)
  : WithinHostExperiment: Container for within-host pathogen variant
  data
- [`show(`*`<ISNVFilter>`*`)`](ISNVFilter-class.md) : ISNVFilter:
  Quality control filter for iSNV analysis
- [`show(`*`<ConcordanceRule>`*`)`](ConcordanceRule-class.md) :
  ConcordanceRule: agreement required between technical replicates
- [`show(`*`<WithinHostSFS>`*`)`](WithinHostSFS-class.md) :
  WithinHostSFS: Site Frequency Spectrum for within-host variants

## Import and conversion

Build the container from VCF files, variant tables or Bioconductor
objects, combine objects, and convert back.

- [`WithinHostExperiment()`](WithinHostExperiment.md) : Create a
  WithinHostExperiment from components
- [`readWithinHost()`](readWithinHost.md) : Import iSNV data from VCF
  files
- [`readWithinHostTable()`](readWithinHostTable.md) : Import iSNV data
  from tabular files
- [`coerce-WHE`](coerce-WHE.md)
  [`coerce,WithinHostExperiment,data.frame-method`](coerce-WHE.md)
  [`coerce,WithinHostExperiment,VRanges-method`](coerce-WHE.md)
  [`coerce,VRanges,WithinHostExperiment-method`](coerce-WHE.md)
  [`coerce,VCF,WithinHostExperiment-method`](coerce-WHE.md) : Convert
  between WithinHostExperiment, data.frame, VRanges and VCF
- [`combineRows(`*`<WithinHostExperiment>`*`)`](WHE-combine.md)
  [`combineCols(`*`<WithinHostExperiment>`*`)`](WHE-combine.md) :
  Combine WithinHostExperiment objects by rows or columns
- [`example_whe`](example_whe.md) : Example WithinHostExperiment dataset

## Accessors

The QC audit log, the filters applied, the calls that passed, and the
transmission pairs described by the sample metadata.

- [`qcLog()`](WHE-accessors.md) [`qcFilters()`](WHE-accessors.md)
  [`passedISNV()`](WHE-accessors.md)
  [`transmissionPairs()`](WHE-accessors.md)
  [`` `qcLog<-`() ``](WHE-accessors.md)
  [`` `qcFilters<-`() ``](WHE-accessors.md) : Access the QC audit log

## Quality control

Choose filters and flag, without deleting, calls that fail depth,
frequency, strand, primer or replicate criteria; report likely
cross-contamination and summarise what QC removed.

- [`ISNVFilter()`](ISNVFilter.md) : Create an ISNVFilter object
- [`ConcordanceRule()`](ConcordanceRule.md) : Create a ConcordanceRule
  object
- [`adaptiveFilter()`](adaptiveFilter.md) : Create a
  technology-appropriate ISNVFilter
- [`flagISNV()`](flagISNV.md) : Flag iSNV sites based on quality
  criteria
- [`flagPrimerSites()`](flagPrimerSites.md) : Flag iSNVs in primer
  binding regions
- [`flagReplicateDiscordance()`](flagReplicateDiscordance.md) : Flag
  iSNVs with replicate discordance
- [`sciStrandBias()`](sciStrandBias.md) : Strand Coverage Imbalance
  (SCI) score
- [`flagStrandBiasSCI()`](flagStrandBiasSCI.md) : Flag strand-biased
  iSNVs using the SCI metric
- [`detectCrossContamination()`](detectCrossContamination.md) : Detect
  potential cross-contamination between samples
- [`qcSummary()`](qcSummary.md) : Generate a QC summary report

## Annotation

Label sites with genes and with the codon and amino acid change each
substitution makes.

- [`annotateGenes()`](annotateGenes.md) : Label variant sites with the
  genes they overlap
- [`annotateCodonChange()`](annotateCodonChange.md) : Codon and amino
  acid consequence of each SNV

## Diversity

Per-sample diversity indices, their comparison between groups, a
plausibility check, and a model of diversity against covariates.

- [`calcDiversity()`](calcDiversity.md) : Compute within-host diversity
  indices
- [`compareDiversity()`](compareDiversity.md) : Compare diversity
  between sample groups
- [`checkDiversityPlausibility()`](checkDiversityPlausibility.md) :
  Check diversity estimates against heuristic plausibility ranges
- [`shannonISNV()`](shannonISNV.md) : Shannon entropy from iSNV
  frequencies
- [`piISNV()`](piISNV.md) : Nucleotide diversity (pi) from iSNV
  frequencies
- [`wattersonISNV()`](wattersonISNV.md) : Watterson's theta from iSNV
  count
- [`simpsonISNV()`](simpsonISNV.md) : Simpson diversity index for iSNV
  frequencies
- [`chao1ISNV()`](chao1ISNV.md) : Chao1 richness estimator for iSNV data
- [`fitDiversityModel()`](fitDiversityModel.md) : Fit a linear model to
  diversity data

## Site frequency spectrum

Bin allele frequencies into a spectrum, correct it for the detection
threshold, compare spectra and compute statistics from one.

- [`buildSFS()`](buildSFS.md) : Build a site frequency spectrum from a
  WithinHostExperiment
- [`buildSFSList()`](buildSFSList.md) : Build SFS for all samples in a
  WithinHostExperiment
- [`correctSFSBias()`](correctSFSBias.md) : Correct ascertainment bias
  in the within-host SFS
- [`compareSFS()`](compareSFS.md) : Compare two within-host SFS objects
- [`neutralityFromSFS()`](neutralityFromSFS.md) : Compute all neutrality
  statistics from a WithinHostSFS

## Neutrality and selection

Tajima’s D and Fu’s Fs, the within-host pN/pS ratio, consensus sequences
for between-host dN/dS, and selection coefficients from frequency
change.

- [`calcNeutralityTests()`](calcNeutralityTests.md) : Compute neutrality
  test statistics
- [`tajimaD()`](tajimaD.md) : Tajima's D statistic
- [`fusFs()`](fusFs.md) : Fu's Fs statistic
- [`pNpSWithinHost()`](pNpSWithinHost.md) : Within-host pN/pS ratio
- [`consensusSequences()`](consensusSequences.md) : Consensus sequence
  of each sample
- [`estimateSelectionCoefficient()`](estimateSelectionCoefficient.md) :
  Estimate selection coefficient from frequency trajectory

## Transmission

Shared and private variants of a donor-recipient pair, bottleneck size
estimates, and exports for other bottleneck tools.

- [`calcSharedVariants()`](calcSharedVariants.md) : Calculate shared
  variant statistics
- [`quickBottleneck()`](quickBottleneck.md) : Quick bottleneck size
  estimation (beta-binomial approximate MLE)
- [`exactBottleneck()`](exactBottleneck.md) : Exact bottleneck size
  estimation
- [`exportPairFrequencies()`](exportPairFrequencies.md) : Export
  donor-recipient frequency table
- [`asViralBottleneckInput()`](asViralBottleneckInput.md) : Export data
  for ViralBottleneck package

## Longitudinal

Follow variant frequencies through the samples of one host over time,
and flag sites seen at too few timepoints.

- [`trackFrequency()`](trackFrequency.md) : Track variant frequency
  trajectories over time
- [`detectEmergingVariants()`](detectEmergingVariants.md) : Detect
  emerging variants across time points
- [`flagTemporalInconsistency()`](flagTemporalInconsistency.md) : Flag
  temporally inconsistent iSNVs in longitudinal data

## Multi-caller

Combine and compare the calls of several variant callers on the same
samples.

- [`consensusISNV()`](consensusISNV.md) : Multi-caller consensus iSNV
  detection
- [`calcCallerConcordance()`](calcCallerConcordance.md) : Pairwise
  caller concordance matrix
- [`benchmarkCallers()`](benchmarkCallers.md) : Benchmark variant
  callers against truth or each other

## Statistics

Permutation tests and multiple-testing correction for the summaries
above.

- [`permutationTest()`](permutationTest.md) : Permutation test for
  two-group comparison
- [`multitestCorrection()`](multitestCorrection.md) : Multiple testing
  correction for p-values

## Visualisation

Frequency spectra of an object or of WithinHostSFS objects,
donor-recipient scatter plots, a QC dashboard and frequency
trajectories.

- [`plotFrequencySpectrum()`](plotFrequencySpectrum.md) : Plot iSNV
  frequency spectrum
- [`plotPairScatter()`](plotPairScatter.md) : Plot donor-recipient
  frequency scatter
- [`plotQCDashboard()`](plotQCDashboard.md) : QC dashboard plot
- [`plotSFS()`](plotSFS.md) : Plot a within-host site frequency spectrum
- [`plotFrequencyTrajectory()`](plotFrequencyTrajectory.md) : Plot
  allele frequency trajectories over time
