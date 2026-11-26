# Detect potential cross-contamination between samples

Identifies iSNV sites where a high-frequency variant in one sample
appears at suspiciously low frequency in another sample from the same
sequencing run, consistent with index hopping or sample
cross-contamination.

## Usage

``` r
detectCrossContamination(
  whe,
  runCol = NULL,
  sourceMinFreq = 0.2,
  sinkMaxFreq = 0.05,
  sinkMinFreq = 0.005,
  usePassedOnly = TRUE
)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md) with multiple
  samples.

- runCol:

  Character scalar (optional). Column name in `colData` identifying the
  sequencing run/batch. If `NULL`, all samples are compared against each
  other. A name that is not a column of `colData` is an error. Samples
  whose run is `NA` are left out of every comparison, with a warning
  giving their number.

- sourceMinFreq:

  Numeric scalar. Minimum frequency in the "source" sample for a variant
  to be considered a potential contaminant origin (default: 0.20).

- sinkMaxFreq:

  Numeric scalar. Maximum frequency in the "sink" sample for a variant
  to be flagged as potential contamination (default: 0.05).

- sinkMinFreq:

  Numeric scalar. Minimum frequency in the "sink" sample (default:
  0.005). Variants below this are likely noise, not contamination.

- usePassedOnly:

  Logical. If `TRUE` (default for source), only consider QC-passed
  variants in the source sample. A source call that
  [`flagISNV`](flagISNV.md) failed for being above its `maxFreq` (0.5 by
  default) is still used when it has the read depth, alternative reads
  and strand balance that every filter in `qcFilters(whe)` asks for: a
  high frequency is the signature looked for, and the frequency ceiling
  is there for the iSNV analyses. Other QC steps do not record why a
  call failed, so a call above `maxFreq` that one of them also flagged
  is used too.

## Value

A `DataFrame` with columns:

- source_sample:

  Sample where the variant is at high frequency.

- sink_sample:

  Sample where the variant appears at suspiciously low frequency.

- chrom:

  Chromosome/segment.

- position:

  Genomic position.

- ref:

  Reference allele.

- alt:

  Alternative allele.

- source_freq:

  Frequency in the source sample.

- sink_freq:

  Frequency in the sink sample.

- freq_ratio:

  sink_freq / source_freq – lower values are more suspicious.

- run:

  Sequencing run (if `runCol` provided).

## Details

Index hopping occurs when adapters are exchanged between library
molecules during pooled sequencing, causing reads from one sample to be
misassigned to another. The signature: a variant at high frequency
(e.g., 80%) in sample A appears at very low frequency (e.g., 0.5-2%) in
sample B sequenced on the same flow cell, at a ratio consistent with the
expected hopping rate (typically 0.1-2% for Illumina platforms).

This function does not require external metadata beyond what is in
`colData`. For finer control, use the output `DataFrame` to implement
study-specific contamination criteria.

## References

Bendall EE et al. (2022). SARS-CoV-2 genomic diversity in households
highlights the challenges of sequence-based transmission inference.
*mSphere* 7:e00400-22.
[doi:10.1128/msphere.00400-22](https://doi.org/10.1128/msphere.00400-22)

Costello M et al. (2018). Characterization and remediation of sample
index swaps by non-redundant dual indexing on massively parallel
sequencing platforms. *BMC Genomics* 19:332.
[doi:10.1186/s12864-018-4703-0](https://doi.org/10.1186/s12864-018-4703-0)

## Examples

``` r
gr <- GenomicRanges::GRanges("seg",
    IRanges::IRanges(c(100, 200, 300), width = 1))
S4Vectors::mcols(gr)$ref <- c("A", "C", "G")
S4Vectors::mcols(gr)$alt <- c("T", "G", "A")
## Sample A has pos 100 at 80%; sample B has same pos at 1%
freq <- matrix(c(0.80, 0.05, 0.30,   # sample A
                 0.01, 0.40, 0.00),   # sample B
               nrow = 3, ncol = 2)
whe <- WithinHostExperiment(
    assays = list(altFreq = freq),
    rowRanges = gr,
    colData = S4Vectors::DataFrame(
        sample_id = c("A", "B"),
        run = c("run1", "run1")))
detectCrossContamination(whe, runCol = "run")
#> DataFrame with 2 rows and 10 columns
#>   source_sample sink_sample       chrom  position         ref         alt
#>     <character> <character> <character> <integer> <character> <character>
#> 1             A           B         seg       100           A           T
#> 2             B           A         seg       200           C           G
#>   source_freq sink_freq freq_ratio         run
#>     <numeric> <numeric>  <numeric> <character>
#> 1         0.8      0.01     0.0125        run1
#> 2         0.4      0.05     0.1250        run1
```
