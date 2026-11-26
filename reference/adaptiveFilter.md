# Create a technology-appropriate ISNVFilter

Returns an [`ISNVFilter`](ISNVFilter.md) with default thresholds tuned
for the error profile and read characteristics of the specified
sequencing technology.

## Usage

``` r
adaptiveFilter(technology = c("illumina", "ont", "pacbio_hifi", "generic"))
```

## Arguments

- technology:

  Character scalar. One of `"illumina"`, `"ont"`, `"pacbio_hifi"`, or
  `"generic"`.

## Value

An [`ISNVFilter`](ISNVFilter.md) object.

## Details

Each technology has distinct error modes that motivate different default
thresholds:

- **Illumina** (`minFreq = 0.03, minDepth = 100`):

  Illumina short reads have low per-base error rates (~0.1–0.5%) but
  systematic context-dependent errors. A 3% frequency threshold safely
  exceeds the noise floor while retaining true low-frequency variants.
  Grubaugh et al. (2019) measured iSNVs above 3% accurately from at
  least 1,000 RNA copies sequenced to at least 400x, so 100x is a
  permissive floor; low-frequency false positives inflate diversity
  estimates markedly (McCrone & Lauring 2016).

- **Oxford Nanopore (ONT)**
  (`minFreq = 0.10, minDepth = 60, maxStrandBias = Inf`):

  ONT reads have higher per-base error rates (~1–5% with R10 chemistry),
  requiring a more conservative frequency threshold of 10% to avoid
  false positives; Bull et al. (2020) found that ONT sequencing does not
  accurately detect variants at low read-count frequencies, and reported
  over 99% sensitivity and precision for single-nucleotide variants
  above about 60-fold depth, which sets the depth threshold. Strand bias
  filtering is disabled (`Inf`) rather than given an untested threshold;
  supply `maxStrandBias` explicitly when the data support one.

- **PacBio HiFi** (`minFreq = 0.05, minDepth = 50`):

  HiFi consensus reads reach about 99.8% read accuracy via circular
  consensus sequencing (Wenger et al. 2019), allowing a frequency
  threshold below ONT but above Illumina. Depth of 50x reflects the
  lower throughput but high per-read accuracy.

- **Generic**:

  Falls back to the standard [`ISNVFilter`](ISNVFilter.md) defaults
  (`minFreq = 0.03, minDepth = 100`), suitable for Illumina or when the
  technology is unknown.

## References

McCrone JT, Lauring AS (2016). Measurements of intrahost viral diversity
are extremely sensitive to systematic errors in variant calling. *J
Virol* 90:6884-6895.
[doi:10.1128/JVI.00667-16](https://doi.org/10.1128/JVI.00667-16)

Grubaugh ND et al. (2019). An amplicon-based sequencing framework for
accurately measuring intrahost virus diversity using PrimalSeq and iVar.
*Genome Biology* 20:8.
[doi:10.1186/s13059-018-1618-7](https://doi.org/10.1186/s13059-018-1618-7)

Bull RA et al. (2020). Analytical validity of nanopore sequencing for
rapid SARS-CoV-2 genome analysis. *Nature Communications* 11:6272.
[doi:10.1038/s41467-020-20075-6](https://doi.org/10.1038/s41467-020-20075-6)

Wenger AM et al. (2019). Accurate circular consensus long-read
sequencing improves variant detection and assembly of a human genome.
*Nature Biotechnology* 37:1155–1162.
[doi:10.1038/s41587-019-0217-9](https://doi.org/10.1038/s41587-019-0217-9)

## See also

[`ISNVFilter`](ISNVFilter.md) for the underlying constructor,
[`flagISNV`](flagISNV.md) for applying filters.

## Examples

``` r
# Illumina defaults
adaptiveFilter("illumina")
#> ISNVFilter object
#>   minDepth:      100 
#>   freq range:   [0.03, 0.5]
#>   minAltReads:   10 
#>   maxStrandBias: 10 
#>   replicateConc: none

# ONT with strand bias disabled
filt <- adaptiveFilter("ont")
slot(filt, "maxStrandBias")  # Inf
#> [1] Inf

# PacBio HiFi
adaptiveFilter("pacbio_hifi")
#> ISNVFilter object
#>   minDepth:      50 
#>   freq range:   [0.05, 0.5]
#>   minAltReads:   10 
#>   maxStrandBias: 10 
#>   replicateConc: none

# Generic (same as ISNVFilter())
adaptiveFilter("generic")
#> ISNVFilter object
#>   minDepth:      100 
#>   freq range:   [0.03, 0.5]
#>   minAltReads:   10 
#>   maxStrandBias: 10 
#>   replicateConc: none
```
