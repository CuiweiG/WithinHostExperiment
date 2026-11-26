# Flag temporally inconsistent iSNVs in longitudinal data

For longitudinal samples from the same host, identifies iSNVs that
appear at only one timepoint (transient) versus those detected at
multiple consecutive timepoints (persistent). Transient iSNVs are more
likely to be sequencing artefacts; persistent ones are more likely
genuine within-host variants.

## Usage

``` r
flagTemporalInconsistency(
  whe,
  hostCol = "host_id",
  timeCol = "timepoint",
  minTimepoints = 2L,
  minConsecutive = NULL,
  usePassedOnly = TRUE
)
```

## Arguments

- whe:

  A [`WithinHostExperiment`](WithinHostExperiment.md) with longitudinal
  samples from the same host.

- hostCol:

  Character scalar. Column name in `colData` for host identifier
  (default: `"host_id"`).

- timeCol:

  Character scalar. Column name in `colData` for time point (default:
  `"timepoint"`). It must hold numbers, dates, date-times, time
  differences or a factor whose levels are in time order; text is read
  as numbers when every value is one, and is otherwise an error, since
  text such as `"day10"` and `"day2"` has no time order. A blank or `NA`
  is a missing timepoint.

- minTimepoints:

  Integer scalar. Minimum number of distinct timepoints at which an iSNV
  must be detected to be considered persistent (default: 2). Replicate
  samples taken at one timepoint count as one timepoint.

- minConsecutive:

  Integer scalar. Minimum number of **consecutive** timepoints required,
  counted over the distinct timepoints at which the host was sampled. If
  `NULL` (default), any `minTimepoints` detections suffice regardless of
  consecutiveness. Set to 2 or 3 for stricter filtering.

- usePassedOnly:

  Logical. If `TRUE` (default), only consider QC-passed entries; a cell
  whose `qcPass` is `NA` counts as not passed.

## Value

A [`WithinHostExperiment`](WithinHostExperiment.md) with updated
`qcPass` assay and a new column `temporal_class` in `mcols(rowRanges)`.
Sites are classified within each host; a transient site is flagged only
in the samples of the hosts where it is transient. The site-level
`temporal_class` is `"persistent"` if the site is persistent in at least
one host, `"transient"` if it is transient in every host where it was
detected, and `NA` where it was never detected.

## Details

This check applies only where a host was sampled at more than one
timepoint. Detections are counted per distinct timepoint: a site called
in any sample taken at a timepoint, including any of its technical
replicates, is detected at that timepoint, so replicates alone cannot
make a site persistent. A host sampled at a single timepoint, however
many replicates it has, is not classified. Samples whose timepoint is
`NA` are not counted, with a warning, but are flagged with the rest of
their host. Transient sites are flagged by setting `qcPass` to `FALSE`
where it was `TRUE`; `NA` cells are left as they are.

The biological rationale: genuine within-host variants arise from viral
replication and should persist (or change smoothly) across consecutive
timepoints. An iSNV that appears once at 5% and is absent at all other
timepoints is more likely a sequencing or library preparation artefact.

The function records the filtering step in `qcLog` with step name
`"temporal_consistency"`.

## References

Farjo M et al. (2024). Within-host evolutionary dynamics and tissue
compartmentalization during acute SARS-CoV-2 infection. *J. Virol.*
98:e01618-23.
[doi:10.1128/jvi.01618-23](https://doi.org/10.1128/jvi.01618-23)

## Examples

``` r
gr <- GenomicRanges::GRanges("seg1",
    IRanges::IRanges(c(100, 200, 300), width = 1))
S4Vectors::mcols(gr)$ref <- c("A", "C", "G")
S4Vectors::mcols(gr)$alt <- c("T", "G", "A")
freq <- matrix(c(0.1, NA,  0.2,   # var1: 2 timepoints
                 0.05, NA,  NA,    # var2: 1 timepoint (transient)
                 0.15, 0.2, 0.1),  # var3: 3 timepoints
               nrow = 3, ncol = 3)
whe <- WithinHostExperiment(
    assays = list(altFreq = freq),
    rowRanges = gr,
    colData = S4Vectors::DataFrame(
        sample_id = c("t1", "t2", "t3"),
        host_id = c("H1", "H1", "H1"),
        timepoint = c(1, 2, 3)))
whe <- flagTemporalInconsistency(whe, minTimepoints = 2L)
S4Vectors::mcols(SummarizedExperiment::rowRanges(whe))$temporal_class
#> [1] "persistent" "transient"  "persistent"
```
