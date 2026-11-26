# Import iSNV data from tabular files

Reads pre-computed iSNV data from TSV or CSV files and constructs a
[`WithinHostExperiment`](WithinHostExperiment.md).

## Usage

``` r
readWithinHostTable(files, format = c("ivar", "generic"), colData = NULL)
```

## Arguments

- files:

  Character vector. Paths to files (one per sample).

- format:

  Character. `"ivar"` for iVar TSV output (columns: REGION, POS, REF,
  ALT, REF_DP, ALT_DP, ALT_FREQ, TOTAL_DP, ...), or `"generic"` for CSV
  with columns: chrom, pos, ref, alt, aaf, depth.

- colData:

  `DataFrame` (optional). Sample metadata with `sample_id`. If `NULL`,
  sample IDs are derived from file basenames.

## Value

A [`WithinHostExperiment`](WithinHostExperiment.md).

## Details

Only the columns named under `format` are read.

iVar's `PASS` column, the outcome of its own test of each call, is not
read, so a row with `PASS` `FALSE` is imported as a call, and
[`flagISNV`](flagISNV.md) judges it by its depth, frequency and
read-count thresholds alone. To drop such rows, remove them from the
file first, reading every column as text so that the file is written
back unchanged:


    calls <- read.delim(path, colClasses = "character")
    passed <- tempfile(fileext = ".tsv")
    write.table(calls[calls$PASS == "TRUE", ], passed, sep = "\t",
        quote = FALSE, row.names = FALSE)
    whe <- readWithinHostTable(passed, format = "ivar")

iVar's per-strand counts (`REF_RV`, `ALT_RV`) are not imported. The
strand-bias step of [`flagISNV()`](flagISNV.md) runs only when the
assays `fwdAltCount` and `revAltCount` are present, so on an imported
object it is skipped and leaves no `qcLog` row, although
[`qcFilters()`](WHE-accessors.md) records the filter's `maxStrandBias`;
[`flagStrandBiasSCI`](flagStrandBiasSCI.md) stops for want of those
assays.

Insertions and deletions are imported as rows of width 1 at their
position, with the caller's allele notation (iVar writes `+TTG` for an
insertion and `-T` for a deletion), and nothing tells them apart from
substitutions: they pass [`flagISNV()`](flagISNV.md) like any other call
and enter the diversity statistics. A VCF read by
[`readWithinHost`](readWithinHost.md) writes the same one-base deletion
as, for instance, `AT` to `A`, so the two routes do not match it by
allele. To keep substitutions only, subset to rows whose `ref` and `alt`
are one base each:


    rr <- SummarizedExperiment::rowRanges(whe)
    whe <- whe[nchar(rr$ref) == 1L & nchar(rr$alt) == 1L, ]

## Examples

``` r
tsv <- system.file("extdata", "test_donor_ivar.tsv",
    package = "WithinHostExperiment")
whe <- readWithinHostTable(tsv, format = "ivar")
whe
#> WithinHostExperiment with 15 variant sites and 1 samples
#>   assays(5): altFreq, totalDepth, refCount, altCount, qcPass
#>   qcLog: 0 filtering step(s)
#>   transmission pairs: 0 

## Only the rows iVar passed
calls <- read.delim(tsv, colClasses = "character")
passed <- tempfile(fileext = ".tsv")
write.table(calls[calls$PASS == "TRUE", ], passed, sep = "\t",
    quote = FALSE, row.names = FALSE)
readWithinHostTable(passed, format = "ivar")
#> WithinHostExperiment with 15 variant sites and 1 samples
#>   assays(5): altFreq, totalDepth, refCount, altCount, qcPass
#>   qcLog: 0 filtering step(s)
#>   transmission pairs: 0 
```
