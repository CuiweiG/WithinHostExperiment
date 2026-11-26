# Import iSNV data from VCF files

Reads one or more VCF files produced by variant callers (iVar through a
VCF converter, LoFreq, Freebayes, or others) with
[`VariantAnnotation::readVcf()`](https://rdrr.io/pkg/VariantAnnotation/man/readVcf-methods.html)
and constructs a [`WithinHostExperiment`](WithinHostExperiment.md). A
VCF already read into R is converted with
`as(vcf, "WithinHostExperiment")`, and a `VRanges` with
`as(vr, "WithinHostExperiment")`; see [`coerce-WHE`](coerce-WHE.md).

## Usage

``` r
readWithinHost(vcfFiles, colData, caller = "auto", genome = NULL)
```

## Arguments

- vcfFiles:

  Character vector. Paths to VCF files (one per sample).

- colData:

  `DataFrame` or `data.frame`. Sample metadata. Must contain
  `sample_id`. May also contain `host_id`, `role`, `pair_id`,
  `replicate_group`, `replicate_id`, etc.

- caller:

  Character scalar. Which field layout the VCF uses: `"ivar"`,
  `"lofreq"`, `"freebayes"`, `"auto"` (default) or `"generic"`. See
  Details.

- genome:

  Optional `DNAStringSet` holding the reference sequence, as returned by
  [`Biostrings::readDNAStringSet()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html).
  The length of each sequence is recorded in `seqlengths(rowRanges())`,
  whether or not a call lies on it. It must name the sequences the calls
  lie on and be long enough to hold them; see Details.

## Value

A [`WithinHostExperiment`](WithinHostExperiment.md).

## Details

Each VCF file corresponds to one sample. The `colData` argument provides
sample metadata and must contain a `sample_id` column whose length
matches the number of files. A file holding more than one sample is an
error; read it with `readVcf()` and convert the object instead.

`caller = "ivar"` refers to a VCF carrying iVar's per-sample fields
(`REF_DP`, `ALT_DP`, `ALT_FREQ`, `ALT_QUAL`). `ivar variants` itself
writes a TSV, not a VCF, so such files come from a converter, for
instance `ivar_variants_to_vcf.py` in the nf-core/viralrecon pipeline.
An iVar TSV is read directly by
[`readWithinHostTable`](readWithinHostTable.md)`(format = "ivar")`.

**Caller.** With `caller = "auto"` the layout is chosen from the VCF
header: a header naming iVar, LoFreq or Freebayes (as the `##source`
line of each does) selects that layout, and any other gives `"generic"`.
The layouts read: iVar, FORMAT `ALT_FREQ`, `ALT_DP` and `REF_DP`, with
depth their sum (INFO `DP` where they are missing); LoFreq, INFO `AF`
and `DP`, with counts derived from them; Freebayes, FORMAT `AO`, `RO`
and `DP` (INFO `DP` where FORMAT has none), with frequency `AO / DP`;
generic, the first of those field names that the file carries.

**Several alternative alleles.** A record with several alternative
alleles is split into one site per allele with
`VariantAnnotation::expand()`. A per-allele value is split correctly
only when the header declares it `Number=A`, as the VCF specification
requires; `readVcf()` keeps just the first value of a field declared
otherwise. When such a record carries one of the per-allele fields read
here (INFO `AF`; FORMAT `AO`, `AF`, `ALT_DP`, `ALT_FREQ`, `ALT_QUAL`)
under another `Number`, the function stops rather than give the second
allele the first allele's value. Correct the header, or split the
records first (for example with `bcftools norm -m-`). On a biallelic
record a field declared `Number=1` or `Number=.` is read as it is. A
field declared `Number=R` stops the import whatever the records, since
it gives the reference allele's value first and that value would be read
as the alternative allele's; declare it `Number=A`, with one value per
alternative allele.

**Sequence lengths.** Lengths are recorded in `seqlengths(rowRanges())`
from the VCF header's `##contig` lines and from `genome`; two sources
that disagree are an error. Every sequence a header line or `genome`
names is recorded, including those without a call, since the genome
length sums over all of them. Functions that need the genome length read
it from there, so `genome` is checked before anything is recorded: it is
an error when it names none of the sequences the calls lie on, or when a
call lies past the end of its sequence, since either way it is not the
reference the calls were made against. A sequence the calls lie on that
`genome` does not name is left without a length, with a warning.

**What is not read.** The `FILTER` column is not read, so a record the
caller rejected (LoFreq's `sb_fdr` or `min_snvqual`, for instance) is
imported as a call, and [`flagISNV`](flagISNV.md) judges it by its
depth, frequency and read-count thresholds alone. To drop such records,
read the file with `readVcf()`, keep those whose `FILTER` is `PASS` or
`.`, write them to a file and import that:


    vcf <- VariantAnnotation::readVcf(path)
    keep <- VariantAnnotation::filt(vcf) %in% c("PASS", ".")
    passed <- tempfile(fileext = ".vcf")
    VariantAnnotation::writeVcf(vcf[keep, ], passed)
    whe <- readWithinHost(passed, colData, caller = "lofreq")

This keeps every argument of `readWithinHost()`. Converting the subset
directly, with `as(vcf[keep, ], "WithinHostExperiment")`, also works,
but the caller is then taken from the `##source` header line (the
generic layout when it names none), the sample names from the VCF, and
no `genome` can be given.

Per-strand counts, such as LoFreq's `DP4` or Freebayes's `SAF` and
`SAR`, are not imported. The strand-bias step of
[`flagISNV()`](flagISNV.md) runs only when the assays `fwdAltCount` and
`revAltCount` are present, so on an imported object it is skipped and
leaves no `qcLog` row, although [`qcFilters()`](WHE-accessors.md)
records the filter's `maxStrandBias`;
[`flagStrandBiasSCI`](flagStrandBiasSCI.md) stops for want of those
assays and of `fwdTotalCount` and `revTotalCount`. Assays added with
`assay(x, "fwdAltCount") <-` and so on before
[`flagISNV()`](flagISNV.md) are used.

Insertions and deletions are imported as rows of width 1 at their `POS`,
with `ref` and `alt` as the VCF writes them (`CTA` and `C` for a
deletion of two bases), and nothing tells them apart from substitutions:
they pass [`flagISNV()`](flagISNV.md) like any other call and enter the
diversity statistics. An iVar TSV read by
[`readWithinHostTable()`](readWithinHostTable.md) writes the same
deletion in its own notation, so the two routes do not match it by
allele. To keep substitutions only, subset to rows whose `ref` and `alt`
are one base each:


    rr <- SummarizedExperiment::rowRanges(whe)
    whe <- whe[nchar(rr$ref) == 1L & nchar(rr$alt) == 1L, ]

When technical replicates are not available, Roder et al. (2023, *mBio*
14:e01046-23) recommend using a combination of multiple variant callers
with stringent cutoffs to reduce false positives. To implement
multi-caller consensus, import VCFs from each caller separately and
intersect the resulting objects using standard `GRanges` overlap
operations.

## Examples

``` r
vcf <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(vcf,
    colData = S4Vectors::DataFrame(sample_id = "donor"),
    caller = "ivar")
whe
#> WithinHostExperiment with 15 variant sites and 1 samples
#>   assays(5): altFreq, totalDepth, refCount, altCount, qcPass
#>   qcLog: 0 filtering step(s)
#>   transmission pairs: 0 

## The same file read into R first, then converted
v <- VariantAnnotation::readVcf(vcf)
as(v, "WithinHostExperiment")
#> WithinHostExperiment with 15 variant sites and 1 samples
#>   assays(5): altFreq, totalDepth, refCount, altCount, qcPass
#>   qcLog: 0 filtering step(s)
#>   transmission pairs: 0 

## Only the records the caller passed (FILTER is PASS or ".")
keep <- VariantAnnotation::filt(v) %in% c("PASS", ".")
as(v[keep, ], "WithinHostExperiment")
#> WithinHostExperiment with 15 variant sites and 1 samples
#>   assays(5): altFreq, totalDepth, refCount, altCount, qcPass
#>   qcLog: 0 filtering step(s)
#>   transmission pairs: 0 
```
