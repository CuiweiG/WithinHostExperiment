# Label variant sites with the genes they overlap

Adds a `GENE` column to `mcols(rowRanges(x))` giving, for each variant
site, the label of the gene range it falls in. Genes are supplied as a
`GRanges`, which is how `rtracklayer::import()` returns a GFF3 or GTF
file, or as a `TxDb`, from which
[`GenomicFeatures::genes()`](https://rdrr.io/pkg/GenomicFeatures/man/transcripts.html)
extracts them. The overlap itself is
[`GenomicRanges::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)
with strand ignored, since variant sites carry no strand.

## Usage

``` r
annotateGenes(x, genes, ...)

# S4 method for class 'WithinHostExperiment,GRanges'
annotateGenes(x, genes, nameCol = "Name", collapse = c("longest", "all"), ...)

# S4 method for class 'WithinHostExperiment,TxDb'
annotateGenes(
  x,
  genes,
  nameCol = "gene_id",
  collapse = c("longest", "all"),
  ...
)
```

## Arguments

- x:

  A [`WithinHostExperiment`](WithinHostExperiment.md).

- genes:

  A `GRanges` of gene ranges, or a `TxDb`.

- ...:

  Not used; supplying anything is an error.

- nameCol:

  Character scalar. The column of `mcols(genes)` holding the label to
  copy. Defaults to `"Name"`, the GFF3 attribute that
  `rtracklayer::import()` returns as a column; for a `TxDb` it defaults
  to `"gene_id"`. An unknown column is an error that lists the columns
  available.

- collapse:

  Either `"longest"` (default) or `"all"`; see Details.

## Value

`x` with a `GENE` column in `mcols(rowRanges(x))`: a character vector
for `collapse = "longest"`, `NA` for sites in no gene, or a
`CharacterList` for `collapse = "all"`, with an empty element for sites
in no gene.

## Details

Before overlapping, the sequence names in use on both sides are
compared. If they share none, the function stops and prints both sets:
the same reference is often named differently by different sources (for
SARS-CoV-2, `NC_045512.2` in RefSeq and `MN908947.3` in GenBank), and an
overlap across mismatched names would otherwise return a column that is
`NA` everywhere. If they share some but not all, a warning names the
variant sequences that have no gene range.

A site can fall in more than one gene where genes overlap.
`collapse = "longest"` keeps the label of the longest overlapping range
(the first in `genes` on a tie) and reports with a message how many
sites carried more than one distinct label; `collapse = "all"` keeps
every distinct label in a `CharacterList`, one element per site, in the
order of `genes`. Several ranges with the same label, such as the two
CDS parts of a gene split by a ribosomal frameshift, count as one gene.

Ranges whose label is missing are dropped with a warning; if none is
left, the function stops. For a `TxDb`,
[`GenomicFeatures::genes()`](https://rdrr.io/pkg/GenomicFeatures/man/transcripts.html)
leaves out, with a message, genes whose transcripts lie on more than one
strand or sequence, so sites in them are not labelled.

## See also

[`annotateCodonChange`](annotateCodonChange.md) for the codon and amino
acid consequence of each substitution.

## Examples

``` r
vcf <- system.file("extdata", "test_donor.vcf",
    package = "WithinHostExperiment")
whe <- readWithinHost(vcf,
    colData = S4Vectors::DataFrame(sample_id = "donor"),
    caller = "ivar")

## Gene ranges built directly
genes <- GenomicRanges::GRanges("test_segment",
    IRanges::IRanges(c(1, 501), c(500, 1000)), strand = "+",
    Name = c("geneA", "geneB"))
whe <- annotateGenes(whe, genes)
table(S4Vectors::mcols(SummarizedExperiment::rowRanges(whe))$GENE,
    useNA = "ifany")
#> 
#> geneA geneB 
#>     7     8 

## The same genes read from a GFF3 file
if (requireNamespace("rtracklayer", quietly = TRUE)) {
    gff <- tempfile(fileext = ".gff3")
    writeLines(c("##gff-version 3",
        "test_segment\t.\tgene\t1\t500\t.\t+\t.\tID=gA;Name=geneA",
        "test_segment\t.\tgene\t501\t1000\t.\t+\t.\tID=gB;Name=geneB"),
        gff)
    features <- rtracklayer::import(gff)
    whe <- annotateGenes(whe, features[features$type == "gene"])
}

## Sequence names that do not match are an error, not a column of NA
renamed <- GenomeInfoDb::renameSeqlevels(genes,
    c(test_segment = "segment_1"))
try(annotateGenes(whe, renamed))
#> Error : The variant sites and 'genes' share no sequence names, so no site could be annotated.
#>   variant seqlevels: test_segment
#>   'genes' seqlevels: segment_1
#> Rename one of them to match, for example with GenomeInfoDb::renameSeqlevels().
```
