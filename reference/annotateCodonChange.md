# Codon and amino acid consequence of each SNV

Annotates every single-nucleotide substitution in a
[`WithinHostExperiment`](WithinHostExperiment.md) with the codon it
changes and the resulting amino acid, by passing the sites to
[`VariantAnnotation::predictCoding()`](https://rdrr.io/pkg/VariantAnnotation/man/predictCoding-methods.html).
The coding sequences are supplied either as a `TxDb` or as the CDS
records of an annotation in a `GRanges` or `GRangesList`; the reference
sequence is a `DNAStringSet`, as returned by
[`Biostrings::readDNAStringSet()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html).

## Usage

``` r
annotateCodonChange(x, cds, ref, ...)

# S4 method for class 'WithinHostExperiment,TxDb,DNAStringSet'
annotateCodonChange(x, cds, ref, collapse = c("longest", "all"), ...)

# S4 method for class 'WithinHostExperiment,GRangesList,DNAStringSet'
annotateCodonChange(
  x,
  cds,
  ref,
  geneCol = NULL,
  collapse = c("longest", "all"),
  ...
)

# S4 method for class 'WithinHostExperiment,GRanges,DNAStringSet'
annotateCodonChange(
  x,
  cds,
  ref,
  txCol = "ID",
  geneCol = "gene",
  collapse = c("longest", "all"),
  ...
)

# S4 method for class 'WithinHostExperiment,ANY,character'
annotateCodonChange(x, cds, ref, ...)

# S4 method for class 'WithinHostExperiment,ANY,XString'
annotateCodonChange(x, cds, ref, ...)

# S4 method for class 'WithinHostExperiment,character,DNAStringSet'
annotateCodonChange(x, cds, ref, ...)
```

## Arguments

- x:

  A [`WithinHostExperiment`](WithinHostExperiment.md) with `ref` and
  `alt` columns in `mcols(rowRanges(x))`.

- cds:

  The coding sequences: a `TxDb`; a `GRanges` of CDS records; or a
  `GRangesList` with one element per transcript, holding its CDS parts.

- ref:

  A `DNAStringSet` with the reference sequence of each seqlevel carrying
  coding sequences, for example from
  [`Biostrings::readDNAStringSet()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html).

- ...:

  Not used; supplying anything is an error.

- collapse:

  Either `"longest"` (default) or `"all"`; see Details.

- geneCol:

  The column of gene labels reported as `GENEID`. For a `GRanges` it is
  a column of the CDS records, `mcols(cds)`, defaulting to `"gene"`; for
  a `GRangesList` it is a column of the list's own metadata, one value
  per transcript, defaulting to `NULL`. `NULL` reports the transcript
  name instead of a gene label.

- txCol:

  For a `GRanges` `cds`: the column identifying the transcript each CDS
  record belongs to (default `"ID"`).

## Value

`x`, with these columns added to `mcols(rowRanges(x))` (`NA`, or an
empty element, for sites not annotated):

- TXNAME:

  Transcript the row describes.

- GENEID:

  Gene of that transcript.

- CONSEQUENCE:

  `"synonymous"`, `"nonsynonymous"` or `"nonsense"`; see Details. A site
  whose reference codon holds an ambiguous base is `NA` in every column.

- CDSLOC:

  Position of the site within the coding sequence.

- PROTEINLOC:

  Position of the affected codon in the protein.

- REFCODON, VARCODON:

  Reference and variant codon, on the transcript's strand.

- REFAA, VARAA:

  Reference and variant amino acid, `"*"` for a stop codon.

## Details

**Coding sequences.** `predictCoding()` reads coding sequences from a
`TxDb`. A `GRanges` of CDS records, for example the rows of
`rtracklayer::import()` with `type == "CDS"`, is grouped into
transcripts by the column `txCol`: in GFF3 files from NCBI and Ensembl
the CDS records of one protein share an `ID`. Each transcript's parts
are joined in order of transcription, so a CDS split across a ribosomal
frameshift, such as SARS-CoV-2 ORF1ab, is translated as one sequence.
The `TxDb` is then built with
[`txdbmaker::makeTxDb()`](https://rdrr.io/pkg/txdbmaker/man/makeTxDb.html),
which is why that route needs the txdbmaker package.

Building a `TxDb` directly from a GFF3 file with
[`txdbmaker::makeTxDbFromGFF()`](https://rdrr.io/pkg/txdbmaker/man/makeTxDbFromGFF.html)
is not equivalent for every annotation. For the NCBI SARS-CoV-2
annotation (NC_045512.2) it drops the ORF1ab gene, whose two CDS parts
overlap by one nucleotide, with a warning that its exon ranks could not
be inferred; sites in ORF1a and ORF1ab are then left unannotated. A
`TxDb` gives no record of what was dropped while it was built, so this
function cannot detect it. Pass the CDS records as a `GRanges` when in
doubt.

A GFF3 `phase` column, when present, is honoured: the first part of each
transcript is shortened by its phase, the number of bases before the
first complete codon, and the phase of every later part must agree with
the length of the parts before it. A later part whose phase disagrees
would shift the reading frame of the joined sequence, so the function
stops.

**Reference.** Sequence names in `ref` are cut at the first space, so a
FASTA header carrying a description still matches the annotation. The
function stops if the variant sites and the coding sequences share no
sequence name (see [`annotateGenes`](annotateGenes.md)), or if a
sequence carrying coding sequences is missing from `ref`. Sites whose
`ref` allele differs from the reference base at that position are left
unannotated, with a warning, since they were called against a different
sequence. The length of each reference sequence is recorded in
`seqlengths(rowRanges(x))`.

**Sites annotated.** Only single-nucleotide substitutions with alleles
A, C, G or T are annotated; insertions, deletions and multi-nucleotide
substitutions get `NA`. Alleles are read on the plus strand and
complemented by `predictCoding()` for transcripts on the minus strand.

**Classification.** `CONSEQUENCE` is taken from `predictCoding()`
unchanged: `"synonymous"` when the amino acid is kept, `"nonsense"` when
a stop codon is gained, and `"nonsynonymous"` for every other change,
including the loss of a stop codon. A codon holding an ambiguous base is
not translated; such a site is left unannotated, with a warning.

`predictCoding()` applies the initiator rule of
[`Biostrings::translate()`](https://rdrr.io/pkg/Biostrings/man/translate.html)
only when the substitution is the first base of a coding sequence: there
an alternative start codon (`CTG` or `TTG`) is reported as `M`, so `ATG`
to `CTG` at that one position is synonymous. A change at the second or
third base of the first codon is translated with the standard table,
giving `L` for `CTG`.

**Sites in several transcripts.** `predictCoding()` returns one row per
site and transcript, so a site where coding sequences overlap (ORF1a and
ORF1ab in SARS-CoV-2) has several. With `collapse = "longest"` (the
default) the row of the transcript with the longest coding sequence is
kept, ties going to the transcript that appears first in `cds`, and a
message reports how many sites had more than one transcript and at how
many of them the consequence differed between transcripts. With
`collapse = "all"` every row is kept: each column becomes a
`CharacterList` or `IntegerList` with one element per transcript.

A nucleotide that one transcript reads twice, the nucleotide shared by
the two parts of a CDS split across a -1 ribosomal frameshift
(SARS-CoV-2 ORF1ab at 13468), gives that transcript two rows, one for
each reading. With `collapse = "longest"` the first reading, at the
lower CDS position, is kept, and a message says how many sites this
affected and at how many the two readings differ in consequence;
`collapse = "all"` keeps both. [`pNpSWithinHost`](pNpSWithinHost.md)
counts such a site by the first of its readings in a sense codon.

## References

Obenchain V, Lawrence M, Carey V, Gogarten S, Shannon P, Morgan M
(2014). VariantAnnotation: a Bioconductor package for exploration and
annotation of genetic variants. *Bioinformatics* 30:2076-2078.
[doi:10.1093/bioinformatics/btu168](https://doi.org/10.1093/bioinformatics/btu168)

## See also

[`annotateGenes`](annotateGenes.md) for gene labels without codons;
[`VariantAnnotation::predictCoding()`](https://rdrr.io/pkg/VariantAnnotation/man/predictCoding-methods.html).

## Examples

``` r
## A 15-nucleotide coding sequence, ATG CTG AAA GGG TAA
ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
cds <- GenomicRanges::GRanges("seg1", IRanges::IRanges(1, 15),
    strand = "+", ID = "cds-A", gene = "geneA")
gr <- GenomicRanges::GRanges("seg1",
    IRanges::IRanges(c(1, 6, 13), width = 1))
S4Vectors::mcols(gr)$ref <- c("A", "G", "T")
S4Vectors::mcols(gr)$alt <- c("G", "A", "C")
whe <- WithinHostExperiment(
    assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
    rowRanges = gr,
    colData = S4Vectors::DataFrame(sample_id = "S1"))
if (requireNamespace("txdbmaker", quietly = TRUE)) {
    whe <- annotateCodonChange(whe, cds, ref)
    mc <- S4Vectors::mcols(SummarizedExperiment::rowRanges(whe))
    mc[, c("GENEID", "REFCODON", "VARCODON", "REFAA", "VARAA",
           "CONSEQUENCE")]
}
#> DataFrame with 3 rows and 6 columns
#>        GENEID    REFCODON    VARCODON       REFAA       VARAA   CONSEQUENCE
#>   <character> <character> <character> <character> <character>   <character>
#> 1       geneA         ATG         GTG           M           V nonsynonymous
#> 2       geneA         CTG         CTA           L           L    synonymous
#> 3       geneA         TAA         CAA           *           Q nonsynonymous
```
