# Within-host pN/pS ratio

Computes, for each sample, the ratio of nonsynonymous to synonymous
polymorphism per site among its iSNVs. This is a within-host
polymorphism ratio, pN/pS: it counts variants segregating inside one
host, not substitutions fixed between two sequences, so it is not the
dN/dS of a between-sequence comparison and carries no multiple-hit
correction. Sites are counted with the method of Nei and Gojobori
(1986). For dN/dS between hosts or between consensus sequences, use the
models in the MSA2dist package.

## Usage

``` r
pNpSWithinHost(x, cds, ref, ...)

# S4 method for class 'WithinHostExperiment,TxDb,DNAStringSet'
pNpSWithinHost(x, cds, ref, usePassedOnly = TRUE, ...)

# S4 method for class 'WithinHostExperiment,GRangesList,DNAStringSet'
pNpSWithinHost(x, cds, ref, usePassedOnly = TRUE, geneCol = NULL, ...)

# S4 method for class 'WithinHostExperiment,GRanges,DNAStringSet'
pNpSWithinHost(
  x,
  cds,
  ref,
  usePassedOnly = TRUE,
  txCol = "ID",
  geneCol = "gene",
  ...
)

# S4 method for class 'WithinHostExperiment,ANY,character'
pNpSWithinHost(x, cds, ref, ...)

# S4 method for class 'WithinHostExperiment,ANY,XString'
pNpSWithinHost(x, cds, ref, ...)

# S4 method for class 'WithinHostExperiment,character,DNAStringSet'
pNpSWithinHost(x, cds, ref, ...)
```

## Arguments

- x:

  A [`WithinHostExperiment`](WithinHostExperiment.md) with `ref` and
  `alt` columns in `mcols(rowRanges(x))`.

- cds:

  The coding sequences, as for
  [`annotateCodonChange`](annotateCodonChange.md): a `TxDb`, a `GRanges`
  of CDS records, or a `GRangesList` with one element per transcript.

- ref:

  A `DNAStringSet` with the reference sequence of each seqlevel carrying
  coding sequences, for example from
  [`Biostrings::readDNAStringSet()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html).

- ...:

  Not used; supplying anything is an error.

- usePassedOnly:

  Logical. If `TRUE` (default), only iSNVs with `qcPass == TRUE` are
  counted.

- txCol, geneCol:

  As for [`annotateCodonChange`](annotateCodonChange.md): the transcript
  and gene columns of a `GRanges` or `GRangesList` `cds`.

## Value

A `DataFrame` with one row per sample and gene carrying at least one
counted iSNV (one row with `gene = NA` when a sample has none):

- sample_id:

  Sample identifier.

- nS, nN:

  Synonymous and nonsynonymous iSNVs in the sample.

- S_sites, N_sites:

  Synonymous and nonsynonymous sites over all coding sequences.

- pS, pN:

  The proportions \\n_S/S\\ and \\n_N/N\\.

- pNpS:

  \\p_N/p_S\\; `NA` when \\n_S = 0\\.

- gene, gene_nS, gene_nN, gene_S_sites, gene_N_sites, gene_pNpS:

  The same quantities for each gene.

## Details

With \\n_S\\ synonymous and \\n_N\\ nonsynonymous iSNVs present in a
sample, and \\S\\ synonymous and \\N\\ nonsynonymous sites in the coding
sequences, \$\$p_S = n_S / S, \quad p_N = n_N / N, \quad \mathrm{pN/pS}
= p_N / p_S.\$\$ Under the assumptions of the counting model (every
substitution equally likely, no selection) the expected ratio is 1. A
ratio below 1 is what purifying selection against amino acid changes
would give, but a within-host ratio rests on few iSNVs, is highly
variable, and should always be read with its counts. It is `NA` when no
synonymous iSNV is present.

**Classification.** The counts use only the codons and the genetic code
that the sites are counted over. This function calls
[`annotateCodonChange`](annotateCodonChange.md) on `x` with the same
`cds` and `ref`, so annotation columns already present in `x` are not
used, and takes from it, for each transcript reading an iSNV, the
reference and variant codon. A reading counts only when its reference
codon is a sense codon: an iSNV in a stop codon, including the loss of a
stop codon, is counted in neither \\n_S\\ nor \\n_N\\, since a stop
codon carries no sites. The reading is then classified by translating
both codons with the standard genetic code: synonymous when the amino
acid is kept, counting towards \\n_S\\, and nonsynonymous otherwise, a
change to a stop codon included, counting towards \\n_N\\. The initiator
rule is not applied, so `ATG` to `CTG` or `TTG` at the first base of a
coding sequence is nonsynonymous here, as `ATG` carries no synonymous
site, although `annotateCodonChange` reports it as synonymous.

The totals take each iSNV by its first counted reading, in the order of
`collapse = "longest"`: the longest transcript that reads it in a sense
codon, and where that transcript reads the nucleotide twice, the first
of those readings in a sense codon. Each gene's counts take it by the
first counted reading of that gene's own transcripts, so an iSNV in a
codon two genes share counts in both genes, as the codon's sites do, but
once in the total. Sites outside every coding sequence, insertions and
deletions, and codons holding an ambiguous base are not counted. An iSNV
is present in a sample when its frequency is above zero and, with
`usePassedOnly = TRUE`, it passed QC; it is counted once whatever its
frequency.

**Sites.** \\S\\ and \\N\\ are summed over the distinct codons of the
genome that some coding sequence reads: a codon is identified by the
genomic positions of its three bases, so where two coding sequences read
the same bases in the same frame (ORF1a inside ORF1ab in SARS-CoV-2) the
shared codons are counted once, while coding sequences that overlap in
different frames (ORF7a and ORF7b) contribute their codons separately.
Each sense codon carries three sites, split between synonymous and
nonsynonymous under the standard genetic code; stop codons and codons
with an ambiguous base carry none, so \\S + N\\ is three times the
number of distinct sense codons. When coding sequences of different
transcripts overlap, a warning names them with the nucleotides and
codons they share. The sites of a gene are summed in the same way over
the codons of that gene's own transcripts, so a codon shared by two
genes is counted in each of them but once in the total.

**Relation to MSA2dist.**
[`MSA2dist::codonmat2pnps()`](https://rdrr.io/pkg/MSA2dist/man/codonmat2pnps.html)
and
[`MSA2dist::dnastring2kaks()`](https://rdrr.io/pkg/MSA2dist/man/dnastring2kaks.html)
compare pairs of aligned sequences codon by codon, so they cannot take a
set of per-site iSNV calls, each present in a sample at its own
frequency; that is why the sites are counted here. Codon for codon the
counts agree: each of the 61 sense codons carries the same Nei-Gojobori
synonymous and nonsynonymous sites here as in
[`MSA2dist::codonmat2pnps()`](https://rdrr.io/pkg/MSA2dist/man/codonmat2pnps.html).
They differ only at stop codons, which carry no sites here but a third
or two thirds of a synonymous site there, so over a coding sequence
holding stop codons the two totals differ.

## References

Nei M, Gojobori T (1986). Simple methods for estimating the numbers of
synonymous and nonsynonymous nucleotide substitutions. *Mol Biol Evol*
3:418-426.
[doi:10.1093/oxfordjournals.molbev.a040410](https://doi.org/10.1093/oxfordjournals.molbev.a040410)

## See also

[`annotateCodonChange`](annotateCodonChange.md) for the reference and
variant codons of each iSNV, from which the classification here is made.

## Examples

``` r
## A 15-nucleotide coding sequence, ATG CTG AAA GGG TAA
ref <- Biostrings::DNAStringSet(c(seg1 = "ATGCTGAAAGGGTAA"))
cds <- GenomicRanges::GRanges("seg1", IRanges::IRanges(1, 15),
    strand = "+", ID = "cds-A", gene = "geneA")
## 1 A>G is ATG -> GTG (nonsynonymous); 6 G>A is CTG -> CTA (synonymous)
gr <- GenomicRanges::GRanges("seg1",
    IRanges::IRanges(c(1, 6), width = 1))
S4Vectors::mcols(gr)$ref <- c("A", "G")
S4Vectors::mcols(gr)$alt <- c("G", "A")
whe <- WithinHostExperiment(
    assays = list(altFreq = matrix(c(0.1, 0.2), ncol = 1)),
    rowRanges = gr,
    colData = S4Vectors::DataFrame(sample_id = "S1"))
if (requireNamespace("txdbmaker", quietly = TRUE)) {
    pNpSWithinHost(whe, cds, ref)
}
#> DataFrame with 1 row and 14 columns
#>     sample_id        nS        nN   S_sites   N_sites        pS        pN
#>   <character> <integer> <integer> <numeric> <numeric> <numeric> <numeric>
#> 1          S1         1         1   2.66667   9.33333     0.375  0.107143
#>        pNpS        gene   gene_nS   gene_nN gene_S_sites gene_N_sites gene_pNpS
#>   <numeric> <character> <integer> <integer>    <numeric>    <numeric> <numeric>
#> 1  0.285714       geneA         1         1      2.66667      9.33333  0.285714
```
