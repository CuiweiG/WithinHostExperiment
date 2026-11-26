# Consensus sequence of each sample

Builds, for every sample of a
[`WithinHostExperiment`](WithinHostExperiment.md), the reference
sequence with each substitution present at or above `minFreq` written
in, and returns the sequences as a `DNAStringSet`. Its purpose is
comparison between hosts: the within-host ratio of
[`pNpSWithinHost`](pNpSWithinHost.md) describes polymorphism inside one
host, whereas dN/dS between hosts is a divergence between their
consensus sequences, which the models of
[`MSA2dist::dnastring2kaks()`](https://rdrr.io/pkg/MSA2dist/man/dnastring2kaks.html)
estimate.

## Usage

``` r
consensusSequences(x, ref, ...)

# S4 method for class 'WithinHostExperiment,DNAStringSet'
consensusSequences(
  x,
  ref,
  minFreq = 0.5,
  usePassedOnly = TRUE,
  seqname = NULL,
  region = NULL,
  ...
)

# S4 method for class 'WithinHostExperiment,character'
consensusSequences(x, ref, ...)

# S4 method for class 'WithinHostExperiment,XString'
consensusSequences(x, ref, ...)
```

## Arguments

- x:

  A [`WithinHostExperiment`](WithinHostExperiment.md) with `ref` and
  `alt` columns in `mcols(rowRanges(x))`.

- ref:

  A `DNAStringSet` holding the reference sequence the variants were
  called against, for example from
  [`Biostrings::readDNAStringSet()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html).
  Sequence names are cut at the first space.

- ...:

  Not used; supplying anything is an error.

- minFreq:

  Numeric scalar in (0, 1\]. The frequency at or above which a
  substitution is written in (default 0.5).

- usePassedOnly:

  Logical. If `TRUE` (default), only sites with `qcPass == TRUE` in that
  sample are written in.

- seqname:

  The sequence of `ref` to build; required when `ref` holds more than
  one sequence and `region` does not name one.

- region:

  Optional `GRanges` on `seqname`, all on one strand: the part of the
  sequence to return; see Details.

## Value

A `DNAStringSet` with one sequence per sample, named by sample. Its
[`mcols()`](https://rdrr.io/pkg/S4Vectors/man/Vector-class.html) hold
`n_substitutions`, the number of substitutions written into the returned
part of each sequence.

## Details

Only substitutions are written in: single-nucleotide changes and
multi-nucleotide changes whose reference and alternative alleles have
the same length. Insertions and deletions are not applied, with a
warning, so every consensus has the length of the reference and stays
aligned to it and to every other sample's consensus; that is what
[`MSA2dist::dnastring2kaks()`](https://rdrr.io/pkg/MSA2dist/man/dnastring2kaks.html)
needs. Where a deletion or insertion is fixed in a host, the returned
sequence therefore differs from that host's true consensus.

A site enters a sample's consensus when its `altFreq` is at least
`minFreq` and, with `usePassedOnly = TRUE`, it passed QC. With `minFreq`
at or below 0.5 two alleles can qualify at one position; the more
frequent is written in, the first in the object on a tie, and a warning
reports how many were set aside. A site whose `ref` allele differs from
the reference at that position is never written in, with a warning,
since it was called against another sequence.

`region` restricts the output to part of the sequence, such as one
coding sequence: its ranges are taken in order along the sequence and
joined, and the result is reverse-complemented when they lie on the
minus strand, so a coding sequence split into parts comes back as one
sequence read in its own direction. Positions are those of the
reference, since no length-changing edit is made.

## See also

[`pNpSWithinHost`](pNpSWithinHost.md) for the within-host ratio;
[`MSA2dist::dnastring2kaks()`](https://rdrr.io/pkg/MSA2dist/man/dnastring2kaks.html)
for dN/dS between sequences.

## Examples

``` r
## A coding sequence of 15 codons and three hosts, each carrying one
## synonymous and one nonsynonymous change at high frequency:
## host1 6 G>A (CTG -> CTA) and 7 A>G (AAA -> GAA),
## host2 12 G>A (GGG -> GGA) and 13 T>C (TTT -> CTT),
## host3 18 C>A (CCC -> CCA) and 19 A>C (AAA -> CAA)
ref <- Biostrings::DNAStringSet(c(
    seg1 = "ATGCTGAAAGGGTTTCCCAAACTGATGGTTCTGAAAGGGTTTTAA"))
gr <- GenomicRanges::GRanges("seg1",
    IRanges::IRanges(c(6, 7, 12, 13, 18, 19), width = 1))
S4Vectors::mcols(gr)$ref <- c("G", "A", "G", "T", "C", "A")
S4Vectors::mcols(gr)$alt <- c("A", "G", "A", "C", "A", "C")
freqs <- cbind(host1 = c(0.95, 0.91, 0.02, 0.03, 0.01, 0.02),
               host2 = c(0.03, 0.02, 0.88, 0.93, 0.04, 0.01),
               host3 = c(0.01, 0.04, 0.02, 0.01, 0.97, 0.90))
whe <- WithinHostExperiment(
    assays = list(altFreq = freqs),
    rowRanges = gr,
    colData = S4Vectors::DataFrame(sample_id = colnames(freqs),
                                   row.names = colnames(freqs)))
cons <- consensusSequences(whe, ref)
cons
#> DNAStringSet object of length 3:
#>     width sequence                                          names               
#> [1]    45 ATGCTAGAAGGGTTTCCCAAACTGATGGTTCTGAAAGGGTTTTAA     host1
#> [2]    45 ATGCTGAAAGGACTTCCCAAACTGATGGTTCTGAAAGGGTTTTAA     host2
#> [3]    45 ATGCTGAAAGGGTTTCCACAACTGATGGTTCTGAAAGGGTTTTAA     host3

## dN/dS between the hosts' consensus coding sequences: every pair
## differs by two synonymous and two nonsynonymous changes
if (requireNamespace("MSA2dist", quietly = TRUE)) {
    kaks <- suppressMessages(MSA2dist::dnastring2kaks(cons, model = "Li"))
    kaks[, c("seq1", "seq2", "Ka", "Ks", "Ka/Ks")]
}
#>    seq1  seq2         Ka        Ks     Ka/Ks
#> 1 host1 host2 0.07421000 0.1786266 0.4154477
#> 2 host1 host3 0.06400214 0.2538504 0.2521254
#> 3 host2 host3 0.06436570 0.2384559 0.2699270
```
