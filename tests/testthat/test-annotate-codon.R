# tests/testthat/test-annotate-codon.R

test_that("annotateCodonChange adds AA annotation from GFF+FASTA", {
    ## Build a minimal CDS GFF3 over the test_ref.fa sequence
    gff_file <- tempfile(fileext = ".gff3")
    writeLines(c(
        "##gff-version 3",
        "test_segment\t.\tCDS\t1\t999\t.\t+\t0\tgene=testGene"
    ), gff_file)

    ref <- system.file("extdata", "test_ref.fa",
        package = "WithinHostExperiment")
    vcf <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")

    whe2 <- annotateCodonChange(whe, gff_file, ref)
    mc <- S4Vectors::mcols(SummarizedExperiment::rowRanges(whe2))

    expect_true("REF_AA" %in% colnames(mc))
    expect_true("ALT_AA" %in% colnames(mc))
    expect_true("AA_CLASS" %in% colnames(mc))
    expect_true("GFF_FEATURE" %in% colnames(mc))
    expect_true("REF_CODON" %in% colnames(mc))
    expect_true("ALT_CODON" %in% colnames(mc))

    ## At least some variants should be annotated
    annotated <- !is.na(mc$REF_AA)
    expect_true(sum(annotated) > 0)

    ## AA_CLASS should be Synonymous, Nonsynonymous, or Nonsense
    classes <- mc$AA_CLASS[annotated]
    expect_true(all(classes %in% c("Synonymous", "Nonsynonymous",
                                    "Nonsense")))

    unlink(gff_file)
})

test_that("annotateCodonChange rejects missing ref/alt columns", {
    gr <- GenomicRanges::GRanges("seg",
        IRanges::IRanges(1:3, width = 1))
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(runif(3), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    gff_file <- tempfile(fileext = ".gff3")
    writeLines(c("##gff-version 3",
        "seg\t.\tCDS\t1\t100\t.\t+\t0\tgene=g"), gff_file)
    ref <- system.file("extdata", "test_ref.fa",
        package = "WithinHostExperiment")
    expect_error(annotateCodonChange(whe, gff_file, ref),
                 "ref.*alt")
    unlink(gff_file)
})

test_that("codon table translates correctly", {
    expect_equal(WithinHostExperiment:::.translate_codon("ATG"), "M")
    expect_equal(WithinHostExperiment:::.translate_codon("TAA"), "*")
    expect_equal(WithinHostExperiment:::.translate_codon("GCT"), "A")
    expect_true(is.na(WithinHostExperiment:::.translate_codon("XY")))
})

test_that("annotateCodonChange honours the GFF3 phase", {
    ## One leading base before the first complete codon (phase 1)
    ref <- tempfile(fileext = ".fa")
    writeLines(c(">seg1", "CATGCTGAAA"), ref)
    gff <- tempfile(fileext = ".gff3")
    writeLines(c("##gff-version 3",
        "seg1\t.\tCDS\t1\t10\t.\t+\t1\tgene=geneA"), gff)
    gr <- GenomicRanges::GRanges("seg1", IRanges::IRanges(c(2, 7, 1), width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "G", "C")
    S4Vectors::mcols(gr)$alt <- c("G", "A", "T")
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    mc <- S4Vectors::mcols(SummarizedExperiment::rowRanges(
        annotateCodonChange(whe, gff, ref)))
    expect_equal(mc$REF_CODON[1], "ATG")
    expect_equal(mc$ALT_CODON[1], "GTG")
    expect_equal(mc$REF_CODON[2], "CTG")
    expect_equal(mc$AA_CLASS[2], "Synonymous")
    expect_true(is.na(mc$REF_CODON[3]))
})

test_that("FASTA reader ignores whitespace and empty final records", {
    fa <- tempfile(fileext = ".fa")
    writeLines(c(">seg1 description", "ATG CTG ", "AAA\t", "", ">empty"), fa)
    seqs <- WithinHostExperiment:::.read_fasta_simple(fa)
    expect_identical(seqs[["seg1"]], "ATGCTGAAA")
    expect_identical(seqs[["empty"]], "")
})
