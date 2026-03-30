## Helper: build a WHE with amino acid annotations
.make_aa_whe <- function() {
    gr <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200, 300, 400, 500), width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "C", "G", "T", "A")
    S4Vectors::mcols(gr)$alt <- c("T", "G", "A", "C", "G")
    S4Vectors::mcols(gr)$REF_AA <- c("M", "L", "A", "V", "D")
    S4Vectors::mcols(gr)$ALT_AA <- c("I", "L", "V", "V", "E")
    freq <- matrix(c(0.10, 0.20, 0.15, 0.30, 0.05), ncol = 1)
    WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
}

test_that("dndsWithinHost classifies syn/nonsyn correctly", {
    whe <- .make_aa_whe()
    result <- dndsWithinHost(whe)
    expect_s4_class(result, "DataFrame")
    ## Sites: M->I (nonsyn), L->L (syn), A->V (nonsyn),
    ##        V->V (syn), D->E (nonsyn)
    expect_equal(result$nS[1], 2L)
    expect_equal(result$nN[1], 3L)
})

test_that("dndsWithinHost computes correct dN/dS ratio", {
    whe <- .make_aa_whe()
    result <- dndsWithinHost(whe)
    expect_equal(result$dNdS[1], 3 / 2, tolerance = 1e-10)
})

test_that("dndsWithinHost returns NA when no synonymous sites", {
    gr <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200), width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "C")
    S4Vectors::mcols(gr)$alt <- c("T", "G")
    S4Vectors::mcols(gr)$REF_AA <- c("M", "A")
    S4Vectors::mcols(gr)$ALT_AA <- c("I", "V")
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.1, 0.2), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    result <- dndsWithinHost(whe)
    expect_true(is.na(result$dNdS[1]))
})

test_that("dndsWithinHost errors with missing AA columns", {
    gr <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(100, width = 1))
    S4Vectors::mcols(gr)$ref <- "A"
    S4Vectors::mcols(gr)$alt <- "T"
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(0.1, ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    expect_error(dndsWithinHost(whe), "REF_AA")
})

test_that("dndsWithinHost includes per-gene breakdown", {
    gr <- GenomicRanges::GRanges("seg1",
        IRanges::IRanges(c(100, 200, 300), width = 1))
    S4Vectors::mcols(gr)$ref <- c("A", "C", "G")
    S4Vectors::mcols(gr)$alt <- c("T", "G", "A")
    S4Vectors::mcols(gr)$REF_AA <- c("M", "L", "A")
    S4Vectors::mcols(gr)$ALT_AA <- c("I", "L", "V")
    S4Vectors::mcols(gr)$GFF_FEATURE <- c("ORF1a", "ORF1a", "S")
    whe <- WithinHostExperiment(
        assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
        rowRanges = gr,
        colData = S4Vectors::DataFrame(sample_id = "S1"))
    result <- dndsWithinHost(whe)
    expect_true("gene" %in% colnames(result))
    expect_true("gene_nS" %in% colnames(result))
    expect_true("gene_dNdS" %in% colnames(result))
})

test_that("estimateSelectionCoefficient basic computation", {
    ## logit(0.5) = 0, logit(~0.731) = 1
    p1 <- 0.5
    p2 <- exp(1) / (1 + exp(1))  # logit(p2) = 1
    s <- estimateSelectionCoefficient(c(p1, p2))
    expect_equal(s, 1.0, tolerance = 1e-10)
})

test_that("estimateSelectionCoefficient with custom times", {
    freqs <- c(0.1, 0.3, 0.5)
    s <- estimateSelectionCoefficient(freqs, times = c(0, 5, 10))
    expect_length(s, 2L)
    expect_true(all(is.finite(s)))
})

test_that("estimateSelectionCoefficient errors on invalid input", {
    expect_error(estimateSelectionCoefficient(0.5), "at least 2")
    expect_error(estimateSelectionCoefficient(c(0, 0.5)),
                 "strictly between 0 and 1")
    expect_error(estimateSelectionCoefficient(c(0.5, NA)),
                 "NA")
})

test_that("estimateSelectionCoefficient respects generation_time", {
    freqs <- c(0.2, 0.4)
    s1 <- estimateSelectionCoefficient(freqs, generation_time = 1)
    s2 <- estimateSelectionCoefficient(freqs, generation_time = 2)
    expect_equal(s2, s1 * 2, tolerance = 1e-10)
})
