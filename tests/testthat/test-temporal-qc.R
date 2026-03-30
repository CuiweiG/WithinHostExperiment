# tests/testthat/test-temporal-qc.R

library(GenomicRanges)
library(S4Vectors)

test_that("flagTemporalInconsistency flags transient iSNVs", {
    gr <- GRanges("seg", IRanges::IRanges(c(100, 200, 300), width = 1))
    mcols(gr)$ref <- c("A", "C", "G")
    mcols(gr)$alt <- c("T", "G", "A")
    freq <- matrix(c(0.1, NA,  0.2,
                     0.05, NA,  NA,
                     0.15, 0.2, 0.1),
                   nrow = 3, ncol = 3)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = DataFrame(
            sample_id = c("t1", "t2", "t3"),
            host_id = c("H1", "H1", "H1"),
            timepoint = c(1, 2, 3)))

    whe2 <- flagTemporalInconsistency(whe, minTimepoints = 2L)
    tc <- mcols(SummarizedExperiment::rowRanges(whe2))$temporal_class
    expect_equal(tc[1], "persistent")   # 2 timepoints
    expect_equal(tc[2], "transient")    # 1 timepoint
    expect_equal(tc[3], "persistent")   # 3 timepoints
    # transient should be flagged
    qc <- SummarizedExperiment::assay(whe2, "qcPass")
    expect_true(all(!qc[2, ]))  # row 2 all FALSE
})

test_that("flagTemporalInconsistency with minConsecutive", {
    gr <- GRanges("seg", IRanges::IRanges(c(100, 200), width = 1))
    mcols(gr)$ref <- c("A", "C")
    mcols(gr)$alt <- c("T", "G")
    # var1: days 1,3 (not consecutive), var2: days 1,2,3
    freq <- matrix(c(0.1,  0.15,
                     NA,   0.2,
                     0.12, 0.1),
                   nrow = 2, ncol = 3)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = DataFrame(
            sample_id = c("t1", "t2", "t3"),
            host_id = c("H1", "H1", "H1"),
            timepoint = c(1, 2, 3)))

    whe2 <- flagTemporalInconsistency(whe, minTimepoints = 2L,
                                       minConsecutive = 2L)
    tc <- mcols(SummarizedExperiment::rowRanges(whe2))$temporal_class
    expect_equal(tc[1], "transient")    # 2 detections but not consecutive
    expect_equal(tc[2], "persistent")   # 3 consecutive
})

test_that("flagTemporalInconsistency logs the step", {
    gr <- GRanges("seg", IRanges::IRanges(100, width = 1))
    freq <- matrix(c(0.1, NA, 0.1), nrow = 1, ncol = 3)
    whe <- WithinHostExperiment(
        assays = list(altFreq = freq),
        rowRanges = gr,
        colData = DataFrame(
            sample_id = c("t1", "t2", "t3"),
            host_id = c("H1", "H1", "H1"),
            timepoint = c(1, 2, 3)))
    whe2 <- flagTemporalInconsistency(whe)
    log <- qcLog(whe2)
    expect_true("temporal_consistency" %in% log$step)
})
