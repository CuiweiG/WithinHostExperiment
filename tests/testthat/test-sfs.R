# tests/testthat/test-sfs.R

library(GenomicRanges)
library(S4Vectors)

## Helper
.make_whe <- function(n = 10) {
    gr <- GRanges("seg", IRanges::IRanges(seq_len(n) * 100L, width = 1))
    freq <- matrix(runif(n, 0.03, 0.45), ncol = 1)
    depth <- matrix(as.integer(runif(n, 500, 3000)), ncol = 1)
    WithinHostExperiment(
        assays = list(altFreq = freq, totalDepth = depth),
        rowRanges = gr,
        colData = DataFrame(sample_id = "S1"))
}

# ---- buildSFS ----

test_that("buildSFS creates a valid WithinHostSFS", {
    whe <- .make_whe()
    sfs <- buildSFS(whe, genomeLength = 1000L)
    expect_s4_class(sfs, "WithinHostSFS")
    expect_true(validObject(sfs))
    expect_equal(sfs@sampleId, "S1")
    expect_equal(sfs@genomeLength, 1000L)
    expect_true(sfs@folded)
    expect_false(sfs@corrected)
    expect_true(sum(sfs@counts) > 0)
})

test_that("buildSFS fold=FALSE gives unfolded SFS", {
    whe <- .make_whe()
    sfs <- buildSFS(whe, fold = FALSE, genomeLength = 1000L)
    expect_false(sfs@folded)
    expect_true(max(sfs@breaks) > 0.5)
})

test_that("buildSFS respects threshold", {
    whe <- .make_whe()
    sfs <- buildSFS(whe, threshold = 0.10, genomeLength = 1000L)
    expect_equal(sfs@threshold, 0.10)
    expect_true(min(sfs@breaks) >= 0.10)
})

test_that("buildSFSList returns one SFS per sample", {
    data(example_whe)
    sfs_list <- buildSFSList(example_whe, genomeLength = 1000L)
    expect_equal(length(sfs_list), ncol(example_whe))
    expect_true(all(vapply(sfs_list, is, logical(1), "WithinHostSFS")))
})

# ---- correctSFSBias ----

test_that("correctSFSBias truncation method works", {
    whe <- .make_whe()
    sfs <- buildSFS(whe, genomeLength = 1000L)
    sfs_c <- correctSFSBias(sfs, method = "truncation")
    expect_true(sfs_c@corrected)
    expect_true(sum(sfs_c@counts) >= sum(sfs@counts))
})

test_that("correctSFSBias binomial method works", {
    whe <- .make_whe()
    sfs <- buildSFS(whe, genomeLength = 1000L)
    sfs_c <- correctSFSBias(sfs, method = "binomial")
    expect_true(sfs_c@corrected)
})

# ---- neutralityFromSFS ----

test_that("neutralityFromSFS returns all statistics", {
    whe <- .make_whe(20)
    sfs <- buildSFS(whe, fold = FALSE, genomeLength = 1000L)
    res <- neutralityFromSFS(sfs)
    expect_true(is.list(res))
    expect_true(all(c("S", "n", "pi", "theta_W", "theta_H",
                       "tajimaD", "fusFs", "fayWuH") %in% names(res)))
    expect_true(is.numeric(res$tajimaD))
    expect_true(is.numeric(res$fayWuH))
})

test_that("neutralityFromSFS handles empty SFS", {
    sfs <- new("WithinHostSFS",
        counts = integer(20), breaks = seq(0.03, 0.5, length.out = 21),
        folded = TRUE, sampleId = "empty", genomeLength = 1000L,
        nSites = 0L, meanDepth = 1000, threshold = 0.03,
        corrected = FALSE)
    res <- neutralityFromSFS(sfs)
    expect_equal(res$S, 0L)
    expect_equal(res$tajimaD, 0)
})

# ---- compareSFS ----

test_that("compareSFS compares two SFS objects", {
    whe <- .make_whe(20)
    sfs1 <- buildSFS(whe, genomeLength = 1000L)
    sfs2 <- buildSFS(whe, genomeLength = 1000L)
    res <- compareSFS(sfs1, sfs2)
    expect_true(is.list(res))
    expect_true("chisq_p" %in% names(res))
})

.sfs_with_counts <- function(counts) {
    new("WithinHostSFS",
        counts = as.integer(counts),
        breaks = seq(0.03, 0.5, length.out = length(counts) + 1L),
        folded = TRUE, sampleId = "s", genomeLength = 1000L,
        nSites = as.integer(sum(counts)), meanDepth = 1000,
        threshold = 0.03, corrected = FALSE)
}

test_that("compareSFS reports small expected counts instead of warning", {
    sfs1 <- .sfs_with_counts(c(3, 1, 0, 2))
    sfs2 <- .sfs_with_counts(c(1, 2, 1, 0))
    expect_no_warning(res <- compareSFS(sfs1, sfs2))
    mat <- rbind(c(3, 1, 0, 2), c(1, 2, 1, 0))
    reference <- suppressWarnings(stats::chisq.test(mat))
    expect_equal(res$chisq_stat, unname(reference$statistic))
    expect_equal(res$chisq_p, reference$p.value)
    expected <- outer(rowSums(mat), colSums(mat)) / sum(mat)
    expect_equal(res$min_expected, min(expected))
    expect_false(res$chisq_approx_ok)
})

test_that("compareSFS flags the approximation as reliable for large counts", {
    sfs1 <- .sfs_with_counts(c(40, 30, 20, 10))
    sfs2 <- .sfs_with_counts(c(35, 35, 15, 15))
    res <- compareSFS(sfs1, sfs2)
    expect_gte(res$min_expected, 5)
    expect_true(res$chisq_approx_ok)
})

test_that("compareSFS returns NA diagnostics when fewer than two bins remain", {
    sfs1 <- .sfs_with_counts(c(4, 0, 0))
    sfs2 <- .sfs_with_counts(c(2, 0, 0))
    res <- compareSFS(sfs1, sfs2)
    expect_true(is.na(res$chisq_p))
    expect_true(is.na(res$min_expected))
    expect_true(is.na(res$chisq_approx_ok))
})

test_that("compareSFS rejects mismatched bins", {
    whe <- .make_whe()
    sfs1 <- buildSFS(whe, nBins = 10, genomeLength = 1000L)
    sfs2 <- buildSFS(whe, nBins = 20, genomeLength = 1000L)
    expect_error(compareSFS(sfs1, sfs2), "same number of bins")
})

# ---- show ----

test_that("show method prints without error", {
    whe <- .make_whe()
    sfs <- buildSFS(whe, genomeLength = 1000L)
    expect_output(show(sfs), "WithinHostSFS")
})
