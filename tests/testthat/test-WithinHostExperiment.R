# tests/testthat/test-WithinHostExperiment.R

library(GenomicRanges)
library(SummarizedExperiment)
library(S4Vectors)

## Helper: build a minimal WHE for testing
.make_test_whe <- function(n_sites = 5, n_samples = 2, with_pairs = FALSE) {
    gr <- GRanges("seg", IRanges::IRanges(seq_len(n_sites) * 100L, width = 1))
    freq_mat  <- matrix(runif(n_sites * n_samples, 0.03, 0.45),
                        nrow = n_sites, ncol = n_samples)
    depth_mat <- matrix(as.integer(runif(n_sites * n_samples, 500, 3000)),
                        nrow = n_sites, ncol = n_samples)
    cd <- if (with_pairs) {
        DataFrame(
            sample_id = c("donor", "recipient"),
            host_id   = c("P1", "P2"),
            role      = c("donor", "recipient"),
            pair_id   = c("pair_1", "pair_1")
        )
    } else {
        DataFrame(sample_id = paste0("S", seq_len(n_samples)))
    }
    WithinHostExperiment(
        assays    = list(altFreq = freq_mat, totalDepth = depth_mat),
        rowRanges = gr,
        colData   = cd
    )
}

# ---- Construction ----

test_that("WithinHostExperiment can be constructed from components", {
    whe <- .make_test_whe()
    expect_s4_class(whe, "WithinHostExperiment")
    expect_true(is(whe, "RangedSummarizedExperiment"))
    expect_equal(nrow(whe), 5L)
    expect_equal(ncol(whe), 2L)
})

test_that("qcPass assay is auto-created as all TRUE", {
    whe <- .make_test_whe()
    expect_true("qcPass" %in% assayNames(whe))
    expect_true(all(assay(whe, "qcPass")))
})

test_that("WithinHostExperiment rejects missing sample_id", {
    gr <- GRanges("seg", IRanges::IRanges(1, width = 1))
    expect_error(
        WithinHostExperiment(
            assays    = list(altFreq = matrix(0.1, ncol = 1)),
            rowRanges = gr,
            colData   = DataFrame(name = "X")
        ),
        "sample_id"
    )
})

test_that("colnames are set from sample_id", {
    whe <- .make_test_whe()
    expect_equal(colnames(assay(whe, "altFreq")), c("S1", "S2"))
})

# ---- Subsetting (inherited) ----

test_that("row subsetting preserves class", {
    whe <- .make_test_whe()
    sub <- whe[1:2, ]
    expect_s4_class(sub, "WithinHostExperiment")
    expect_equal(nrow(sub), 2L)
    expect_equal(ncol(sub), 2L)
})

test_that("column subsetting preserves class", {
    whe <- .make_test_whe()
    sub <- whe[, 1]
    expect_s4_class(sub, "WithinHostExperiment")
    expect_equal(ncol(sub), 1L)
})

# ---- Accessors ----

test_that("qcLog starts empty", {
    whe <- .make_test_whe()
    log <- qcLog(whe)
    expect_s4_class(log, "DataFrame")
    expect_equal(nrow(log), 0L)
})

test_that("qcFilters starts empty", {
    whe <- .make_test_whe()
    expect_equal(length(qcFilters(whe)), 0L)
})

test_that("passedISNV returns all rows when qcPass is all TRUE", {
    whe    <- .make_test_whe()
    passed <- passedISNV(whe)
    expect_equal(nrow(passed), nrow(whe))
})

test_that("passedISNV filters rows with all-FALSE qcPass", {
    whe    <- .make_test_whe(n_sites = 3, n_samples = 1)
    qc_mat <- assay(whe, "qcPass")
    qc_mat[2, ] <- FALSE
    assay(whe, "qcPass") <- qc_mat
    passed <- passedISNV(whe)
    expect_equal(nrow(passed), 2L)
})

# ---- transmissionPairs ----

test_that("transmissionPairs extracts donor-recipient pairs", {
    whe   <- .make_test_whe(n_sites = 3, n_samples = 2, with_pairs = TRUE)
    pairs <- transmissionPairs(whe)
    expect_s4_class(pairs, "DataFrame")
    expect_equal(nrow(pairs), 1L)
    expect_equal(as.character(pairs$donor),     "donor")
    expect_equal(as.character(pairs$recipient), "recipient")
    expect_equal(as.character(pairs$pair_id),   "pair_1")
})

test_that("transmissionPairs returns empty when no role column", {
    whe   <- .make_test_whe()
    pairs <- transmissionPairs(whe)
    expect_equal(nrow(pairs), 0L)
})

test_that("transmissionPairs handles missing pair_id gracefully", {
    gr  <- GRanges("seg", IRanges::IRanges(1:2, width = 1))
    whe <- WithinHostExperiment(
        assays    = list(altFreq = matrix(c(0.1, 0.2), ncol = 1)),
        rowRanges = gr,
        colData   = DataFrame(sample_id = "X", role = "donor",
                              pair_id = NA_character_)
    )
    pairs <- transmissionPairs(whe)
    expect_equal(nrow(pairs), 0L)
})

# ---- show ----

test_that("show method produces expected output", {
    whe <- .make_test_whe()
    expect_output(show(whe),
                  "WithinHostExperiment with 5 variant sites and 2 samples")
    expect_output(show(whe), "assays")
    expect_output(show(whe), "qcLog: 0")
    expect_output(show(whe), "transmission pairs: 0")
})

test_that("show reports transmission pairs when present", {
    whe <- .make_test_whe(n_sites = 3, n_samples = 2, with_pairs = TRUE)
    expect_output(show(whe), "transmission pairs: 1")
})

# ---- .addQcStep ----

test_that(".addQcStep appends to qcLog", {
    whe  <- .make_test_whe()
    whe2 <- WithinHostExperiment:::.addQcStep(
        whe, "depth_filter", "minDepth", "100", 2L, 3L)
    expect_equal(nrow(qcLog(whe2)), 1L)
    expect_equal(qcLog(whe2)$step,      "depth_filter")
    expect_equal(qcLog(whe2)$n_flagged, 2L)

    whe3 <- WithinHostExperiment:::.addQcStep(
        whe2, "freq_filter", "minFreq", "0.03", 1L, 2L)
    expect_equal(nrow(qcLog(whe3)), 2L)
})

# ---- Subsetting preserves custom slots ----

test_that("subsetting preserves qcLog and qcFilters", {
    whe <- .make_test_whe(n_sites = 5, n_samples = 2)
    whe <- flagISNV(whe, ISNVFilter())

    ## Row subset
    sub_r <- whe[1:3, ]
    expect_s4_class(sub_r, "WithinHostExperiment")
    expect_equal(nrow(qcLog(sub_r)), nrow(qcLog(whe)))
    expect_equal(length(qcFilters(sub_r)), length(qcFilters(whe)))
    expect_equal(nrow(sub_r), 3L)

    ## Column subset
    sub_c <- whe[, 1]
    expect_s4_class(sub_c, "WithinHostExperiment")
    expect_equal(nrow(qcLog(sub_c)), nrow(qcLog(whe)))
    expect_equal(ncol(sub_c), 1L)

    ## Combined
    sub_rc <- whe[1:2, 1]
    expect_s4_class(sub_rc, "WithinHostExperiment")
    expect_equal(nrow(sub_rc), 2L)
    expect_equal(ncol(sub_rc), 1L)
    expect_equal(nrow(qcLog(sub_rc)), nrow(qcLog(whe)))
    expect_equal(length(qcFilters(sub_rc)), length(qcFilters(whe)))
})

# ---- Validity ----

test_that("WithinHostExperiment validity rejects invalid qcPass type", {
    gr <- GRanges("seg", IRanges::IRanges(1:3, width = 1))
    whe <- WithinHostExperiment(
        assays    = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
        rowRanges = gr,
        colData   = DataFrame(sample_id = "S1"))
    ## Corrupt qcPass to numeric (preserve dimnames)
    bad_mat <- matrix(1.0, nrow = 3, ncol = 1,
                      dimnames = list(NULL, "S1"))
    assay(whe, "qcPass", withDimnames = FALSE) <- bad_mat
    expect_error(validObject(whe), "qcPass")
})

test_that("WithinHostExperiment validity rejects non-ISNVFilter in qcFilters", {
    gr <- GRanges("seg", IRanges::IRanges(1:3, width = 1))
    whe <- WithinHostExperiment(
        assays    = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
        rowRanges = gr,
        colData   = DataFrame(sample_id = "S1"))
    whe@qcFilters <- list("not_a_filter")
    expect_error(validObject(whe), "ISNVFilter")
})

# ---- Replacement methods ----

test_that("qcLog<- and qcFilters<- replacement methods work", {
    whe <- .make_test_whe(n_sites = 3, n_samples = 1)

    new_log <- DataFrame(
        step = "test", parameter = "x", value = "1",
        n_flagged = 0L, n_passed = 3L, timestamp = "now")
    qcLog(whe) <- new_log
    expect_equal(nrow(qcLog(whe)), 1L)
    expect_equal(qcLog(whe)$step, "test")

    qcFilters(whe) <- list(ISNVFilter())
    expect_equal(length(qcFilters(whe)), 1L)
    expect_s4_class(qcFilters(whe)[[1]], "ISNVFilter")
})
