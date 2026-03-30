# tests/testthat/test-import-multicallers.R
# Tests for LoFreq and Freebayes VCF format support

library(S4Vectors)
library(SummarizedExperiment)

# ---- LoFreq ----

test_that("readWithinHost parses LoFreq VCF correctly", {
    vcf <- system.file("extdata", "test_donor_lofreq.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"),
        caller = "lofreq")
    expect_s4_class(whe, "WithinHostExperiment")
    expect_equal(nrow(whe), 15L)
    freq <- assay(whe, "altFreq")[, 1]
    expect_true(all(freq >= 0 & freq <= 1, na.rm = TRUE))
})

test_that("LoFreq frequencies match iVar for same data", {
    vcf_lf <- system.file("extdata", "test_donor_lofreq.vcf",
                          package = "WithinHostExperiment")
    vcf_iv <- system.file("extdata", "test_donor.vcf",
                          package = "WithinHostExperiment")
    whe_lf <- readWithinHost(vcf_lf,
        colData = DataFrame(sample_id = "d"), caller = "lofreq")
    whe_iv <- readWithinHost(vcf_iv,
        colData = DataFrame(sample_id = "d"), caller = "ivar")
    freq_lf <- sort(assay(whe_lf, "altFreq")[, 1])
    freq_iv <- sort(assay(whe_iv, "altFreq")[, 1])
    expect_equal(freq_lf, freq_iv, tolerance = 0.01)
})

test_that("readWithinHost auto-detects LoFreq", {
    vcf <- system.file("extdata", "test_donor_lofreq.vcf",
                       package = "WithinHostExperiment")
    expect_message(
        whe <- readWithinHost(vcf,
            colData = DataFrame(sample_id = "d"),
            caller = "auto"),
        "lofreq"
    )
    expect_s4_class(whe, "WithinHostExperiment")
})

# ---- Freebayes ----

test_that("readWithinHost parses Freebayes VCF correctly", {
    vcf <- system.file("extdata", "test_donor_freebayes.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"),
        caller = "freebayes")
    expect_s4_class(whe, "WithinHostExperiment")
    expect_equal(nrow(whe), 15L)
    freq <- assay(whe, "altFreq")[, 1]
    expect_true(all(freq >= 0 & freq <= 1, na.rm = TRUE))
})

test_that("readWithinHost auto-detects Freebayes", {
    vcf <- system.file("extdata", "test_donor_freebayes.vcf",
                       package = "WithinHostExperiment")
    expect_message(
        whe <- readWithinHost(vcf,
            colData = DataFrame(sample_id = "d"),
            caller = "auto"),
        "freebayes"
    )
})

test_that("Freebayes depths are correct", {
    vcf <- system.file("extdata", "test_donor_freebayes.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"),
        caller = "freebayes")
    dp <- assay(whe, "totalDepth")[, 1]
    expect_true(all(dp > 0, na.rm = TRUE))
    expect_true(is.integer(dp))
})
