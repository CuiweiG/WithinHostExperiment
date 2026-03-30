library(S4Vectors)
library(SummarizedExperiment)
library(GenomicRanges)

# ---- as(whe, "data.frame") ----

test_that("coerce WHE to data.frame produces long format", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "donor"),
        caller = "ivar")
    df <- as(whe, "data.frame")
    expect_true(is.data.frame(df))
    expect_true(all(c("chrom", "position", "ref", "alt",
                      "sample_id", "altFreq") %in% colnames(df)))
    expect_equal(nrow(df), 15L)
    expect_true(all(df$altFreq >= 0 & df$altFreq <= 1))
})

test_that("coerce multi-sample WHE to data.frame", {
    vcf1 <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(vcf1, vcf2),
        colData = DataFrame(sample_id = c("D", "R")),
        caller = "ivar")
    df <- as(whe, "data.frame")
    expect_true(is.data.frame(df))
    expect_equal(nrow(df), 25L)
    expect_true(all(c("D", "R") %in% df$sample_id))
})

test_that("coerce WHE to data.frame includes qcPass", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "d"),
        caller = "ivar")
    whe <- flagISNV(whe, ISNVFilter(minDepth = 500L))
    df <- as(whe, "data.frame")
    expect_true("qcPass" %in% colnames(df))
    expect_true(any(!df$qcPass))
})

# ---- as(whe, "VRanges") ----

test_that("coerce WHE to VRanges", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = DataFrame(sample_id = "donor"),
        caller = "ivar")
    vr <- as(whe, "VRanges")
    expect_s4_class(vr, "VRanges")
    expect_equal(length(vr), 15L)
    expect_true(all(!is.na(VariantAnnotation::ref(vr))))
    expect_true(all(!is.na(VariantAnnotation::alt(vr))))
})

test_that("coerce multi-sample WHE to VRanges preserves sample names", {
    vcf1 <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(vcf1, vcf2),
        colData = DataFrame(sample_id = c("D", "R")),
        caller = "ivar")
    vr <- as(whe, "VRanges")
    expect_s4_class(vr, "VRanges")
    expect_true(all(c("D", "R") %in% sampleNames(vr)))
    expect_equal(length(vr), 25L)
})
