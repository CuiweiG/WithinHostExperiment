test_that("shannonISNV computes correct value", {
    expect_equal(shannonISNV(0.5), log(2), tolerance = 1e-10)
    expect_equal(shannonISNV(numeric(0)), 0)
    expect_equal(shannonISNV(c(NA, NA)), 0)
    expect_equal(shannonISNV(c(0, 0.5, 1)), log(2), tolerance = 1e-10)
})

test_that("piISNV computes correct value", {
    expect_equal(piISNV(0.5, 1000), 0.0005, tolerance = 1e-10)
    expect_equal(piISNV(numeric(0), 1000), 0)
})

test_that("piISNV with finite-sample correction", {
    raw <- piISNV(0.5, 1000)
    corrected <- piISNV(0.5, 1000, meanDepth = 100)
    ## Correction factor = 100/99
    expect_equal(corrected, raw * 100 / 99, tolerance = 1e-10)
    ## Higher depth -> correction closer to 1
    deep <- piISNV(0.5, 1000, meanDepth = 10000)
    expect_true(abs(deep - raw) < abs(corrected - raw))
})

test_that("piISNV rejects invalid genomeLength", {
    expect_error(piISNV(0.5, -1))
    expect_error(piISNV(0.5, 0))
})

test_that("wattersonISNV returns positive value", {
    w <- wattersonISNV(15, 13588, 2000)
    expect_true(is.numeric(w))
    expect_true(w > 0)
})

test_that("wattersonISNV returns 0 for 0 sites", {
    expect_equal(wattersonISNV(0, 1000, 500), 0)
})

test_that("calcDiversity works on WHE", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    div <- calcDiversity(whe, genomeLength = 1000L)
    expect_s4_class(div, "DataFrame")
    expect_true("shannon" %in% colnames(div))
    expect_true("pi" %in% colnames(div))
    expect_true("richness" %in% colnames(div))
    expect_true(div$shannon > 0)
    expect_true(div$richness > 0)
})

test_that("calcDiversity requires genomeLength for pi", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "d"),
        caller = "ivar")
    expect_error(calcDiversity(whe, indices = "pi"), "genomeLength")
})

test_that("calcDiversity warns when totalDepth missing for watterson", {
    gr <- GenomicRanges::GRanges("s",
        IRanges::IRanges(1:3, width = 1))
    whe <- WithinHostExperiment(
        assays    = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
        rowRanges = gr,
        colData   = S4Vectors::DataFrame(sample_id = "S1"))
    expect_warning(
        calcDiversity(whe, indices = "watterson",
                      genomeLength = 1000L),
        "totalDepth"
    )
})

test_that("calcDiversity respects QC pass filter", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "d"),
        caller = "ivar")
    div_before <- calcDiversity(whe, indices = "richness")

    whe <- flagISNV(whe, ISNVFilter(minDepth = 500L))
    div_after <- calcDiversity(whe, indices = "richness")
    expect_true(div_after$richness <= div_before$richness)
})

# --- simpsonISNV tests ---

test_that("simpsonISNV computes correct value", {
    ## Single site at p=0.5: 1 - (0.25 + 0.25) = 0.5
    expect_equal(simpsonISNV(0.5), 0.5, tolerance = 1e-10)
    ## Empty input
    expect_equal(simpsonISNV(numeric(0)), 0)
    expect_equal(simpsonISNV(c(NA, NA)), 0)
    ## 0 and 1 are excluded, only 0.5 counted
    expect_equal(simpsonISNV(c(0, 0.5, 1)), 0.5, tolerance = 1e-10)
})

test_that("simpsonISNV is higher for more even single site", {
    ## p=0.5 is maximally diverse, p=0.1 is less diverse
    expect_true(simpsonISNV(0.5) > simpsonISNV(0.1))
})

# --- chao1ISNV tests ---

test_that("chao1ISNV computes correct value", {
    ## 2 singletons, 1 doubleton: 5 + 4/2 = 7
    expect_equal(chao1ISNV(c(1, 1, 2, 5, 10), detected = 5), 7)
})

test_that("chao1ISNV handles zero doubletons", {
    ## f1=2, f2=0: bias-corrected = 5 + 2*1/2 = 6
    expect_equal(chao1ISNV(c(1, 1, 3, 5, 10), detected = 5), 6)
})

test_that("chao1ISNV returns detected when no singletons", {
    expect_equal(chao1ISNV(c(3, 4, 5), detected = 3), 3)
})

test_that("chao1ISNV handles empty counts", {
    expect_equal(chao1ISNV(integer(0), detected = 5), 5)
})

# --- calcDiversity with simpson and chao1 ---

test_that("calcDiversity supports simpson index", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    div <- calcDiversity(whe, indices = "simpson")
    expect_s4_class(div, "DataFrame")
    expect_true("simpson" %in% colnames(div))
    expect_true(is.numeric(div$simpson))
})

test_that("calcDiversity supports chao1 index", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    div <- calcDiversity(whe, indices = "chao1")
    expect_s4_class(div, "DataFrame")
    expect_true("chao1" %in% colnames(div))
    expect_true(div$chao1 >= 0)
})

test_that("calcDiversity warns when altCount missing for chao1", {
    gr <- GenomicRanges::GRanges("s",
        IRanges::IRanges(1:3, width = 1))
    whe <- WithinHostExperiment(
        assays    = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
        rowRanges = gr,
        colData   = S4Vectors::DataFrame(sample_id = "S1"))
    expect_warning(
        calcDiversity(whe, indices = "chao1"),
        "altCount"
    )
})
