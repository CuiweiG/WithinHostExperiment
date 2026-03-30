# tests/testthat/test-adaptive-filter.R

test_that("adaptiveFilter('illumina') returns correct defaults", {
    filt <- adaptiveFilter("illumina")
    expect_s4_class(filt, "ISNVFilter")
    expect_equal(slot(filt, "minFreq"), 0.03)
    expect_equal(slot(filt, "minDepth"), 100L)
    # Non-specified slots keep ISNVFilter defaults
    expect_equal(slot(filt, "maxStrandBias"), 10.0)
    expect_equal(slot(filt, "minAltReads"), 10L)
})

test_that("adaptiveFilter('ont') sets high minFreq and disables strand bias", {
    filt <- adaptiveFilter("ont")
    expect_s4_class(filt, "ISNVFilter")
    expect_equal(slot(filt, "minFreq"), 0.10)
    expect_equal(slot(filt, "minDepth"), 30L)
    expect_equal(slot(filt, "maxStrandBias"), Inf)
})

test_that("adaptiveFilter('pacbio_hifi') returns intermediate thresholds", {
    filt <- adaptiveFilter("pacbio_hifi")
    expect_s4_class(filt, "ISNVFilter")
    expect_equal(slot(filt, "minFreq"), 0.05)
    expect_equal(slot(filt, "minDepth"), 50L)
    # Strand bias remains at default (HiFi preserves strand info)
    expect_equal(slot(filt, "maxStrandBias"), 10.0)
})

test_that("adaptiveFilter('generic') matches ISNVFilter() defaults", {
    filt <- adaptiveFilter("generic")
    default <- ISNVFilter()
    expect_s4_class(filt, "ISNVFilter")
    expect_equal(slot(filt, "minFreq"), slot(default, "minFreq"))
    expect_equal(slot(filt, "minDepth"), slot(default, "minDepth"))
    expect_equal(slot(filt, "maxFreq"), slot(default, "maxFreq"))
    expect_equal(slot(filt, "minAltReads"), slot(default, "minAltReads"))
    expect_equal(slot(filt, "maxStrandBias"), slot(default, "maxStrandBias"))
})
