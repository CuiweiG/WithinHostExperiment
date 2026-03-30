.make_pair_whe <- function() {
    vcf1 <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_recipient.vcf",
                        package = "WithinHostExperiment")
    readWithinHost(c(vcf1, vcf2),
        colData = S4Vectors::DataFrame(
            sample_id = c("donor", "recipient"),
            role      = c("donor", "recipient"),
            pair_id   = c("p1", "p1")),
        caller = "ivar")
}

test_that("plotFrequencySpectrum returns ggplot", {
    skip_if_not_installed("ggplot2")
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "d"),
        caller = "ivar")
    p <- plotFrequencySpectrum(whe)
    expect_s3_class(p, "ggplot")
})

test_that("plotPairScatter returns ggplot", {
    skip_if_not_installed("ggplot2")
    whe <- .make_pair_whe()
    p   <- plotPairScatter(whe, pairId = "p1")
    expect_s3_class(p, "ggplot")
})

test_that("plotQCDashboard returns patchwork", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "d"),
        caller = "ivar")
    whe <- flagISNV(whe, ISNVFilter())
    p   <- plotQCDashboard(whe)
    expect_s3_class(p, c("patchwork", "ggplot"))
})
