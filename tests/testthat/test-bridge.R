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

test_that("exportPairFrequencies returns correct structure", {
    whe <- .make_pair_whe()
    pf  <- exportPairFrequencies(whe, pairId = "p1")
    expect_s4_class(pf, "DataFrame")
    expect_true(all(c("donor_freq", "recipient_freq", "shared") %in%
                        colnames(pf)))
    expect_true(nrow(pf) > 0L)
})

test_that("exportPairFrequencies identifies shared variants", {
    whe <- .make_pair_whe()
    pf  <- exportPairFrequencies(whe)
    expect_equal(sum(pf$shared), 5L)
})

test_that("asViralBottleneckInput returns data.frame", {
    whe <- .make_pair_whe()
    vb  <- asViralBottleneckInput(whe, pairId = "p1")
    expect_true(is.data.frame(vb))
    expect_true(all(c("pos", "donor_freq", "recipient_freq") %in%
                        colnames(vb)))
    expect_true(all(vb$donor_freq >= 0.03))
    expect_true(!any(is.na(vb$recipient_freq)))
})

test_that("calcSharedVariants returns correct counts", {
    whe <- .make_pair_whe()
    sv  <- calcSharedVariants(whe, pairId = "p1")
    expect_equal(sv$n_shared, 5L)
    expect_equal(sv$n_donor_only, 10L)
    expect_equal(sv$n_recipient_only, 5L)
    expect_equal(sv$n_total, 20L)
})

test_that("exportPairFrequencies errors without pairs", {
    vcf <- system.file("extdata", "test_donor.vcf",
                       package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "d"),
        caller = "ivar")
    expect_error(exportPairFrequencies(whe), "No transmission pairs")
})
