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

test_that("quickBottleneck returns valid result", {
    whe <- .make_pair_whe()
    nb  <- quickBottleneck(whe, pairId = "p1", maxNb = 50L)
    expect_true(is.list(nb))
    expect_true(is.finite(nb$Nb))
    expect_true(nb$Nb >= 1)
    expect_equal(length(nb$ci), 2L)
    expect_true(nb$ci[1] <= nb$Nb)
    expect_true(nb$ci[2] >= nb$Nb)
    expect_true(nb$n_variants > 0L)
})

test_that(".bbSiteLogLik returns finite value", {
    ll <- WithinHostExperiment:::.bbSiteLogLik(0.3, 0.25, 5L, 0.03)
    expect_true(is.finite(ll))
    # Log-density can be positive for beta distribution
    expect_true(is.numeric(ll))
})

test_that(".logSumExp is numerically stable", {
    x      <- c(-1000, -1001, -999)
    result <- WithinHostExperiment:::.logSumExp(x)
    expect_true(is.finite(result))
    expect_true(result > -1000)
})

test_that("quickBottleneck handles no shared variants gracefully", {
    vcf  <- system.file("extdata", "test_donor.vcf",
                        package = "WithinHostExperiment")
    vcf2 <- system.file("extdata", "test_independent.vcf",
                        package = "WithinHostExperiment")
    whe <- readWithinHost(c(vcf, vcf2),
        colData = S4Vectors::DataFrame(
            sample_id = c("d", "ind"),
            role      = c("donor", "recipient"),
            pair_id   = c("p2", "p2")),
        caller = "ivar")
    nb <- quickBottleneck(whe, pairId = "p2", maxNb = 20L)
    expect_true(is.list(nb))
})
