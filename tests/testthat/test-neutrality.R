# tests/testthat/test-neutrality.R
# Neutrality test suite

test_that("tajimaD returns 0 when pi equals theta_W", {
    ## When pi_hat * L == S / a1, the numerator is zero => D = 0
    n <- 100L
    S <- 10L
    L <- 10000L
    a1 <- sum(1 / seq_len(n - 1L))
    theta_W <- S / a1
    ## Set pi_hat so that pi_hat * L == theta_W
    pi_hat <- theta_W / L
    D <- tajimaD(S = S, n = n, pi_hat = pi_hat, genomeLength = L)
    expect_equal(D, 0, tolerance = 1e-10)
})

test_that("tajimaD returns negative for excess rare variants", {
    ## Many segregating sites but low pi => excess rare variants => D < 0
    n <- 200L
    S <- 50L
    L <- 10000L
    ## pi_hat much smaller than theta_W / L
    pi_hat <- 1e-6
    D <- tajimaD(S = S, n = n, pi_hat = pi_hat, genomeLength = L)
    expect_true(is.numeric(D))
    expect_true(D < 0)
})

test_that("fusFs returns numeric", {
    result <- fusFs(S = 10, n = 50, theta_pi = 5)
    expect_true(is.numeric(result))
    expect_length(result, 1L)
})

test_that("calcNeutralityTests works on WHE from test VCF", {
    vcf <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    res <- calcNeutralityTests(whe, genomeLength = 1000L)
    expect_s4_class(res, "DataFrame")
    expect_true("sample_id" %in% colnames(res))
    expect_true("tajimaD" %in% colnames(res))
    expect_true("fusFs" %in% colnames(res))
    expect_equal(nrow(res), 1L)
    expect_true(is.numeric(res$tajimaD))
    expect_true(is.numeric(res$fusFs))
})

test_that("calcNeutralityTests errors on missing genomeLength", {
    vcf <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    whe <- readWithinHost(vcf,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    expect_error(calcNeutralityTests(whe),
        "genomeLength.*required")
})

test_that("tajimaD returns 0 when S is 0", {
    D <- tajimaD(S = 0, n = 1000, pi_hat = 0, genomeLength = 10000)
    expect_equal(D, 0)
})

test_that("fusFs returns 0 when S is 0", {
    result <- fusFs(S = 0, n = 100, theta_pi = 5)
    expect_equal(result, 0)
})

test_that("tajimaD validates inputs", {
    expect_error(tajimaD(S = -1, n = 10, pi_hat = 0.001,
        genomeLength = 1000), "non-negative")
    expect_error(tajimaD(S = 5, n = 1, pi_hat = 0.001,
        genomeLength = 1000), "must be a single number > 1")
    expect_error(tajimaD(S = 5, n = 100, pi_hat = 0.001,
        genomeLength = 0), "positive")
})

test_that("fusFs returns NA with a warning when S + 1 exceeds n", {
    expect_warning(res <- fusFs(S = 100, n = 100, theta_pi = 5),
                   "undefined")
    expect_true(is.na(res))
    expect_true(is.finite(fusFs(S = 10, n = 100, theta_pi = 5)))
})
