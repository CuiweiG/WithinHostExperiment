## Helper: load two caller WHEs from test VCFs
.load_caller_whes <- function() {
    vcf_ivar <- system.file("extdata", "test_donor.vcf",
        package = "WithinHostExperiment")
    vcf_lofreq <- system.file("extdata", "test_donor_lofreq.vcf",
        package = "WithinHostExperiment")
    vcf_fb <- system.file("extdata", "test_donor_freebayes.vcf",
        package = "WithinHostExperiment")

    whe_ivar <- readWithinHost(vcf_ivar,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "ivar")
    whe_lofreq <- readWithinHost(vcf_lofreq,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "lofreq")
    whe_fb <- readWithinHost(vcf_fb,
        colData = S4Vectors::DataFrame(sample_id = "donor"),
        caller = "freebayes")

    list(ivar = whe_ivar, lofreq = whe_lofreq, freebayes = whe_fb)
}

test_that("consensusISNV returns WHE with all shared sites from test VCFs", {
    whes <- .load_caller_whes()
    result <- consensusISNV(whes[c("ivar", "lofreq")], min_callers = 2L)
    expect_s4_class(result, "WithinHostExperiment")
    ## Both VCFs have the same 15 sites with identical alleles
    expect_equal(nrow(result), 15L)
    mc <- S4Vectors::mcols(SummarizedExperiment::rowRanges(result))
    expect_true(all(mc$n_callers == 2L))
    expect_true(all(grepl("ivar", mc$caller_names)))
    expect_true(all(grepl("lofreq", mc$caller_names)))
})

test_that("consensusISNV freq_spread is zero for identical-frequency callers", {
    whes <- .load_caller_whes()
    result <- consensusISNV(whes[c("ivar", "lofreq")], min_callers = 2L)
    mc <- S4Vectors::mcols(SummarizedExperiment::rowRanges(result))
    ## ivar and lofreq report the same AF values, so spread should be ~0
    expect_true(all(mc$freq_spread < 1e-6))
})

test_that("consensusISNV works with three callers and min_callers = 3", {
    whes <- .load_caller_whes()
    result <- consensusISNV(whes, min_callers = 3L)
    expect_s4_class(result, "WithinHostExperiment")
    ## All three callers share the same 15 sites
    expect_equal(nrow(result), 15L)
    mc <- S4Vectors::mcols(SummarizedExperiment::rowRanges(result))
    expect_true(all(mc$n_callers == 3L))
})

test_that("calcCallerConcordance returns correct pairwise stats from VCFs", {
    whes <- .load_caller_whes()
    result <- calcCallerConcordance(whes)
    expect_s4_class(result, "DataFrame")
    ## 3 callers -> 3 pairs
    expect_equal(nrow(result), 3L)
    expect_true(all(c("caller1", "caller2", "n_shared",
                      "n_only1", "n_only2", "jaccard") %in%
                    colnames(result)))
    ## All callers share all 15 sites -> jaccard = 1
    expect_true(all(result$jaccard == 1.0))
    expect_true(all(result$n_shared == 15L))
    expect_true(all(result$n_only1 == 0L))
    expect_true(all(result$n_only2 == 0L))
})

test_that("consensusISNV rejects unnamed or single-element lists", {
    whes <- .load_caller_whes()
    unnamed <- whes[c("ivar", "lofreq")]
    names(unnamed) <- NULL
    expect_error(consensusISNV(unnamed), "named list")
    expect_error(consensusISNV(whes[1]), "at least 2")
})
