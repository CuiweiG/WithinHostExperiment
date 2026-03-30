test_that("fitDiversityModel returns expected structure", {
    df <- data.frame(
        shannon = c(1.2, 1.5, 0.8, 1.1, 1.9, 0.7),
        timepoint = c(1, 2, 3, 1, 2, 3)
    )
    res <- fitDiversityModel(df, shannon ~ timepoint)
    expect_type(res, "list")
    expect_named(res, c("model", "summary_df", "residuals"))
    expect_s3_class(res$model, "lm")
    expect_s3_class(res$summary_df, "data.frame")
    expect_true("term" %in% colnames(res$summary_df))
    expect_true("p_value" %in% colnames(res$summary_df))
    expect_length(res$residuals, nrow(df))
})

test_that("fitDiversityModel adds group_col as fixed effect", {
    df <- data.frame(
        shannon = c(1.2, 1.5, 0.8, 1.1, 1.9, 0.7),
        timepoint = c(1, 2, 3, 1, 2, 3),
        group = c("A", "A", "A", "B", "B", "B")
    )
    res <- fitDiversityModel(df, shannon ~ timepoint, group_col = "group")
    ## group term should appear in coefficients
    expect_true(any(grepl("group", res$summary_df$term)))
})

test_that("multitestCorrection returns adjusted p-values", {
    raw_p <- c(0.001, 0.04, 0.03, 0.5)
    adj <- multitestCorrection(raw_p)
    expect_length(adj, length(raw_p))
    expect_true(all(adj >= raw_p))
    ## Bonferroni should be more conservative than BH
    adj_bon <- multitestCorrection(raw_p, method = "bonferroni")
    expect_true(all(adj_bon >= adj))
})

test_that("permutationTest detects clear group difference", {
    set.seed(123)
    vals <- c(rep(0, 20), rep(10, 20))
    grps <- rep(c("A", "B"), each = 20)
    res <- permutationTest(vals, grps, n_perm = 499)
    expect_type(res, "list")
    expect_named(res, c("observed", "p_value", "n_perm"))
    expect_equal(res$observed, -10)
    expect_true(res$p_value < 0.05)
    expect_equal(res$n_perm, 499L)
})

test_that("permutationTest rejects non-two-group input", {
    expect_error(
        permutationTest(1:9, rep(c("A", "B", "C"), each = 3)),
        "exactly 2 unique levels"
    )
    expect_error(
        permutationTest(1:3, c("A", "A", "A")),
        "exactly 2 unique levels"
    )
})
