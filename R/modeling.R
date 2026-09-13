# R/modeling.R
# Statistical modelling helpers for diversity data

#' @importFrom stats lm p.adjust coef residuals as.formula terms update
NULL

#' Fit a linear model to diversity data
#'
#' Fits a linear model (\code{\link[stats]{lm}}) relating a diversity
#' metric to one or more predictors. When \code{group_col} is supplied
#' it is included as a fixed effect (a full mixed-effects model would
#' require the \pkg{nlme} dependency).
#'
#' @param div_df A \code{data.frame} containing the diversity metrics
#'   and predictor columns (e.g.\ output of
#'   \code{as.data.frame(calcDiversity(...))}).
#' @param formula A formula or character string describing the model
#'   (e.g.\ \code{shannon ~ timepoint}).
#' @param group_col Character scalar (optional). Column name in
#'   \code{div_df} representing a grouping variable.
#'   Added as a fixed effect if not already present in \code{formula}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{model}{The fitted \code{lm} object.}
#'     \item{summary_df}{A \code{data.frame} of coefficient estimates,
#'       standard errors, t-values, and p-values.}
#'     \item{residuals}{Numeric vector of model residuals.}
#'   }
#'
#' @export
#' @examples
#' df <- data.frame(
#'     shannon = c(1.2, 1.5, 0.8, 1.1, 1.9, 0.7),
#'     timepoint = c(1, 2, 3, 1, 2, 3),
#'     group = c("A", "A", "A", "B", "B", "B")
#' )
#' res <- fitDiversityModel(df, shannon ~ timepoint)
#' res_grouped <- fitDiversityModel(df, shannon ~ timepoint,
#'     group_col = "group")
fitDiversityModel <- function(div_df, formula, group_col = NULL) {
    if (!is.data.frame(div_df))
        stop("'div_df' must be a data.frame.")
    formula <- as.formula(formula)

    if (!is.null(group_col)) {
        if (!group_col %in% colnames(div_df))
            stop("'", group_col, "' not found in div_df.")
        div_df[[group_col]] <- as.factor(div_df[[group_col]])
        ## Add group as fixed effect if not already in the formula
        fterms <- attr(stats::terms(formula), "term.labels")
        if (!group_col %in% fterms) {
            formula <- stats::update(formula,
                paste(". ~ . +", group_col))
        }
    }

    fit <- stats::lm(formula, data = div_df)
    s <- summary(fit)
    coef_df <- as.data.frame(s$coefficients)
    colnames(coef_df) <- c("estimate", "std_error", "t_value", "p_value")
    coef_df$term <- rownames(coef_df)
    rownames(coef_df) <- NULL
    coef_df <- coef_df[, c("term", "estimate", "std_error",
                            "t_value", "p_value")]

    list(model = fit,
         summary_df = coef_df,
         residuals = stats::residuals(fit))
}

#' Multiple testing correction for p-values
#'
#' Thin wrapper around \code{\link[stats]{p.adjust}} for convenient
#' multiple testing correction of p-value vectors.
#'
#' @param pvalues Numeric vector of raw p-values.
#' @param method Character scalar. Correction method passed to
#'   \code{\link[stats]{p.adjust}} (default: \code{"BH"} for
#'   Benjamini-Hochberg).
#'
#' @return Numeric vector of adjusted p-values (same length as
#'   \code{pvalues}).
#'
#' @export
#' @examples
#' raw_p <- c(0.001, 0.04, 0.03, 0.5)
#' multitestCorrection(raw_p)
#' multitestCorrection(raw_p, method = "bonferroni")
multitestCorrection <- function(pvalues, method = "BH") {
    if (!is.numeric(pvalues))
        stop("'pvalues' must be a numeric vector.")
    stats::p.adjust(pvalues, method = method)
}

#' Permutation test for two-group comparison
#'
#' Non-parametric permutation test comparing a summary statistic
#' between two groups. The observed test statistic is compared to a
#' null distribution obtained by randomly permuting group labels.
#'
#' @param values Numeric vector. The measured values (e.g.\ a diversity
#'   index).
#' @param groups Factor or character vector (same length as
#'   \code{values}). Must contain exactly two unique levels.
#' @param n_perm Integer scalar. Number of permutations (default: 999).
#'   Permutation uses the random number generator, so set a seed for
#'   a reproducible p-value.
#' @param statistic Character scalar. The test statistic to compute:
#'   \code{"mean_diff"} (default) for difference in group means.
#'
#' @return A list with components:
#'   \describe{
#'     \item{observed}{Numeric scalar. The observed test statistic.}
#'     \item{p_value}{Numeric scalar. Two-sided permutation p-value.}
#'     \item{n_perm}{Integer. Number of permutations performed.}
#'   }
#'
#' @export
#' @examples
#' set.seed(42)
#' vals <- c(rnorm(10, 2), rnorm(10, 3))
#' grps <- rep(c("A", "B"), each = 10)
#' permutationTest(vals, grps, n_perm = 499)
permutationTest <- function(values, groups, n_perm = 999L,
                            statistic = "mean_diff") {
    if (!is.numeric(values))
        stop("'values' must be a numeric vector.")
    groups <- as.character(groups)
    if (length(values) != length(groups))
        stop("'values' and 'groups' must have the same length.")
    ugroups <- unique(groups)
    if (length(ugroups) != 2L)
        stop("'groups' must contain exactly 2 unique levels, got ",
             length(ugroups), ".")
    statistic <- match.arg(statistic, choices = "mean_diff")
    n_perm <- as.integer(n_perm)

    .calc_stat <- function(v, g) {
        mean(v[g == ugroups[1L]]) - mean(v[g == ugroups[2L]])
    }

    observed <- .calc_stat(values, groups)

    n <- length(values)
    perm_stats <- vapply(seq_len(n_perm), function(i) {
        .calc_stat(values, sample(groups))
    }, numeric(1L))

    ## Two-sided p-value: proportion of permutations at least as extreme
    p_value <- (sum(abs(perm_stats) >= abs(observed)) + 1L) / (n_perm + 1L)

    list(observed = observed,
         p_value = p_value,
         n_perm = n_perm)
}
