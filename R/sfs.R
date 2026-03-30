# R/sfs.R
# Within-host Site Frequency Spectrum (SFS) as a first-class object
#
# Design rationale:
# The SFS is the fundamental summary statistic for population genetics.
# In within-host pathogen analysis, the SFS from deep sequencing has
# unique properties: (1) the "sample size" n is the read depth, not
# the number of individuals; (2) low-frequency variants are under-
# ascertained due to the detection threshold; (3) near-fixed variants
# (freq > 0.5) represent consensus-level changes, not iSNVs.
# This class makes all of these explicit and correctable.

#' @include AllClasses.R
#' @include methods.R
#' @importFrom methods setClass setValidity new is
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom stats median pnorm dbinom
NULL

# ============================================================
# WithinHostSFS class
# ============================================================

#' WithinHostSFS: Site Frequency Spectrum for within-host variants
#'
#' An S4 class representing the site frequency spectrum (SFS) of
#' within-host variants from deep sequencing data. The SFS bins
#' variant allele frequencies into discrete categories and stores
#' metadata needed for population-genetic inference.
#'
#' @slot counts Integer vector. Number of iSNVs in each frequency
#'   bin. Length equals \code{length(breaks) - 1}.
#' @slot breaks Numeric vector. Bin boundaries (left-closed,
#'   right-open except the last bin). Default: 20 equal bins
#'   from 0 to 1.
#' @slot folded Logical scalar. If \code{TRUE}, the SFS is folded
#'   (minor allele frequency); if \code{FALSE}, unfolded
#'   (derived allele frequency).
#' @slot sampleId Character scalar. Sample identifier.
#' @slot genomeLength Integer scalar. Reference genome length.
#' @slot nSites Integer scalar. Total number of variant sites
#'   used to build the SFS.
#' @slot meanDepth Numeric scalar. Mean sequencing depth.
#' @slot threshold Numeric scalar. Detection threshold applied
#'   (e.g., 0.03). Sites below this are unobservable.
#' @slot corrected Logical scalar. Whether ascertainment bias
#'   correction has been applied.
#'
#' @seealso \code{\link{buildSFS}} for constructing from a
#'   \code{WithinHostExperiment}, \code{\link{neutralityFromSFS}}
#'   for computing all neutrality tests from the SFS.
#'
#' @aliases WithinHostSFS
#' @exportClass WithinHostSFS
.WithinHostSFS <- setClass("WithinHostSFS",
    slots = list(
        counts       = "integer",
        breaks       = "numeric",
        folded       = "logical",
        sampleId     = "character",
        genomeLength = "integer",
        nSites       = "integer",
        meanDepth    = "numeric",
        threshold    = "numeric",
        corrected    = "logical"
    ),
    prototype = list(
        counts       = integer(0),
        breaks       = numeric(0),
        folded       = FALSE,
        sampleId     = NA_character_,
        genomeLength = NA_integer_,
        nSites       = 0L,
        meanDepth    = NA_real_,
        threshold    = 0.03,
        corrected    = FALSE
    )
)

setValidity("WithinHostSFS", function(object) {
    msg <- character()
    if (length(object@counts) != length(object@breaks) - 1L)
        msg <- c(msg, "length(counts) must equal length(breaks) - 1")
    if (any(object@counts < 0L))
        msg <- c(msg, "counts must be non-negative")
    if (is.unsorted(object@breaks))
        msg <- c(msg, "breaks must be sorted")
    if (length(msg) == 0L) TRUE else msg
})

#' @rdname WithinHostSFS-class
#' @param object A \code{WithinHostSFS} object.
#' @export
setMethod("show", "WithinHostSFS", function(object) {
    cat("WithinHostSFS with", sum(object@counts), "variants in",
        length(object@counts), "bins\n")
    cat("  sample:", object@sampleId, "\n")
    cat("  folded:", object@folded, "\n")
    cat("  threshold:", object@threshold, "\n")
    cat("  corrected:", object@corrected, "\n")
    if (!is.na(object@genomeLength))
        cat("  genome:", object@genomeLength, "bp\n")
    if (!is.na(object@meanDepth))
        cat("  mean depth:", round(object@meanDepth), "x\n")
})


# ============================================================
# buildSFS -- construct from WithinHostExperiment
# ============================================================

#' Build a site frequency spectrum from a WithinHostExperiment
#'
#' Constructs a \code{\link{WithinHostSFS}} object from the allele
#' frequencies in a \code{\link{WithinHostExperiment}}. One SFS is
#' built per sample.
#'
#' @param whe A \code{\link{WithinHostExperiment}}.
#' @param sampleIdx Integer or character. Which sample (column) to
#'   use. Default: 1.
#' @param nBins Integer. Number of frequency bins (default: 20).
#' @param fold Logical. If \code{TRUE}, fold the SFS to minor
#'   allele frequency (0 to 0.5). Default: \code{TRUE}.
#' @param usePassedOnly Logical. If \code{TRUE}, only use QC-passed
#'   variants. Default: \code{TRUE}.
#' @param threshold Numeric. Detection threshold (default: 0.03).
#'   The first bin starts at this value, not 0.
#' @param genomeLength Integer (optional). Reference genome length.
#'
#' @return A \code{\link{WithinHostSFS}} object.
#'
#' @export
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' sfs <- buildSFS(whe, genomeLength = 1000L)
#' sfs
buildSFS <- function(whe, sampleIdx = 1L, nBins = 20L,
                     fold = TRUE, usePassedOnly = TRUE,
                     threshold = 0.03, genomeLength = NULL) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")

    freq_mat <- assay(whe, "altFreq")
    if (is.character(sampleIdx))
        sampleIdx <- match(sampleIdx, colnames(freq_mat))
    freq <- freq_mat[, sampleIdx]

    ## Apply QC
    if (usePassedOnly && "qcPass" %in% assayNames(whe)) {
        qc <- assay(whe, "qcPass")[, sampleIdx]
        freq[!qc] <- NA
    }
    freq <- freq[!is.na(freq)]

    ## Get sample metadata
    cd <- colData(whe)
    sid <- as.character(cd$sample_id[sampleIdx])

    ## Mean depth
    md <- if ("totalDepth" %in% assayNames(whe)) {
        median(assay(whe, "totalDepth")[, sampleIdx], na.rm = TRUE)
    } else {
        NA_real_
    }

    ## Fold if requested
    if (fold) {
        freq <- ifelse(freq > 0.5, 1 - freq, freq)
        max_f <- 0.5
    } else {
        max_f <- 1.0
    }

    ## Build bins from threshold to max
    breaks <- seq(threshold, max_f, length.out = nBins + 1L)
    counts <- as.integer(hist(freq[freq >= threshold & freq <= max_f],
                              breaks = breaks, plot = FALSE)$counts)

    gl <- if (!is.null(genomeLength)) as.integer(genomeLength) else NA_integer_

    new("WithinHostSFS",
        counts       = counts,
        breaks       = breaks,
        folded       = fold,
        sampleId     = sid,
        genomeLength = gl,
        nSites       = as.integer(sum(counts)),
        meanDepth    = md,
        threshold    = threshold,
        corrected    = FALSE)
}


#' Build SFS for all samples in a WithinHostExperiment
#'
#' Convenience wrapper that calls \code{\link{buildSFS}} for each
#' sample and returns a named list of \code{\link{WithinHostSFS}}
#' objects.
#'
#' @inheritParams buildSFS
#'
#' @return A named list of \code{\link{WithinHostSFS}} objects.
#'
#' @export
#' @examples
#' data(example_whe)
#' sfs_list <- buildSFSList(example_whe, genomeLength = 1000L)
#' length(sfs_list)
buildSFSList <- function(whe, nBins = 20L, fold = TRUE,
                         usePassedOnly = TRUE, threshold = 0.03,
                         genomeLength = NULL) {
    n <- ncol(whe)
    sfs_list <- lapply(seq_len(n), function(j) {
        buildSFS(whe, sampleIdx = j, nBins = nBins, fold = fold,
                 usePassedOnly = usePassedOnly, threshold = threshold,
                 genomeLength = genomeLength)
    })
    names(sfs_list) <- as.character(colData(whe)$sample_id)
    sfs_list
}


# ============================================================
# Ascertainment bias correction
# ============================================================

#' Correct ascertainment bias in the within-host SFS
#'
#' Deep sequencing imposes a detection threshold (e.g., 3\%): sites
#' with true frequency below this threshold are unobservable. This
#' inflates the apparent fraction of moderate-frequency variants
#' relative to rare variants.
#'
#' This function applies a correction based on the expected number
#' of unobserved sites under a neutral beta distribution model:
#' given the observed SFS and depth-dependent detectability, it
#' estimates the missing mass in the lowest frequency bins.
#'
#' @param sfs A \code{\link{WithinHostSFS}} object.
#' @param method Character. Correction method:
#'   \code{"truncation"} (default) adjusts the monomorphic site
#'   count to account for variants below threshold;
#'   \code{"binomial"} models the probability of detecting a
#'   variant at each frequency given the read depth.
#'
#' @return A corrected \code{\link{WithinHostSFS}} with
#'   \code{corrected = TRUE}.
#'
#' @details
#' The \code{"truncation"} method assumes the SFS below the
#' threshold follows the same shape as the observable portion and
#' scales the monomorphic count accordingly. The \code{"binomial"}
#' method uses the mean depth to compute
#' \eqn{P(\mathrm{detect} | \nu, n) = P(\mathrm{Alt} \ge k_{\min} | n, \nu)}{P(detect|freq,depth)}
#' where \eqn{k_{\min} = \lceil n \times \mathrm{threshold} \rceil}{k_min = ceil(n*threshold)},
#' and reweights each bin by the inverse of its detectability.
#'
#' @export
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' sfs <- buildSFS(whe, genomeLength = 1000L)
#' sfs_corr <- correctSFSBias(sfs, method = "binomial")
#' sfs_corr
correctSFSBias <- function(sfs, method = c("truncation", "binomial")) {
    if (!is(sfs, "WithinHostSFS"))
        stop("'sfs' must be a WithinHostSFS object.")
    method <- match.arg(method)

    counts <- sfs@counts
    breaks <- sfs@breaks
    n_bins <- length(counts)
    thr    <- sfs@threshold
    depth  <- sfs@meanDepth

    if (method == "truncation") {
        ## Estimate the missing fraction: assume uniform density
        ## in the first bin extends below the threshold.
        ## fraction_observable = (max_freq - threshold) / max_freq
        max_f <- max(breaks)
        frac_obs <- (max_f - thr) / max_f
        ## Scale all counts by 1/frac_obs
        counts <- as.integer(round(counts / frac_obs))

    } else if (method == "binomial") {
        if (is.na(depth) || depth <= 0)
            stop("'binomial' correction requires mean depth ",
                 "in the SFS object.")
        n <- as.integer(round(depth))
        k_min <- ceiling(n * thr)

        ## For each bin midpoint, compute P(detect)
        mids <- (breaks[-length(breaks)] + breaks[-1L]) / 2
        p_detect <- vapply(mids, function(nu) {
            ## P(alt_reads >= k_min | n, nu)
            1 - pbinom(k_min - 1L, n, nu)
        }, numeric(1))

        ## Reweight: true_count = observed / P(detect)
        p_detect[p_detect < 0.01] <- 0.01  # floor to avoid Inf
        counts <- as.integer(round(counts / p_detect))
    }

    new("WithinHostSFS",
        counts       = counts,
        breaks       = breaks,
        folded       = sfs@folded,
        sampleId     = sfs@sampleId,
        genomeLength = sfs@genomeLength,
        nSites       = as.integer(sum(counts)),
        meanDepth    = sfs@meanDepth,
        threshold    = sfs@threshold,
        corrected    = TRUE)
}


# ============================================================
# neutralityFromSFS -- all tests from one SFS
# ============================================================

#' Compute all neutrality statistics from a WithinHostSFS
#'
#' Computes Tajima's D, Fu's Fs, and Fay and Wu's H from a single
#' \code{\link{WithinHostSFS}} object. This is the recommended
#' interface for neutrality testing: build the SFS once, then
#' derive all statistics from it.
#'
#' \strong{Fay and Wu's H} measures an excess of high-frequency
#' derived variants, which can indicate positive selection or
#' genetic hitchhiking (Fay and Wu 2000). It is defined as:
#' \deqn{H = \hat{\pi} - \hat{\theta}_H}
#' where \eqn{\hat{\theta}_H = \frac{2}{n(n-1)} \sum_i i^2 \xi_i}{theta_H = 2/(n(n-1)) * sum(i^2 * xi_i)}
#' and \eqn{\xi_i} is the count of sites with \eqn{i} derived
#' alleles. A negative H indicates an excess of high-frequency
#' variants.
#'
#' @param sfs A \code{\link{WithinHostSFS}} object.
#' @param n Integer. Effective sample size for neutrality tests.
#'   For deep sequencing data, do \strong{not} use raw depth.
#'   Default: \code{min(100, meanDepth)} as recommended in the
#'   \code{\link{tajimaD}} documentation.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{S}{Number of segregating sites.}
#'     \item{n}{Effective sample size used.}
#'     \item{pi}{Nucleotide diversity (per site).}
#'     \item{theta_W}{Watterson's theta.}
#'     \item{theta_H}{Fay and Wu's theta_H.}
#'     \item{tajimaD}{Tajima's D.}
#'     \item{fusFs}{Fu's Fs.}
#'     \item{fayWuH}{Fay and Wu's H.}
#'     \item{sfs_corrected}{Whether the input SFS was bias-corrected.}
#'   }
#'
#' @references
#' Fay JC, Wu CI (2000). Hitchhiking under positive Darwinian
#' selection. \emph{Genetics} 155:1405-1413.
#'
#' @export
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' sfs <- buildSFS(whe, genomeLength = 1000L, fold = FALSE)
#' neutralityFromSFS(sfs)
neutralityFromSFS <- function(sfs, n = NULL) {
    if (!is(sfs, "WithinHostSFS"))
        stop("'sfs' must be a WithinHostSFS object.")
    if (sfs@folded)
        message("Note: Fay & Wu's H requires an unfolded SFS. ",
                "Results may not be meaningful for folded spectra.")

    S <- sfs@nSites
    if (is.null(n)) {
        n <- if (!is.na(sfs@meanDepth)) {
            min(100L, max(as.integer(round(sfs@meanDepth)), 2L))
        } else {
            100L
        }
    }
    n <- as.integer(n)
    if (n <= 1L) stop("'n' must be > 1.")

    gl <- sfs@genomeLength
    if (is.na(gl)) gl <- 1L  # per-site

    ## Reconstruct per-site pi from SFS bin midpoints
    breaks <- sfs@breaks
    mids <- (breaks[-length(breaks)] + breaks[-1L]) / 2
    counts <- sfs@counts

    ## Nucleotide diversity from SFS
    pi_per_site <- if (S > 0 && gl > 0) {
        (2 / gl) * sum(counts * mids * (1 - mids))
    } else {
        0
    }

    ## Watterson's theta
    a1 <- sum(1 / seq_len(n - 1L))
    theta_W <- if (S > 0) (S / a1) / gl else 0

    ## Fay and Wu's theta_H
    ## For deep sequencing SFS, use bin midpoints as proxies
    ## for the frequency class i/n
    theta_H <- if (S > 0) {
        freq_classes <- round(mids * n)
        (2 / (n * (n - 1))) * sum(counts * freq_classes^2) / gl
    } else {
        0
    }

    ## Tajima's D
    D <- tajimaD(S, n, pi_per_site, gl)

    ## Fu's Fs
    theta_pi <- pi_per_site * gl
    Fs <- fusFs(S, n, theta_pi)

    ## Fay and Wu's H
    H <- pi_per_site - theta_H

    list(
        S             = S,
        n             = n,
        pi            = pi_per_site,
        theta_W       = theta_W,
        theta_H       = theta_H,
        tajimaD       = D,
        fusFs         = Fs,
        fayWuH        = H,
        sfs_corrected = sfs@corrected
    )
}


# ============================================================
# compareSFS -- compare two SFS objects
# ============================================================

#' Compare two within-host SFS objects
#'
#' Computes summary statistics and a chi-squared test for
#' homogeneity between two \code{\link{WithinHostSFS}} objects
#' (e.g., before vs. after QC, or two timepoints).
#'
#' @param sfs1 A \code{\link{WithinHostSFS}} object.
#' @param sfs2 A \code{\link{WithinHostSFS}} object (same binning).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{chisq_stat}{Chi-squared statistic.}
#'     \item{chisq_p}{Chi-squared p-value.}
#'     \item{df}{Degrees of freedom.}
#'     \item{n1}{Total variants in sfs1.}
#'     \item{n2}{Total variants in sfs2.}
#'     \item{proportions1}{Proportion in each bin, sfs1.}
#'     \item{proportions2}{Proportion in each bin, sfs2.}
#'   }
#'
#' @export
#' @examples
#' data(example_whe)
#' whe1 <- example_whe
#' sfs1 <- buildSFS(whe1, genomeLength = 1000L)
#' sfs2 <- buildSFS(whe1, sampleIdx = 2L, genomeLength = 1000L)
#' compareSFS(sfs1, sfs2)
compareSFS <- function(sfs1, sfs2) {
    if (!is(sfs1, "WithinHostSFS") || !is(sfs2, "WithinHostSFS"))
        stop("Both arguments must be WithinHostSFS objects.")
    if (length(sfs1@counts) != length(sfs2@counts))
        stop("SFS objects must have the same number of bins.")

    c1 <- sfs1@counts
    c2 <- sfs2@counts
    n1 <- sum(c1)
    n2 <- sum(c2)

    ## Remove bins where both are zero
    keep <- (c1 + c2) > 0
    if (sum(keep) < 2L) {
        return(list(chisq_stat = NA_real_, chisq_p = NA_real_,
                    df = NA_integer_, n1 = n1, n2 = n2,
                    proportions1 = c1 / max(n1, 1),
                    proportions2 = c2 / max(n2, 1)))
    }

    mat <- rbind(c1[keep], c2[keep])
    test <- suppressWarnings(chisq.test(mat))

    list(
        chisq_stat   = unname(test$statistic),
        chisq_p      = test$p.value,
        df           = unname(test$parameter),
        n1           = n1,
        n2           = n2,
        proportions1 = c1 / max(n1, 1),
        proportions2 = c2 / max(n2, 1)
    )
}


# ============================================================
# plotSFS -- publication-quality SFS visualization
# ============================================================

#' Plot a within-host site frequency spectrum
#'
#' Creates a publication-quality bar plot of a
#' \code{\link{WithinHostSFS}} object.
#'
#' @param sfs A \code{\link{WithinHostSFS}} object, or a named
#'   list of \code{WithinHostSFS} objects for overlay.
#' @param normalize Logical. If \code{TRUE}, plot proportions
#'   instead of counts. Default: \code{FALSE}.
#'
#' @return A ggplot object.
#'
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     vcf <- system.file("extdata", "test_donor.vcf",
#'         package = "WithinHostExperiment")
#'     whe <- readWithinHost(vcf,
#'         colData = S4Vectors::DataFrame(sample_id = "donor"),
#'         caller = "ivar")
#'     sfs <- buildSFS(whe, genomeLength = 1000L)
#'     plotSFS(sfs)
#' }
plotSFS <- function(sfs, normalize = FALSE) {
    .requireGgplot2()

    if (is(sfs, "WithinHostSFS")) {
        sfs <- list(sfs)
        names(sfs) <- sfs[[1]]@sampleId
    }
    if (!is.list(sfs) || !all(vapply(sfs, is, logical(1), "WithinHostSFS")))
        stop("'sfs' must be a WithinHostSFS or a named list of them.")

    df_list <- lapply(names(sfs), function(nm) {
        s <- sfs[[nm]]
        breaks <- s@breaks
        mids <- (breaks[-length(breaks)] + breaks[-1L]) / 2
        y <- if (normalize) s@counts / max(sum(s@counts), 1) else s@counts
        data.frame(
            freq = mids * 100,
            count = y,
            sample = nm,
            stringsAsFactors = FALSE)
    })
    plot_df <- do.call(rbind, df_list)
    n_samples <- length(sfs)

    ylab <- if (normalize) "Proportion" else "Number of iSNVs"

    p <- ggplot2::ggplot(plot_df,
        ggplot2::aes(x = .data$freq, y = .data$count))

    if (n_samples == 1) {
        p <- p + ggplot2::geom_col(fill = .whe_pal$blue,
            colour = "white", linewidth = 0.2, alpha = 0.90,
            width = diff(sfs[[1]]@breaks)[1] * 100 * 0.85)
    } else {
        p <- p + ggplot2::geom_col(
            ggplot2::aes(fill = .data$sample),
            colour = "white", linewidth = 0.2, alpha = 0.70,
            position = "dodge",
            width = diff(sfs[[1]]@breaks)[1] * 100 * 0.85) +
            ggplot2::scale_fill_manual(name = NULL,
                values = unlist(.whe_pal[seq_len(n_samples)]))
    }

    p <- p +
        ggplot2::scale_x_continuous(
            ifelse(sfs[[1]]@folded,
                   "Minor allele frequency (%)",
                   "Alternative allele frequency (%)")) +
        ggplot2::labs(y = ylab) +
        .theme_pub(base_size = 9)

    p
}
