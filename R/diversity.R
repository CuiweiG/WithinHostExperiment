# R/diversity.R
# Within-host diversity indices

#' @include methods.R
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom S4Vectors DataFrame
#' @importFrom stats median wilcox.test kruskal.test
NULL

# ============================================================
# Standalone math functions (exported, work on plain vectors)
# ============================================================

#' Shannon entropy from iSNV frequencies
#'
#' Calculates per-site biallelic Shannon entropy, summed across
#' all iSNV sites:
#' \deqn{H = -\sum_i [p_i \log(p_i) + (1-p_i) \log(1-p_i)]}
#' where \eqn{p_i} is the alternative allele frequency at site
#' \eqn{i}. This formulation treats each iSNV as an independent
#' biallelic locus, following McCrone & Lauring (2018) and Popa
#' et al. (2020). Higher values indicate greater within-host
#' diversity.
#'
#' @param freq Numeric vector. Alternative allele frequencies
#'   (values in (0, 1); 0, 1, and NA are excluded).
#'
#' @return Numeric scalar. Shannon entropy.
#'
#' @references
#' Lauring AS (2020). Within-Host Viral Diversity: A Window into
#' Viral Evolution. \emph{Annu Rev Virol} 7:63-81.
#' \doi{10.1146/annurev-virology-010320-061642}
#'
#' @export
#' @examples
#' shannonISNV(c(0.05, 0.15, 0.30, 0.45))
#' shannonISNV(0.5) # log(2)
shannonISNV <- function(freq) {
    freq <- as.numeric(freq)
    freq <- freq[!is.na(freq) & freq > 0 & freq < 1]
    if (length(freq) == 0L) return(0)
    -sum(freq * log(freq) + (1 - freq) * log(1 - freq))
}

#' Nucleotide diversity (pi) from iSNV frequencies
#'
#' Estimates nucleotide diversity following Nei & Li (1979),
#' adapted for deep sequencing allele frequencies:
#' \deqn{\pi = \frac{n}{n-1} \cdot \frac{2}{L} \sum_i p_i(1 - p_i)}
#' where \eqn{p_i} is the alternative allele frequency at site
#' \eqn{i}, \eqn{L} is the genome length, and \eqn{n} is the
#' sample size (approximated by mean read depth for deep
#' sequencing data). The \eqn{n/(n-1)} correction removes
#' finite-sample bias (Nei 1987, Eq. 10.5). When
#' \code{meanDepth} is not provided, the uncorrected estimator
#' is returned (equivalent to infinite sample size).
#'
#' @param freq Numeric vector. Alternative allele frequencies.
#' @param genomeLength Integer scalar. Reference genome length in bp.
#' @param meanDepth Numeric scalar (optional). Mean read depth for
#'   finite-sample correction. If \code{NULL}, no correction is
#'   applied.
#'
#' @return Numeric scalar. Nucleotide diversity.
#'
#' @references
#' Nei M, Li WH (1979). Mathematical model for studying genetic
#' variation in terms of restriction endonucleases. \emph{Proc
#' Natl Acad Sci USA} 76:5269-5273.
#' \doi{10.1073/pnas.76.10.5269}
#'
#' Nei M (1987). \emph{Molecular Evolutionary Genetics.}
#' Columbia University Press, New York. Chapter 10.
#'
#' @export
#' @examples
#' piISNV(c(0.1, 0.3, 0.5), genomeLength = 13588)
#' piISNV(c(0.1, 0.3, 0.5), genomeLength = 13588, meanDepth = 1000)
piISNV <- function(freq, genomeLength, meanDepth = NULL) {
    if (!is.numeric(freq))
        stop("'freq' must be a numeric vector.")
    if (!is.numeric(genomeLength) || length(genomeLength) != 1L ||
        genomeLength <= 0)
        stop("'genomeLength' must be a single positive number.")
    freq <- freq[!is.na(freq)]
    if (length(freq) == 0L) return(0)
    raw_pi <- (2 / genomeLength) * sum(freq * (1 - freq))
    ## Apply finite-sample correction n/(n-1) (Nei 1987)
    if (!is.null(meanDepth)) {
        if (!is.numeric(meanDepth) || length(meanDepth) != 1L ||
            meanDepth <= 1)
            stop("'meanDepth' must be a single number > 1.")
        n <- as.integer(round(meanDepth))
        raw_pi <- raw_pi * n / (n - 1L)
    }
    raw_pi
}

#' Watterson's theta from iSNV count
#'
#' Estimates Watterson's theta (Watterson 1975):
#' \deqn{\hat{\theta}_W = S / (a_n \cdot L)}
#' where \eqn{S} is the number of segregating sites,
#' \eqn{a_n = \sum_{i=1}^{n-1} 1/i} is the harmonic number,
#' and \eqn{L} is the genome length. For deep sequencing data,
#' the mean read depth is used as an approximation for the
#' sample size \eqn{n} (the number of sampled viral genomes).
#'
#' @param nSites Integer scalar. Number of segregating (iSNV) sites.
#' @param genomeLength Integer scalar. Genome length in bp.
#' @param meanDepth Numeric scalar. Mean read depth (approximates
#'   sample size n).
#'
#' @return Numeric scalar. Watterson's theta estimate.
#'
#' @references
#' Watterson GA (1975). On the number of segregating sites in
#' genetical models without recombination. \emph{Theor Popul
#' Biol} 7:256-276.
#' \doi{10.1016/0040-5809(75)90020-9}
#'
#' @export
#' @examples
#' wattersonISNV(15, genomeLength = 13588, meanDepth = 2000)
wattersonISNV <- function(nSites, genomeLength, meanDepth) {
    if (!is.numeric(nSites) || length(nSites) != 1L || nSites < 0)
        stop("'nSites' must be a single non-negative number.")
    if (!is.numeric(genomeLength) || length(genomeLength) != 1L ||
        genomeLength <= 0)
        stop("'genomeLength' must be a single positive number.")
    if (!is.numeric(meanDepth) || length(meanDepth) != 1L ||
        meanDepth <= 1)
        stop("'meanDepth' must be a single number > 1.")
    if (nSites == 0) return(0)
    n <- max(as.integer(round(meanDepth)), 2L)
    a_n <- sum(1 / seq_len(n - 1L))
    nSites / (a_n * genomeLength)
}

#' Simpson diversity index for iSNV frequencies
#'
#' Calculates Simpson diversity for biallelic iSNV sites:
#' \deqn{D = 1 - \sum_i [p_i^2 + (1 - p_i)^2]}
#' where \eqn{p_i} is the alternative allele frequency at site
#' \eqn{i}. Each site contributes two allele probabilities
#' (\eqn{p} and \eqn{1-p}), and the index measures the
#' probability that two randomly drawn alleles differ. Higher
#' values indicate greater diversity.
#'
#' @param freq Numeric vector. Alternative allele frequencies
#'   (values in (0, 1); 0, 1, and NA are excluded).
#'
#' @return Numeric scalar. Simpson diversity index.
#'
#' @export
#' @examples
#' simpsonISNV(c(0.05, 0.15, 0.30, 0.45))
#' simpsonISNV(0.5) # maximum per-site diversity = 0.5
simpsonISNV <- function(freq) {
    freq <- as.numeric(freq)
    freq <- freq[!is.na(freq) & freq > 0 & freq < 1]
    if (length(freq) == 0L) return(0)
    1 - sum(freq^2 + (1 - freq)^2)
}

#' Chao1 richness estimator for iSNV data
#'
#' Estimates total species (variant) richness using the Chao1
#' estimator (Chao 1984):
#' \deqn{\hat{S}_{Chao1} = S_{obs} + \frac{f_1^2}{2 f_2}}
#' where \eqn{S_{obs}} is the observed number of variants,
#' \eqn{f_1} is the number of singletons (variants with
#' \code{altCount == 1}), and \eqn{f_2} is the number of
#' doubletons (variants with \code{altCount == 2}). When
#' \eqn{f_2 = 0}, the bias-corrected form is used:
#' \eqn{S_{obs} + f_1 (f_1 - 1) / 2}.
#'
#' @param counts Integer vector. Alternative allele read counts
#'   (\code{altCount}) for each variant.
#' @param detected Integer scalar. Number of observed variants
#'   (\eqn{S_{obs}}).
#'
#' @return Numeric scalar. Chao1 richness estimate.
#'
#' @references
#' Chao A (1984). Nonparametric estimation of the number of
#' classes in a population. \emph{Scand J Statist} 11:265-270.
#'
#' @export
#' @examples
#' chao1ISNV(c(1, 1, 2, 5, 10), detected = 5)
chao1ISNV <- function(counts, detected) {
    counts <- as.integer(counts)
    counts <- counts[!is.na(counts)]
    if (length(counts) == 0L) return(as.numeric(detected))
    f1 <- sum(counts == 1L)
    f2 <- sum(counts == 2L)
    if (f2 == 0L) {
        detected + f1 * (f1 - 1) / 2
    } else {
        detected + f1^2 / (2 * f2)
    }
}

# ============================================================
# High-level interface (operates on WHE)
# ============================================================

#' Compute within-host diversity indices
#'
#' Calculates diversity metrics per sample from a
#' \code{\link{WithinHostExperiment}}. Implements standard
#' population genetics indices adapted for deep sequencing
#' iSNV data, following Lauring (2020).
#'
#' @references
#' Lauring AS (2020). Within-Host Viral Diversity: A Window into
#' Viral Evolution. \emph{Annu Rev Virol} 7:63-81.
#' \doi{10.1146/annurev-virology-010320-061642}
#'
#' @param x A \code{\link{WithinHostExperiment}}.
#' @param indices Character vector. Which indices to compute:
#'   \code{"shannon"}, \code{"pi"}, \code{"watterson"},
#'   \code{"richness"}, \code{"simpson"}, \code{"chao1"}
#'   (default: all six).
#' @param genomeLength Integer (optional). Required for \code{"pi"}
#'   and \code{"watterson"}.
#' @param usePassedOnly Logical. If TRUE (default), only use variants
#'   where \code{qcPass == TRUE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{DataFrame} with one row per sample and columns
#'   for each requested index.
#'
#' @export
#' @rdname calcDiversity
#' @aliases calcDiversity,WithinHostExperiment-method
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' calcDiversity(whe, genomeLength = 1000L)
setMethod("calcDiversity", "WithinHostExperiment",
    function(x,
             indices = c("shannon", "pi",
                         "watterson", "richness",
                         "simpson", "chao1"),
             genomeLength = NULL,
             usePassedOnly = TRUE, ...) {
    whe <- x
    indices <- match.arg(indices, several.ok = TRUE)

    needs_L <- any(c("pi", "watterson") %in% indices)
    if (needs_L && is.null(genomeLength)) {
        stop("'genomeLength' is required when computing '",
             paste(intersect(indices, c("pi", "watterson")),
                   collapse = "', '"),
             "' indices. Provide it as: ",
             "calcDiversity(whe, genomeLength = <genome_length>)")
    }

    freq_mat  <- assay(whe, "altFreq")
    qc_mat    <- if ("qcPass" %in% assayNames(whe)) {
        assay(whe, "qcPass")
    } else {
        NULL
    }
    if ("watterson" %in% indices &&
        !"totalDepth" %in% assayNames(whe)) {
        warning("No 'totalDepth' assay available; using default ",
                "meanDepth = 1000 for Watterson's theta. ",
                "This may produce inaccurate estimates.",
                call. = FALSE)
    }
    depth_mat <- if ("totalDepth" %in% assayNames(whe)) {
        assay(whe, "totalDepth")
    } else {
        NULL
    }
    if ("chao1" %in% indices &&
        !"altCount" %in% assayNames(whe)) {
        warning("No 'altCount' assay available; 'chao1' will ",
                "equal observed richness.",
                call. = FALSE)
    }
    alt_count_mat <- if ("altCount" %in% assayNames(whe)) {
        assay(whe, "altCount")
    } else {
        NULL
    }

    sample_ids <- colnames(freq_mat)
    if (is.null(sample_ids)) {
        sample_ids <- paste0("S", seq_len(ncol(freq_mat)))
    }
    n_samples <- ncol(freq_mat)

    results <- lapply(seq_len(n_samples), function(j) {
        freq_j <- freq_mat[, j]
        ## Apply QC mask
        if (usePassedOnly && !is.null(qc_mat)) {
            freq_j[!qc_mat[, j]] <- NA
        }
        valid <- freq_j[!is.na(freq_j)]

        vals <- list(sample_id = sample_ids[j])
        ## Median depth for this sample (used by pi and watterson)
        md <- if (!is.null(depth_mat)) {
            median(depth_mat[, j], na.rm = TRUE)
        } else {
            NULL
        }

        if ("richness" %in% indices) vals$richness <- length(valid)
        if ("shannon" %in% indices) vals$shannon <- shannonISNV(valid)
        if ("pi" %in% indices) {
            vals$pi <- piISNV(valid, genomeLength,
                              meanDepth = md)
        }
        if ("watterson" %in% indices) {
            wd <- if (!is.null(md)) md else 1000
            vals$watterson <- wattersonISNV(length(valid),
                                            genomeLength, wd)
        }
        if ("simpson" %in% indices) vals$simpson <- simpsonISNV(valid)
        if ("chao1" %in% indices) {
            n_obs <- length(valid)
            if (!is.null(alt_count_mat)) {
                ac_j <- alt_count_mat[, j]
                if (usePassedOnly && !is.null(qc_mat)) {
                    ac_j[!qc_mat[, j]] <- NA
                }
                ac_valid <- ac_j[!is.na(ac_j) & !is.na(freq_j)]
                vals$chao1 <- chao1ISNV(ac_valid, n_obs)
            } else {
                vals$chao1 <- as.numeric(n_obs)
            }
        }
        as.data.frame(vals, stringsAsFactors = FALSE)
    })

    DataFrame(do.call(rbind, results))
})

#' Compare diversity between sample groups
#'
#' @param whe A \code{\link{WithinHostExperiment}}.
#' @param groupBy Character. Column name in \code{colData} to group by.
#' @param index Character. Which index to compare (default: "shannon").
#' @param genomeLength Integer (optional). Required if index is "pi"
#'   or "watterson".
#'
#' @return A \code{DataFrame} with test results.
#'
#' @export
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' ## Single-group, returns diversity table:
#' compareDiversity(whe, groupBy = "sample_id", index = "shannon")
compareDiversity <- function(whe, groupBy, index = "shannon",
                             genomeLength = NULL) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    cd <- colData(whe)
    if (!groupBy %in% colnames(cd))
        stop("'", groupBy, "' not found in colData. ",
             "Available columns: ",
             paste(colnames(cd), collapse = ", "))

    div <- calcDiversity(whe, indices = index, genomeLength = genomeLength)
    groups <- as.character(cd[[groupBy]])
    values <- div[[index]]
    unique_groups <- unique(groups[!is.na(groups)])

    if (length(unique_groups) < 2L) {
        message("Need at least 2 groups for comparison; ",
                "returning diversity table.")
        return(div)
    }

    if (length(unique_groups) == 2L) {
        g1 <- values[groups == unique_groups[1L]]
        g2 <- values[groups == unique_groups[2L]]
        if (length(g1) < 2L || length(g2) < 2L) {
            return(DataFrame(test = "wilcoxon", statistic = NA_real_,
                             p_value = NA_real_,
                             note = "Insufficient samples per group"))
        }
        tt <- wilcox.test(g1, g2)
        DataFrame(test = "wilcoxon", statistic = tt$statistic,
                  p_value = tt$p.value)
    } else {
        tt <- kruskal.test(values ~ factor(groups))
        DataFrame(test = "kruskal", statistic = tt$statistic,
                  p_value = tt$p.value)
    }
}
