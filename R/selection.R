# R/selection.R
# Within-host dN/dS and selection coefficients

#' @include methods.R
#' @importFrom SummarizedExperiment assay assayNames colData rowRanges
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom GenomicRanges start
#' @importFrom GenomeInfoDb seqnames
NULL

# ============================================================
# dndsWithinHost -- dN/dS ratio per sample
# ============================================================

#' Within-host dN/dS ratio
#'
#' Classifies each iSNV as synonymous (REF_AA == ALT_AA) or
#' nonsynonymous and computes the dN/dS ratio per sample. When
#' a \code{GFF_FEATURE} column is present in \code{mcols(rowRanges)},
#' per-gene ratios are also returned.
#'
#' The dN/dS ratio (omega) is the ratio of nonsynonymous to synonymous
#' substitution rates. Under neutrality, dN/dS = 1. Values > 1
#' indicate positive (diversifying) selection; values < 1 indicate
#' purifying (negative) selection. This implementation follows the
#' counting method of Nei & Gojobori (1986).
#'
#' @param whe A \code{\link{WithinHostExperiment}} with amino acid
#'   annotation in \code{mcols(rowRanges)}.
#' @param gffFeatureCol Character scalar. Column name in
#'   \code{mcols(rowRanges)} containing the gene/feature annotation
#'   (default: \code{"GFF_FEATURE"}).
#' @param refAACol Character scalar. Column name in
#'   \code{mcols(rowRanges)} for the reference amino acid
#'   (default: \code{"REF_AA"}).
#' @param altAACol Character scalar. Column name in
#'   \code{mcols(rowRanges)} for the alternative amino acid
#'   (default: \code{"ALT_AA"}).
#' @param usePassedOnly Logical. If TRUE (default), only use variants
#'   where \code{qcPass == TRUE}.
#'
#' @return A \code{DataFrame} with columns:
#'   \describe{
#'     \item{sample_id}{Sample identifier.}
#'     \item{nS}{Number of synonymous iSNVs.}
#'     \item{nN}{Number of nonsynonymous iSNVs.}
#'     \item{dNdS}{dN/dS ratio (NA if nS == 0).}
#'   }
#'   If \code{gffFeatureCol} is present, additional columns
#'   \code{gene}, \code{gene_nS}, \code{gene_nN}, and
#'   \code{gene_dNdS} are included (one row per sample-gene
#'   combination).
#'
#' @references
#' Nei M, Gojobori T (1986). Simple methods for estimating the
#' numbers of synonymous and nonsynonymous nucleotide substitutions.
#' \emph{Mol Biol Evol} 3:418-426.
#' \doi{10.1093/oxfordjournals.molbev.a040410}
#'
#' @export
#' @examples
#' gr <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(100, 200, 300), width = 1))
#' S4Vectors::mcols(gr)$ref <- c("A", "C", "G")
#' S4Vectors::mcols(gr)$alt <- c("T", "G", "A")
#' S4Vectors::mcols(gr)$REF_AA <- c("M", "L", "A")
#' S4Vectors::mcols(gr)$ALT_AA <- c("I", "L", "V")
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(sample_id = "S1"))
#' dndsWithinHost(whe)
dndsWithinHost <- function(whe,
                           gffFeatureCol = "GFF_FEATURE",
                           refAACol = "REF_AA",
                           altAACol = "ALT_AA",
                           usePassedOnly = TRUE) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")

    rr <- rowRanges(whe)
    mc <- mcols(rr)

    if (!refAACol %in% colnames(mc))
        stop("Column '", refAACol, "' not found in mcols(rowRanges). ",
             "Available columns: ",
             paste(colnames(mc), collapse = ", "))
    if (!altAACol %in% colnames(mc))
        stop("Column '", altAACol, "' not found in mcols(rowRanges). ",
             "Available columns: ",
             paste(colnames(mc), collapse = ", "))

    ref_aa <- as.character(mc[[refAACol]])
    alt_aa <- as.character(mc[[altAACol]])
    is_syn <- !is.na(ref_aa) & !is.na(alt_aa) & ref_aa == alt_aa

    freq_mat <- assay(whe, "altFreq")
    qc_mat <- if (usePassedOnly && "qcPass" %in% assayNames(whe)) {
        assay(whe, "qcPass")
    } else {
        NULL
    }

    sample_ids <- colnames(freq_mat)
    if (is.null(sample_ids))
        sample_ids <- paste0("S", seq_len(ncol(freq_mat)))

    has_gene <- gffFeatureCol %in% colnames(mc)
    gene_vals <- if (has_gene) as.character(mc[[gffFeatureCol]]) else NULL

    results <- lapply(seq_len(ncol(freq_mat)), function(j) {
        freq_j <- freq_mat[, j]
        if (!is.null(qc_mat)) freq_j[!qc_mat[, j]] <- NA
        present <- !is.na(freq_j) & freq_j > 0

        nS <- sum(present & is_syn)
        nN <- sum(present & !is_syn)
        dNdS <- if (nS > 0L) nN / nS else NA_real_

        if (!has_gene) {
            return(data.frame(
                sample_id = sample_ids[j],
                nS = nS, nN = nN, dNdS = dNdS,
                stringsAsFactors = FALSE))
        }

        ## Per-gene breakdown
        genes <- unique(gene_vals[present & !is.na(gene_vals)])
        if (length(genes) == 0L) {
            return(data.frame(
                sample_id = sample_ids[j],
                nS = nS, nN = nN, dNdS = dNdS,
                gene = NA_character_,
                gene_nS = NA_integer_, gene_nN = NA_integer_,
                gene_dNdS = NA_real_,
                stringsAsFactors = FALSE))
        }

        gene_rows <- vapply(genes, function(g) {
            in_gene <- present & !is.na(gene_vals) & gene_vals == g
            gnS <- sum(in_gene & is_syn)
            gnN <- sum(in_gene & !is_syn)
            gdNdS <- if (gnS > 0L) gnN / gnS else NA_real_
            c(gnS, gnN, gdNdS)
        }, numeric(3))

        data.frame(
            sample_id = rep(sample_ids[j], length(genes)),
            nS = rep(nS, length(genes)),
            nN = rep(nN, length(genes)),
            dNdS = rep(dNdS, length(genes)),
            gene = genes,
            gene_nS = as.integer(gene_rows[1L, ]),
            gene_nN = as.integer(gene_rows[2L, ]),
            gene_dNdS = gene_rows[3L, ],
            stringsAsFactors = FALSE)
    })

    DataFrame(do.call(rbind, results))
}


# ============================================================
# estimateSelectionCoefficient
# ============================================================

#' Estimate selection coefficient from frequency trajectory
#'
#' Given a vector of allele frequencies observed over equally-spaced
#' time points, estimates the selection coefficient \eqn{s} using
#' the change in logit-transformed frequency per unit time:
#' \deqn{s = \frac{\Delta \mathrm{logit}(p)}{\Delta t}}
#' where \eqn{\mathrm{logit}(p) = \log(p / (1 - p))}.
#'
#' This assumes a simple deterministic selection model. Each
#' consecutive pair of time points yields one estimate of \eqn{s}.
#'
#' @param freqs Numeric vector. Allele frequencies at successive
#'   time points (values must be in (0, 1)).
#' @param times Numeric vector (optional). Time values for each
#'   observation. If \code{NULL}, integer time steps 1, 2, ... are
#'   used.
#' @param generation_time Numeric scalar. Time per generation
#'   (default: 1). Used to scale the selection coefficient.
#'
#' @return Numeric vector of selection coefficient estimates, one
#'   per consecutive time interval. Length is
#'   \code{length(freqs) - 1}.
#'
#' @references
#' Feder AF et al. (2014). More effective drugs lead to harder
#' selective sweeps in the evolution of drug resistance in HIV-1.
#' \emph{eLife} 3:e04492. \doi{10.7554/eLife.04492}
#'
#' @export
#' @examples
#' freqs <- c(0.05, 0.10, 0.25, 0.50)
#' estimateSelectionCoefficient(freqs)
#' estimateSelectionCoefficient(freqs, times = c(0, 3, 7, 14))
estimateSelectionCoefficient <- function(freqs, times = NULL,
                                         generation_time = 1) {
    if (!is.numeric(freqs) || length(freqs) < 2L)
        stop("'freqs' must be a numeric vector with at least 2 elements.")
    if (any(is.na(freqs)))
        stop("'freqs' must not contain NA values.")
    if (any(freqs <= 0 | freqs >= 1))
        stop("'freqs' must contain values strictly between 0 and 1.")
    if (!is.numeric(generation_time) || length(generation_time) != 1L ||
        generation_time <= 0)
        stop("'generation_time' must be a single positive number.")

    if (is.null(times)) {
        times <- seq_along(freqs)
    } else {
        if (!is.numeric(times) || length(times) != length(freqs))
            stop("'times' must be a numeric vector with the same ",
                 "length as 'freqs'.")
    }

    logit_p <- log(freqs / (1 - freqs))
    n <- length(freqs)
    delta_logit <- logit_p[2:n] - logit_p[1:(n - 1L)]
    delta_t <- (times[2:n] - times[1:(n - 1L)]) / generation_time

    if (any(delta_t == 0))
        stop("'times' must contain strictly increasing values.")

    delta_logit / delta_t
}
