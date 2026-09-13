# R/bridge.R
# Bridge to downstream tools (ViralBottleneck, etc.)

#' @include methods.R
#' @importFrom SummarizedExperiment assay assayNames colData rowRanges
#' @importFrom GenomicRanges start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors DataFrame mcols
NULL

#' Export donor-recipient frequency table
#'
#' For each transmission pair, creates a table of variant frequencies
#' in the donor and recipient.
#'
#' @param whe A \code{\link{WithinHostExperiment}} with \code{role}
#'   and \code{pair_id} in \code{colData}.
#' @param pairId Character (optional). Specific pair to export.
#'   If NULL, exports all pairs.
#' @param usePassedOnly Logical. Use only QC-passed variants
#'   (default: TRUE).
#'
#' @return A \code{DataFrame} with columns: chrom, position, ref, alt,
#'   donor_freq, recipient_freq, shared (logical).
#'
#' @export
#' @examples
#' vcf1 <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' vcf2 <- system.file("extdata", "test_recipient.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(c(vcf1, vcf2),
#'     colData = S4Vectors::DataFrame(
#'         sample_id = c("donor", "recipient"),
#'         role = c("donor", "recipient"),
#'         pair_id = c("p1", "p1")),
#'     caller = "ivar")
#' exportPairFrequencies(whe)
exportPairFrequencies <- function(whe, pairId = NULL,
                                  usePassedOnly = TRUE) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    pairs <- transmissionPairs(whe)
    if (nrow(pairs) == 0L) stop("No transmission pairs defined in colData")

    if (!is.null(pairId)) {
        pairs <- pairs[pairs$pair_id == pairId, , drop = FALSE]
        if (nrow(pairs) == 0L) stop("Pair '", pairId, "' not found")
    }

    freq_mat <- assay(whe, "altFreq")
    qc_mat <- if ("qcPass" %in% assayNames(whe) && usePassedOnly) {
        assay(whe, "qcPass")
    } else {
        NULL
    }

    sample_ids <- as.character(colData(whe)$sample_id)
    rr <- rowRanges(whe)
    ref_col <- if ("ref" %in% colnames(mcols(rr))) {
        mcols(rr)$ref
    } else {
        rep(NA_character_, length(rr))
    }
    alt_col <- if ("alt" %in% colnames(mcols(rr))) {
        mcols(rr)$alt
    } else {
        rep(NA_character_, length(rr))
    }

    result_list <- lapply(seq_len(nrow(pairs)), function(p) {
        d_col <- match(as.character(pairs$donor[p]), sample_ids)
        r_col <- match(as.character(pairs$recipient[p]), sample_ids)

        d_freq <- freq_mat[, d_col]
        r_freq <- freq_mat[, r_col]

        ## Apply QC mask
        if (!is.null(qc_mat)) {
            d_freq[!qc_mat[, d_col]] <- NA
            r_freq[!qc_mat[, r_col]] <- NA
        }

        ## Keep rows where at least one has data
        has_data <- !is.na(d_freq) | !is.na(r_freq)

        DataFrame(
            chrom          = as.character(seqnames(rr))[has_data],
            position       = start(rr)[has_data],
            ref            = ref_col[has_data],
            alt            = alt_col[has_data],
            donor_freq     = d_freq[has_data],
            recipient_freq = r_freq[has_data],
            shared         = !is.na(d_freq[has_data]) &
                !is.na(r_freq[has_data]),
            pair_id        = as.character(pairs$pair_id[p])
        )
    })

    do.call(rbind, result_list)
}

#' Export data for ViralBottleneck package
#'
#' Converts a \code{\link{WithinHostExperiment}} to the input format
#' expected by the \pkg{ViralBottleneck} package.
#'
#' @param whe A \code{\link{WithinHostExperiment}}.
#' @param pairId Character (optional). Specific pair.
#' @param threshold Numeric. Variant calling threshold (default: 0.03).
#' @param usePassedOnly Logical (default: TRUE).
#'
#' @return A \code{data.frame} with columns: pos, donor_freq,
#'   recipient_freq -- a simplified format for beta-binomial MLE.
#'   For the full ViralBottleneck input format (per-base read
#'   counts), prepare data externally using the read count
#'   assays (\code{altCount}, \code{refCount}).
#'
#' @details
#' For comprehensive bottleneck analysis with multiple methods
#' (presence-absence, KL divergence, binomial, beta-binomial
#' approximate/exact, and Wright-Fisher), use the
#' \pkg{ViralBottleneck} package (Zheng et al. 2025,
#' \emph{Virus Evolution} 11:veaf071). This function prepares
#' a simplified input suitable for the beta-binomial methods.
#'
#' @export
#' @examples
#' vcf1 <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' vcf2 <- system.file("extdata", "test_recipient.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(c(vcf1, vcf2),
#'     colData = S4Vectors::DataFrame(
#'         sample_id = c("donor", "recipient"),
#'         role = c("donor", "recipient"),
#'         pair_id = c("p1", "p1")),
#'     caller = "ivar")
#' asViralBottleneckInput(whe, pairId = "p1")
asViralBottleneckInput <- function(whe, pairId = NULL, threshold = 0.03,
                                   usePassedOnly = TRUE) {
    pf <- exportPairFrequencies(whe, pairId = pairId,
                                usePassedOnly = usePassedOnly)

    ## Filter to donor variants above threshold
    donor_above <- !is.na(pf$donor_freq) & pf$donor_freq >= threshold
    pf <- pf[donor_above, , drop = FALSE]

    ## Replace NA recipient freq with 0
    recip_freq <- pf$recipient_freq
    recip_freq[is.na(recip_freq)] <- 0

    data.frame(
        pos            = pf$position,
        donor_freq     = pf$donor_freq,
        recipient_freq = recip_freq,
        stringsAsFactors = FALSE
    )
}

#' Calculate shared variant statistics
#'
#' @param whe A \code{\link{WithinHostExperiment}}.
#' @param pairId Character. Which pair to analyse.
#'
#' @return A named list: n_donor_only, n_recipient_only, n_shared,
#'   n_total, sharing_rate.
#'
#' @export
#' @examples
#' vcf1 <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' vcf2 <- system.file("extdata", "test_recipient.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(c(vcf1, vcf2),
#'     colData = S4Vectors::DataFrame(
#'         sample_id = c("donor", "recipient"),
#'         role = c("donor", "recipient"),
#'         pair_id = c("p1", "p1")),
#'     caller = "ivar")
#' calcSharedVariants(whe, pairId = "p1")
calcSharedVariants <- function(whe, pairId) {
    pf <- exportPairFrequencies(whe, pairId = pairId)

    d_only <- sum(!is.na(pf$donor_freq) & is.na(pf$recipient_freq))
    r_only <- sum(is.na(pf$donor_freq) & !is.na(pf$recipient_freq))
    shared <- sum(pf$shared)
    total  <- nrow(pf)

    list(
        n_donor_only     = d_only,
        n_recipient_only = r_only,
        n_shared         = shared,
        n_total          = total,
        sharing_rate     = if (total > 0L) shared / total else NA_real_
    )
}
