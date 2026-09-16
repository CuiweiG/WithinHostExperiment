# R/qc-temporal.R
# Longitudinal QC: flag transient iSNVs across timepoints

#' @include methods.R
#' @importFrom SummarizedExperiment assay assay<- assayNames colData
#' @importFrom S4Vectors DataFrame
NULL

#' Flag temporally inconsistent iSNVs in longitudinal data
#'
#' For longitudinal samples from the same host, identifies iSNVs
#' that appear at only one timepoint (transient) versus those
#' detected at multiple consecutive timepoints (persistent).
#' Transient iSNVs are more likely to be sequencing artefacts;
#' persistent ones are more likely genuine within-host variants.
#'
#' @param whe A \code{\link{WithinHostExperiment}} with longitudinal
#'   samples from the same host.
#' @param hostCol Character scalar. Column name in \code{colData}
#'   for host identifier (default: \code{"host_id"}).
#' @param timeCol Character scalar. Column name in \code{colData}
#'   for time point (default: \code{"timepoint"}).
#' @param minTimepoints Integer scalar. Minimum number of timepoints
#'   at which an iSNV must be detected to be considered persistent
#'   (default: 2).
#' @param minConsecutive Integer scalar. Minimum number of
#'   \strong{consecutive} timepoints required. If \code{NULL}
#'   (default), any \code{minTimepoints} detections suffice
#'   regardless of consecutiveness. Set to 2 or 3 for stricter
#'   filtering.
#' @param usePassedOnly Logical. If \code{TRUE} (default), only
#'   consider QC-passed entries.
#'
#' @return A \code{\link{WithinHostExperiment}} with updated
#'   \code{qcPass} assay and a new column \code{temporal_class} in
#'   \code{mcols(rowRanges)}. Sites are classified within each host;
#'   a transient site is flagged only in the samples of the hosts
#'   where it is transient. The site-level \code{temporal_class} is
#'   \code{"persistent"} if the site is persistent in at least one
#'   host, \code{"transient"} if it is transient in every host where
#'   it was detected, and \code{NA} where it was never detected.
#'
#' @details
#' This check applies only where a host was sampled more than once.
#' The biological rationale:
#' genuine within-host variants arise from viral replication and
#' should persist (or change smoothly) across consecutive
#' timepoints. An iSNV that appears once at 5\% and is absent at
#' all other timepoints is more likely a sequencing or library
#' preparation artefact.
#'
#' The function records the filtering step in \code{qcLog} with
#' step name \code{"temporal_consistency"}.
#'
#' @references
#' Farjo M et al. (2024). Within-host evolutionary dynamics and
#' tissue compartmentalization during acute SARS-CoV-2 infection.
#' \emph{J. Virol.} 98:e01618-23.
#' \doi{10.1128/jvi.01618-23}
#'
#' @export
#' @examples
#' gr <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(100, 200, 300), width = 1))
#' S4Vectors::mcols(gr)$ref <- c("A", "C", "G")
#' S4Vectors::mcols(gr)$alt <- c("T", "G", "A")
#' freq <- matrix(c(0.1, NA,  0.2,   # var1: 2 timepoints
#'                  0.05, NA,  NA,    # var2: 1 timepoint (transient)
#'                  0.15, 0.2, 0.1),  # var3: 3 timepoints
#'                nrow = 3, ncol = 3)
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = freq),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(
#'         sample_id = c("t1", "t2", "t3"),
#'         host_id = c("H1", "H1", "H1"),
#'         timepoint = c(1, 2, 3)))
#' whe <- flagTemporalInconsistency(whe, minTimepoints = 2L)
#' S4Vectors::mcols(SummarizedExperiment::rowRanges(whe))$temporal_class
flagTemporalInconsistency <- function(whe,
                                      hostCol = "host_id",
                                      timeCol = "timepoint",
                                      minTimepoints = 2L,
                                      minConsecutive = NULL,
                                      usePassedOnly = TRUE) {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")

    cd <- colData(whe)
    if (!hostCol %in% colnames(cd))
        stop("Column '", hostCol, "' not found in colData.")
    if (!timeCol %in% colnames(cd))
        stop("Column '", timeCol, "' not found in colData.")

    freq_mat <- assay(whe, "altFreq")
    qc_mat <- if (usePassedOnly && "qcPass" %in% assayNames(whe)) {
        assay(whe, "qcPass")
    } else {
        matrix(TRUE, nrow = nrow(whe), ncol = ncol(whe))
    }

    hosts <- as.character(cd[[hostCol]])
    times <- as.numeric(cd[[timeCol]])
    n_sites <- nrow(freq_mat)
    n_samples <- ncol(freq_mat)

    ## For each site and host, count how many timepoints it is detected
    total_flagged <- 0L

    ## Group by host
    unique_hosts <- unique(hosts[!is.na(hosts)])
    class_by_host <- matrix(NA_character_, nrow = n_sites,
                            ncol = length(unique_hosts),
                            dimnames = list(NULL, unique_hosts))

    for (h in unique_hosts) {
        h_idx <- which(hosts == h)
        if (length(h_idx) < 2L) next

        ## Sort by timepoint
        h_order <- h_idx[order(times[h_idx])]

        for (i in seq_len(n_sites)) {
            ## Get detection status at each timepoint (QC-aware)
            detected <- vapply(h_order, function(j) {
                !is.na(freq_mat[i, j]) && qc_mat[i, j]
            }, logical(1))

            n_detected <- sum(detected)

            if (n_detected == 0L) {
                next  # not present in this host
            }

            ## Check consecutive requirement
            is_persistent <- n_detected >= minTimepoints

            if (is_persistent && !is.null(minConsecutive)) {
                ## Check for a run of minConsecutive TRUE values
                max_run <- 0L
                current_run <- 0L
                for (d in detected) {
                    if (d) {
                        current_run <- current_run + 1L
                        if (current_run > max_run)
                            max_run <- current_run
                    } else {
                        current_run <- 0L
                    }
                }
                is_persistent <- max_run >= minConsecutive
            }

            class_by_host[i, h] <- if (is_persistent) {
                "persistent"
            } else {
                "transient"
            }
        }
    }

    ## Flag transient iSNVs in the samples of the host where they are
    ## transient; a site can be transient in one host and persistent in
    ## another
    qc <- assay(whe, "qcPass")
    for (h in unique_hosts) {
        transient_rows <- which(class_by_host[, h] == "transient")
        if (!length(transient_rows)) next
        host_cols <- which(hosts == h & !is.na(hosts))
        newly_flagged <- qc[transient_rows, host_cols, drop = FALSE]
        qc[transient_rows, host_cols] <- FALSE
        total_flagged <- total_flagged + sum(newly_flagged)
    }
    assay(whe, "qcPass") <- qc

    ## Site-level summary across hosts
    temporal_class <- apply(class_by_host, 1L, function(classes) {
        classes <- classes[!is.na(classes)]
        if (!length(classes)) {
            NA_character_
        } else if (any(classes == "persistent")) {
            "persistent"
        } else {
            "transient"
        }
    })

    ## Add temporal_class to rowRanges
    rr <- SummarizedExperiment::rowRanges(whe)
    S4Vectors::mcols(rr)$temporal_class <- temporal_class
    SummarizedExperiment::rowRanges(whe) <- rr

    ## Log
    whe <- .addQcStep(whe, "temporal_consistency", "minTimepoints",
                      minTimepoints, total_flagged,
                      sum(qc, na.rm = TRUE))

    whe
}
