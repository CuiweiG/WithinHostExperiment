# R/methods.R
# Constructor functions, accessor methods, show methods, internal helpers

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom methods setMethod setAs new callNextMethod is show
#' @importFrom BiocGenerics nrow ncol
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom VariantAnnotation VRanges
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames
#'   colData rowRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
NULL

# ============================================================
# ISNVFilter constructor + show
# ============================================================

#' Create an ISNVFilter object
#'
#' Constructs a quality control filter configuration for use with
#' \code{\link{flagISNV}}.
#'
#' @param minDepth Integer. Minimum read depth (default: 100).
#' @param minFreq Numeric. Minimum alt allele frequency (default: 0.03).
#' @param maxFreq Numeric. Maximum alt allele frequency (default: 0.50).
#' @param minAltReads Integer. Minimum alt allele reads (default: 10).
#' @param maxStrandBias Numeric. Max strand bias fold-diff (default: 10).
#' @param minBaseQual Integer. Minimum base quality (default: 20).
#' @param minMapQual Integer. Minimum mapping quality (default: 20).
#' @param replicateConc Logical. Require replicate concordance
#'   (default: FALSE).
#'
#' @return An \code{\link{ISNVFilter}} object.
#' @export
#' @examples
#' ISNVFilter()
#' ISNVFilter(minDepth = 500L, minFreq = 0.05)
ISNVFilter <- function(minDepth = 100L, minFreq = 0.03, maxFreq = 0.50,
                       minAltReads = 10L, maxStrandBias = 10.0,
                       minBaseQual = 20L, minMapQual = 20L,
                       replicateConc = FALSE) {
    new("ISNVFilter",
        minDepth      = as.integer(minDepth),
        minFreq       = as.numeric(minFreq),
        maxFreq       = as.numeric(maxFreq),
        minAltReads   = as.integer(minAltReads),
        maxStrandBias = as.numeric(maxStrandBias),
        minBaseQual   = as.integer(minBaseQual),
        minMapQual    = as.integer(minMapQual),
        replicateConc = as.logical(replicateConc))
}

#' @rdname ISNVFilter-class
#' @param object An \code{ISNVFilter} object.
#' @export
setMethod("show", "ISNVFilter", function(object) {
    cat("ISNVFilter object\n")
    cat("  minDepth:     ", slot(object, "minDepth"), "\n")
    cat("  freq range:   [", slot(object, "minFreq"), ", ",
        slot(object, "maxFreq"), "]\n", sep = "")
    cat("  minAltReads:  ", slot(object, "minAltReads"), "\n")
    cat("  maxStrandBias:", slot(object, "maxStrandBias"), "\n")
    cat("  minBaseQual:  ", slot(object, "minBaseQual"), "\n")
    cat("  minMapQual:   ", slot(object, "minMapQual"), "\n")
    cat("  replicateConc:", slot(object, "replicateConc"), "\n")
})


# ============================================================
# WithinHostExperiment constructor
# ============================================================

#' Create a WithinHostExperiment from components
#'
#' Constructs a \code{\link{WithinHostExperiment}} from assay matrices,
#' genomic ranges, and sample metadata.
#'
#' @param assays Named list of matrices. Must include at least
#'   \code{altFreq}. If \code{qcPass} is not provided, it is
#'   automatically created as an all-TRUE matrix.
#' @param rowRanges \code{GRanges}. Genomic coordinates for each
#'   variant site.
#' @param colData \code{DataFrame} or \code{data.frame}. Sample
#'   metadata. Must contain a \code{sample_id} column.
#' @param metadata Named list (optional). Additional metadata such as
#'   \code{reference_genome}, \code{genome_length}, \code{calling_info}.
#'
#' @return A \code{\link{WithinHostExperiment}} object.
#'
#' @export
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges("seg1", IRanges::IRanges(c(100, 200, 300), width = 1))
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = matrix(c(0.1, 0.2, 0.3), ncol = 1)),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(sample_id = "S1"))
#' whe
WithinHostExperiment <- function(assays, rowRanges, colData,
                                 metadata = list()) {
    ## Coerce colData
    if (is.data.frame(colData) && !is(colData, "DataFrame")) {
        colData <- DataFrame(colData)
    }

    ## Validate sample_id
    if (!"sample_id" %in% colnames(colData)) {
        stop("'colData' must contain a 'sample_id' column. ",
             "Found columns: ",
             paste(colnames(colData), collapse = ", "))
    }

    n_sites   <- nrow(assays[[1]])
    n_samples <- ncol(assays[[1]])
    sample_ids <- as.character(colData$sample_id)

    ## Set column names on assay matrices from sample_id
    for (nm in names(assays)) {
        if (is.null(colnames(assays[[nm]]))) {
            colnames(assays[[nm]]) <- sample_ids
        }
    }

    ## Auto-create qcPass if missing
    if (!"qcPass" %in% names(assays)) {
        assays[["qcPass"]] <- matrix(TRUE,
            nrow = n_sites, ncol = n_samples,
            dimnames = list(NULL, sample_ids))
    }

    ## Build SE
    se <- SummarizedExperiment(
        assays    = assays,
        rowRanges = rowRanges,
        colData   = colData
    )
    S4Vectors::metadata(se) <- metadata

    ## Construct WHE
    new("WithinHostExperiment",
        se,
        qcLog = DataFrame(
            step      = character(),
            parameter = character(),
            value     = character(),
            n_flagged = integer(),
            n_passed  = integer(),
            timestamp = character()
        ),
        qcFilters = list()
    )
}


# ============================================================
# WithinHostExperiment accessor methods
# ============================================================

#' @rdname WHE-accessors
#' @export
setMethod("qcLog", "WithinHostExperiment", function(x) x@qcLog)

#' @rdname WHE-accessors
#' @export
setMethod("qcFilters", "WithinHostExperiment", function(x) x@qcFilters)

#' @rdname WHE-accessors
#' @export
setReplaceMethod("qcLog", "WithinHostExperiment", function(x, value) {
    x@qcLog <- value
    x
})

#' @rdname WHE-accessors
#' @export
setReplaceMethod("qcFilters", "WithinHostExperiment", function(x, value) {
    x@qcFilters <- value
    x
})

#' @rdname WHE-accessors
#' @export
setMethod("passedISNV", "WithinHostExperiment", function(x) {
    if (!"qcPass" %in% assayNames(x)) {
        return(x)
    }
    qc_mat <- assay(x, "qcPass")
    ## Keep rows where at least one sample passed
    keep <- rowSums(qc_mat, na.rm = TRUE) > 0L
    x[keep, ]
})

#' @rdname WHE-accessors
#' @export
setMethod("transmissionPairs", "WithinHostExperiment", function(x) {
    cd <- colData(x)

    ## Check required columns exist
    if (!"role" %in% colnames(cd) || !"pair_id" %in% colnames(cd)) {
        return(DataFrame(donor     = character(),
                         recipient = character(),
                         pair_id   = character()))
    }

    ## Filter to donor/recipient rows with non-NA pair_id
    has_role <- !is.na(cd$role) & cd$role %in% c("donor", "recipient")
    has_pair <- !is.na(cd$pair_id)
    relevant <- cd[has_role & has_pair, , drop = FALSE]

    if (nrow(relevant) == 0L) {
        return(DataFrame(donor     = character(),
                         recipient = character(),
                         pair_id   = character()))
    }

    ## Group by pair_id
    pair_ids <- unique(as.character(relevant$pair_id))
    result <- lapply(pair_ids, function(pid) {
        sub       <- relevant[relevant$pair_id == pid, , drop = FALSE]
        donors    <- as.character(sub$sample_id[sub$role == "donor"])
        recipients <- as.character(sub$sample_id[sub$role == "recipient"])
        if (length(donors) == 0L || length(recipients) == 0L) {
            return(NULL)
        }
        expand.grid(donor     = donors,
                    recipient = recipients,
                    pair_id   = pid,
                    stringsAsFactors = FALSE)
    })
    result <- do.call(rbind, Filter(Negate(is.null), result))
    if (is.null(result) || nrow(result) == 0L) {
        return(DataFrame(donor     = character(),
                         recipient = character(),
                         pair_id   = character()))
    }
    DataFrame(result)
})


# ============================================================
# WithinHostExperiment show method
# ============================================================

#' @rdname WithinHostExperiment-class
#' @param object A \code{WithinHostExperiment} object.
#' @export
setMethod("show", "WithinHostExperiment", function(object) {
    cat("WithinHostExperiment with",
        nrow(object), "variant sites and",
        ncol(object), "samples\n")
    anames <- assayNames(object)
    if (length(anames) > 0L) {
        cat("  assays(", length(anames), "): ",
            paste(anames, collapse = ", "), "\n", sep = "")
    }
    cat("  qcLog:", nrow(qcLog(object)), "filtering step(s)\n")
    pairs <- transmissionPairs(object)
    cat("  transmission pairs:", nrow(pairs), "\n")
})


# ============================================================
# Internal helper functions
# ============================================================

#' Record a QC step in the audit log
#'
#' @param whe A \code{WithinHostExperiment} object.
#' @param step Character. Filter step name.
#' @param parameter Character. Parameter name.
#' @param value Character. Parameter value (as string).
#' @param n_flagged Integer. Number of entries newly flagged.
#' @param n_passed Integer. Number of entries still passing.
#' @return Updated \code{WithinHostExperiment}.
#' @keywords internal
#' @noRd
.addQcStep <- function(whe, step, parameter, value, n_flagged, n_passed) {
    entry <- DataFrame(
        step      = step,
        parameter = parameter,
        value     = as.character(value),
        n_flagged = as.integer(n_flagged),
        n_passed  = as.integer(n_passed),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    )
    qcLog(whe) <- rbind(qcLog(whe), entry)
    whe
}

#' Internal constructor from pre-validated components
#'
#' @param rowRanges \code{GRanges}. Variant positions.
#' @param assays Named list of matrices.
#' @param colData \code{DataFrame}. Sample metadata.
#' @param callerName Character. Caller name for metadata.
#' @param metadata Named list. Additional metadata.
#' @return A \code{WithinHostExperiment}.
#' @keywords internal
#' @noRd
.makeWHE <- function(rowRanges, assays, colData,
                     callerName = NA_character_, metadata = list()) {
    sample_ids <- as.character(colData$sample_id)
    n_sites    <- length(rowRanges)
    n_samples  <- length(sample_ids)

    ## Ensure colnames on all assays
    for (nm in names(assays)) {
        if (is.null(colnames(assays[[nm]]))) {
            colnames(assays[[nm]]) <- sample_ids
        }
    }

    ## Auto-create qcPass
    if (!"qcPass" %in% names(assays)) {
        assays[["qcPass"]] <- matrix(TRUE,
            nrow = n_sites, ncol = n_samples,
            dimnames = list(NULL, sample_ids))
    }

    ## Build SE
    se <- SummarizedExperiment(
        assays    = assays,
        rowRanges = rowRanges,
        colData   = colData
    )

    ## Store caller in metadata
    if (!is.na(callerName)) {
        metadata[["calling_info"]] <- DataFrame(
            caller  = callerName,
            version = NA_character_
        )
    }
    S4Vectors::metadata(se) <- metadata

    new("WithinHostExperiment",
        se,
        qcLog = DataFrame(
            step      = character(),
            parameter = character(),
            value     = character(),
            n_flagged = integer(),
            n_passed  = integer(),
            timestamp = character()
        ),
        qcFilters = list()
    )
}


# ============================================================
# Subsetting -- preserve custom slots
# ============================================================

#' Subset a WithinHostExperiment
#'
#' Subsetting preserves the \code{qcLog} and \code{qcFilters} slots.
#' These slots are experiment-level metadata (not row/column-dependent)
#' and are carried forward unchanged.
#'
#' @param x A \code{\link{WithinHostExperiment}}.
#' @param i Row indices (variant sites).
#' @param j Column indices (samples).
#' @param ... Additional arguments passed to the parent method.
#' @param drop Logical (ignored for \code{SummarizedExperiment}).
#'
#' @return A \code{\link{WithinHostExperiment}}.
#'
#' @rdname WithinHostExperiment-class
#' @export
setMethod("[", "WithinHostExperiment",
    function(x, i, j, ..., drop = TRUE) {
        result <- callNextMethod()
        new("WithinHostExperiment", result,
            qcLog = qcLog(x),
            qcFilters = qcFilters(x))
    }
)

# ============================================================
# Combine WithinHostExperiment objects
# ============================================================

#' Combine WithinHostExperiment objects by rows or columns
#'
#' \code{combineRows} combines objects vertically (adding variant
#' sites); all objects must have the same samples (columns).
#' \code{combineCols} combines objects horizontally (adding
#' samples); all objects must have the same variant sites (rows).
#' QC logs and filters from all objects are merged.
#'
#' @param ... \code{\link{WithinHostExperiment}} objects to combine.
#'
#' @return A \code{\link{WithinHostExperiment}}.
#'
#' @export
#' @rdname WHE-combine
#' @examples
#' data(example_whe)
#' whe1 <- example_whe[1:10, ]
#' whe2 <- example_whe[11:nrow(example_whe), ]
#' merged <- combineRows(whe1, whe2)
#' nrow(merged) == nrow(example_whe)
combineRows <- function(...) {
    args <- list(...)
    if (!all(vapply(args, is, logical(1), "WithinHostExperiment")))
        stop("All arguments must be WithinHostExperiment objects.")
    ## Merge qcLogs and qcFilters
    all_logs <- .mergeQcLogs(args)
    all_filters <- unlist(lapply(args, qcFilters),
                          recursive = FALSE)

    ## Concatenate rowRanges
    all_rr <- do.call(c, lapply(args, rowRanges))

    ## Get common colData from first object
    cd <- colData(args[[1L]])

    ## Concatenate assay matrices by rows
    anames <- assayNames(args[[1L]])
    merged_assays <- lapply(anames, function(nm) {
        mats <- lapply(args, function(x) assay(x, nm))
        do.call(base::rbind, mats)
    })
    names(merged_assays) <- anames

    ## Build merged WHE
    se <- SummarizedExperiment(
        assays = merged_assays, rowRanges = all_rr, colData = cd)
    S4Vectors::metadata(se) <- S4Vectors::metadata(args[[1L]])
    new("WithinHostExperiment", se,
        qcLog = all_logs, qcFilters = all_filters)
}

#' @export
#' @rdname WHE-combine
#' @examples
#' data(example_whe)
#' whe_d <- example_whe[, 1]
#' whe_r <- example_whe[, 2]
#' merged <- combineCols(whe_d, whe_r)
#' ncol(merged)
combineCols <- function(...) {
    args <- list(...)
    if (!all(vapply(args, is, logical(1), "WithinHostExperiment")))
        stop("All arguments must be WithinHostExperiment objects.")
    all_logs <- .mergeQcLogs(args)
    all_filters <- unlist(lapply(args, qcFilters),
                          recursive = FALSE)

    ## Common rowRanges from first object -- strip all names
    ## to avoid SummarizedExperiment cbind validation conflicts
    rr <- rowRanges(args[[1L]])
    names(rr) <- NULL
    rownames(rr) <- NULL
    if (ncol(S4Vectors::mcols(rr)) > 0L) {
        rownames(S4Vectors::mcols(rr)) <- NULL
    }

    ## Concatenate colData
    all_cd <- do.call(
        function(...) base::rbind(...),
        lapply(args, function(x) as.data.frame(colData(x))))
    all_cd <- DataFrame(all_cd)
    rownames(all_cd) <- NULL

    ## Concatenate assay matrices by columns
    anames <- assayNames(args[[1L]])
    merged_assays <- lapply(anames, function(nm) {
        mats <- lapply(args, function(x) {
            m <- assay(x, nm)
            rownames(m) <- NULL
            m
        })
        do.call(base::cbind, mats)
    })
    names(merged_assays) <- anames

    se <- SummarizedExperiment(
        assays = merged_assays, rowRanges = rr, colData = all_cd)
    S4Vectors::metadata(se) <- S4Vectors::metadata(args[[1L]])
    new("WithinHostExperiment", se,
        qcLog = all_logs, qcFilters = all_filters)
}

#' @keywords internal
#' @noRd
.mergeQcLogs <- function(whe_list) {
    logs_list <- lapply(whe_list, qcLog)
    ## Convert to data.frame for safe rbind, then back
    dfs <- lapply(logs_list, as.data.frame)
    merged <- do.call(base::rbind, dfs)
    if (nrow(merged) == 0L) {
        return(DataFrame(
            step = character(), parameter = character(),
            value = character(), n_flagged = integer(),
            n_passed = integer(), timestamp = character()))
    }
    DataFrame(merged)
}

# ============================================================
# Coercion methods -- interoperability
# ============================================================

#' Coerce WithinHostExperiment to data.frame or VRanges
#'
#' @name coerce-WHE
#' @rdname coerce-WHE
#' @aliases coerce,WithinHostExperiment,data.frame-method
#'   coerce,WithinHostExperiment,VRanges-method
#'   coerce,VRanges,WithinHostExperiment-method
#'
#' @return A \code{data.frame} or \code{\link[VariantAnnotation]{VRanges}}
#'   object containing the variant data from the
#'   \code{WithinHostExperiment} object.
#'
#' @examples
#' data(example_whe)
#' df <- as(example_whe, "data.frame")
#' head(df)
#'
#' vr <- as(example_whe, "VRanges")
#' vr
NULL

setAs("WithinHostExperiment", "data.frame", function(from) {
    rr       <- rowRanges(from)
    freq_mat <- assay(from, "altFreq")
    sample_ids <- colnames(freq_mat)
    if (is.null(sample_ids)) {
        sample_ids <- paste0("S", seq_len(ncol(freq_mat)))
    }

    n_sites   <- nrow(from)
    n_samples <- ncol(from)

    ref_vals <- if ("ref" %in% colnames(mcols(rr))) {
        mcols(rr)$ref
    } else {
        rep(NA_character_, n_sites)
    }
    alt_vals <- if ("alt" %in% colnames(mcols(rr))) {
        mcols(rr)$alt
    } else {
        rep(NA_character_, n_sites)
    }

    result <- data.frame(
        chrom     = rep(as.character(seqnames(rr)), times = n_samples),
        position  = rep(start(rr), times = n_samples),
        ref       = rep(ref_vals, times = n_samples),
        alt       = rep(alt_vals, times = n_samples),
        sample_id = rep(sample_ids, each = n_sites),
        altFreq   = as.numeric(freq_mat),
        stringsAsFactors = FALSE
    )

    for (aname in setdiff(assayNames(from), c("altFreq", "qcPass"))) {
        result[[aname]] <- as.numeric(assay(from, aname))
    }

    if ("qcPass" %in% assayNames(from)) {
        result$qcPass <- as.logical(assay(from, "qcPass"))
    }

    result <- result[!is.na(result$altFreq), , drop = FALSE]
    rownames(result) <- NULL
    result
})

setAs("WithinHostExperiment", "VRanges", function(from) {
    rr       <- rowRanges(from)
    freq_mat <- assay(from, "altFreq")
    sample_ids <- colnames(freq_mat)
    if (is.null(sample_ids)) {
        sample_ids <- paste0("S", seq_len(ncol(freq_mat)))
    }

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

    depth_mat <- if ("totalDepth" %in% assayNames(from)) {
        assay(from, "totalDepth")
    } else {
        matrix(NA_integer_, nrow = nrow(from), ncol = ncol(from))
    }
    alt_mat <- if ("altCount" %in% assayNames(from)) {
        assay(from, "altCount")
    } else {
        matrix(NA_integer_, nrow = nrow(from), ncol = ncol(from))
    }
    ref_mat <- if ("refCount" %in% assayNames(from)) {
        assay(from, "refCount")
    } else {
        matrix(NA_integer_, nrow = nrow(from), ncol = ncol(from))
    }

    n_sites   <- nrow(from)
    n_samples <- ncol(from)

    all_seqnames    <- rep(as.character(seqnames(rr)), times = n_samples)
    all_start       <- rep(start(rr), times = n_samples)
    all_ref         <- rep(ref_col, times = n_samples)
    all_alt         <- rep(alt_col, times = n_samples)
    all_totalDepth  <- as.integer(depth_mat)
    all_altDepth    <- as.integer(alt_mat)
    all_refDepth    <- as.integer(ref_mat)
    all_sampleNames <- rep(sample_ids, each = n_sites)
    all_freq        <- as.numeric(freq_mat)

    keep <- !is.na(all_freq)

    VRanges(
        seqnames    = all_seqnames[keep],
        ranges      = IRanges(all_start[keep], width = 1L),
        ref         = all_ref[keep],
        alt         = all_alt[keep],
        totalDepth  = all_totalDepth[keep],
        altDepth    = all_altDepth[keep],
        refDepth    = all_refDepth[keep],
        sampleNames = all_sampleNames[keep]
    )
})


# ============================================================
# VRanges -> WithinHostExperiment coercion
# ============================================================

setAs("VRanges", "WithinHostExperiment", function(from) {
    ## Extract components
    sample_names <- as.character(
        VariantAnnotation::sampleNames(from))
    unique_samples <- unique(sample_names)
    n_samples <- length(unique_samples)

    ## Build site key for grouping
    chrom_vec <- as.character(GenomeInfoDb::seqnames(from))
    start_vec <- GenomicRanges::start(from)
    ref_vec <- as.character(VariantAnnotation::ref(from))
    alt_vec <- as.character(VariantAnnotation::alt(from))
    site_key <- paste(chrom_vec, start_vec, ref_vec, alt_vec,
                      sep = ":")
    unique_sites <- unique(site_key)
    n_sites <- length(unique_sites)

    ## Build GRanges for unique sites
    site_parts <- strsplit(unique_sites, ":")
    gr <- GRanges(
        seqnames = vapply(site_parts, `[`, character(1), 1L),
        ranges = IRanges(
            as.integer(vapply(site_parts, `[`, character(1), 2L)),
            width = 1L)
    )
    mcols(gr)$ref <- vapply(site_parts, `[`, character(1), 3L)
    mcols(gr)$alt <- vapply(site_parts, `[`, character(1), 4L)

    ## Initialize matrices
    alt_dp_mat <- matrix(NA_integer_, nrow = n_sites,
                         ncol = n_samples,
                         dimnames = list(NULL, unique_samples))
    ref_dp_mat <- matrix(NA_integer_, nrow = n_sites,
                         ncol = n_samples,
                         dimnames = list(NULL, unique_samples))
    total_dp_mat <- matrix(NA_integer_, nrow = n_sites,
                           ncol = n_samples,
                           dimnames = list(NULL, unique_samples))

    ## Fill matrices
    site_idx <- match(site_key, unique_sites)
    sample_idx <- match(sample_names, unique_samples)
    ad <- as.integer(VariantAnnotation::altDepth(from))
    rd <- as.integer(VariantAnnotation::refDepth(from))
    td <- as.integer(VariantAnnotation::totalDepth(from))
    for (k in seq_along(from)) {
        i <- site_idx[k]
        j <- sample_idx[k]
        alt_dp_mat[i, j] <- ad[k]
        ref_dp_mat[i, j] <- rd[k]
        total_dp_mat[i, j] <- td[k]
    }

    ## Compute frequency
    freq_mat <- alt_dp_mat / total_dp_mat

    WithinHostExperiment(
        assays = list(
            altFreq = freq_mat,
            totalDepth = total_dp_mat,
            altCount = alt_dp_mat,
            refCount = ref_dp_mat
        ),
        rowRanges = gr,
        colData = DataFrame(sample_id = unique_samples)
    )
})
