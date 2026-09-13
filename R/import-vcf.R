# R/import-vcf.R
# Import iSNV data from VCF files

#' @include methods.R
#' @importFrom utils read.delim
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom stats setNames
NULL

# ============================================================
# readWithinHost -- main entry point
# ============================================================

#' Import iSNV data from VCF files
#'
#' Reads one or more VCF files produced by variant callers (iVar,
#' LoFreq, Freebayes, or generic) and constructs a
#' \code{\link{WithinHostExperiment}}.
#'
#' Each VCF file corresponds to one sample. The \code{colData}
#' argument provides sample metadata and must contain a
#' \code{sample_id} column whose length matches the number of files.
#'
#' @param vcfFiles Character vector. Paths to VCF files (one per sample).
#' @param colData \code{DataFrame} or \code{data.frame}. Sample
#'   metadata. Must contain \code{sample_id}. May also contain
#'   \code{host_id}, \code{role}, \code{pair_id},
#'   \code{replicate_group}, \code{replicate_id}, etc.
#' @param caller Character scalar. Which caller produced the VCF:
#'   \code{"ivar"}, \code{"lofreq"}, \code{"auto"} (default),
#'   or \code{"generic"}.
#' @param genome Character scalar (optional). Path to reference FASTA
#'   (stored in metadata, not currently used for annotation).
#'
#' @return A \code{\link{WithinHostExperiment}}.
#'
#' @details
#' When technical replicates are not available, Roder et al. (2023,
#' \emph{mBio} 14:e01046-23) recommend using a combination of
#' multiple variant callers with stringent cutoffs to reduce false
#' positives. To implement multi-caller consensus, import VCFs from
#' each caller separately and intersect the resulting objects using
#' standard \code{GRanges} overlap operations.
#'
#' @export
#' @examples
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' whe
readWithinHost <- function(vcfFiles, colData, caller = "auto",
                           genome = NULL) {
    ## ---- Validate inputs ----
    if (!is.character(vcfFiles) || length(vcfFiles) < 1L)
        stop("'vcfFiles' must be a non-empty character vector.")
    missing_f <- vcfFiles[!file.exists(vcfFiles)]
    if (length(missing_f) > 0L) {
        stop("VCF file(s) not found: ", paste(missing_f, collapse = ", "))
    }

    if (is.data.frame(colData) && !is(colData, "DataFrame")) {
        colData <- DataFrame(colData)
    }
    if (!is(colData, "DataFrame"))
        stop("'colData' must be a DataFrame (S4Vectors::DataFrame).")

    if (!"sample_id" %in% colnames(colData)) {
        stop("'colData' must contain a 'sample_id' column. ",
             "Found columns: ",
             paste(colnames(colData), collapse = ", "))
    }
    if (nrow(colData) != length(vcfFiles)) {
        stop("Length of 'vcfFiles' (", length(vcfFiles),
             ") must match nrow(colData) (", nrow(colData), ")")
    }

    ## ---- Detect caller ----
    valid_callers <- c("ivar", "lofreq", "freebayes", "auto", "generic")
    caller <- match.arg(caller, valid_callers)
    if (caller == "auto") {
        caller <- .detectCaller(vcfFiles[1L])
        message("Auto-detected caller: ", caller)
    }

    ## ---- Parse each VCF ----
    parsed <- lapply(vcfFiles, function(f) {
        .parseVcf(f, caller)
    })

    ## ---- Merge and build WHE ----
    meta <- list()
    if (!is.null(genome)) meta[["reference_genome"]] <- genome

    .mergeToWHE(parsed, colData, callerName = caller, metadata = meta)
}


# ============================================================
# Internal: caller detection
# ============================================================

#' @keywords internal
.detectCaller <- function(vcfFile) {
    con <- file(vcfFile, "r")
    on.exit(close(con))
    lines <- character()
    repeat {
        line <- readLines(con, n = 1L)
        if (length(line) == 0L || !startsWith(line, "##")) break
        lines <- c(lines, line)
    }
    header_text <- paste(lines, collapse = "\n")

    if (grepl("iVar", header_text, ignore.case = TRUE)) return("ivar")
    if (grepl("lofreq", header_text, ignore.case = TRUE)) return("lofreq")
    if (grepl("freebayes", header_text, ignore.case = TRUE)) return("freebayes")
    "generic"
}


# ============================================================
# Internal: VCF parsing (unified entry)
# ============================================================

#' @keywords internal
.parseVcf <- function(vcfFile, caller) {
    ## Read data lines (skip ## and # header)
    all_lines <- readLines(vcfFile)
    data_lines <- all_lines[!startsWith(all_lines, "#")]

    if (length(data_lines) == 0L) {
        return(list(
            chrom = character(), pos = integer(), ref = character(),
            alt = character(), altFreq = numeric(), totalDepth = integer(),
            refCount = integer(), altCount = integer()
        ))
    }

    ## Parse tab-separated fields
    fields <- strsplit(data_lines, "\t")
    n_cols <- length(fields[[1L]])
    ## Validate consistent column count
    col_counts <- vapply(fields, length, integer(1))
    if (any(col_counts != n_cols)) {
        warning("Inconsistent column counts in VCF: ",
                "expected ", n_cols, " but found rows with ",
                paste(unique(col_counts[col_counts != n_cols]),
                      collapse = ", "),
                " columns. Using first row format.", call. = FALSE)
    }
    chrom  <- vapply(fields, `[`, character(1), 1L)
    pos    <- as.integer(vapply(fields, `[`, character(1), 2L))
    ref    <- vapply(fields, `[`, character(1), 4L)
    alt    <- vapply(fields, `[`, character(1), 5L)
    info   <- vapply(fields, `[`, character(1), 8L)

    ## Check if FORMAT/SAMPLE columns exist (LoFreq VCFs have only 8 cols)
    has_format <- n_cols >= 10L

    ## Helper to extract INFO field values
    .get_info_val <- function(tag, as_type = "character") {
        pattern <- paste0("(?:^|;)", tag, "=([^;]+)")
        vals <- vapply(info, function(inf) {
            m <- regmatches(inf, regexpr(pattern, inf, perl = TRUE))
            if (length(m) == 0L || nchar(m) == 0L) {
                NA_character_
            } else {
                gsub(paste0("^;?", tag, "="), "", m)
            }
        }, character(1), USE.NAMES = FALSE)
        if (as_type == "numeric") return(as.numeric(vals))
        if (as_type == "integer") return(as.integer(vals))
        vals
    }

    if (has_format) {
        fmt  <- vapply(fields, `[`, character(1), 9L)
        samp <- vapply(fields, `[`, character(1), 10L)
        ## Assume consistent FORMAT across rows within a single
        ## caller's output. This is true for iVar, LoFreq, and
        ## Freebayes single-sample VCFs.
        fmt_names  <- strsplit(fmt[1L], ":")[[1L]]
        samp_split <- strsplit(samp, ":")

        .get_fmt_field <- function(field_name) {
            idx <- match(field_name, fmt_names)
            if (is.na(idx)) {
                return(rep(NA_character_, length(samp_split)))
            }
            vapply(samp_split, function(s) {
                if (length(s) >= idx) s[idx] else NA_character_
            }, character(1))
        }
    }

    ## Extract values based on caller
    if (caller == "ivar" && has_format) {
        alt_freq <- as.numeric(.get_fmt_field("ALT_FREQ"))
        alt_dp   <- as.integer(.get_fmt_field("ALT_DP"))
        ref_dp   <- as.integer(.get_fmt_field("REF_DP"))
        total_dp <- alt_dp + ref_dp
        info_dp  <- .get_info_val("DP", "integer")
        total_dp <- ifelse(is.na(total_dp), info_dp, total_dp)
    } else if (caller == "lofreq") {
        ## LoFreq: AF and DP in INFO field (no FORMAT/SAMPLE columns)
        alt_freq <- .get_info_val("AF", "numeric")
        total_dp <- .get_info_val("DP", "integer")
        alt_dp   <- as.integer(round(total_dp * alt_freq))
        ref_dp   <- total_dp - alt_dp
    } else if (caller == "freebayes" && has_format) {
        ## Freebayes: uses RO/AO in FORMAT
        alt_dp   <- as.integer(.get_fmt_field("AO"))
        ref_dp   <- as.integer(.get_fmt_field("RO"))
        total_dp <- as.integer(.get_fmt_field("DP"))
        if (all(is.na(total_dp))) {
            total_dp <- .get_info_val("DP", "integer")
        }
        alt_freq <- ifelse(!is.na(alt_dp) & !is.na(total_dp) &
                               total_dp > 0L,
                           alt_dp / total_dp, NA_real_)
    } else if (has_format) {
        ## Generic with FORMAT: try common field names
        alt_freq_raw <- .get_fmt_field("ALT_FREQ")
        if (all(is.na(alt_freq_raw))) {
            alt_freq_raw <- .get_fmt_field("AF")
        }
        alt_freq <- as.numeric(alt_freq_raw)

        alt_dp_raw <- .get_fmt_field("ALT_DP")
        if (all(is.na(alt_dp_raw))) alt_dp_raw <- .get_fmt_field("AO")
        alt_dp <- as.integer(alt_dp_raw)

        ref_dp_raw <- .get_fmt_field("REF_DP")
        if (all(is.na(ref_dp_raw))) ref_dp_raw <- .get_fmt_field("RO")
        ref_dp <- as.integer(ref_dp_raw)

        total_dp <- as.integer(.get_fmt_field("DP"))
        if (all(is.na(total_dp))) {
            info_dp  <- .get_info_val("DP", "integer")
            total_dp <- ifelse(!is.na(alt_dp) & !is.na(ref_dp),
                               alt_dp + ref_dp, info_dp)
        }
    } else {
        ## No FORMAT columns -- parse INFO only (like LoFreq)
        alt_freq <- .get_info_val("AF", "numeric")
        total_dp <- .get_info_val("DP", "integer")
        alt_dp   <- as.integer(round(total_dp * alt_freq))
        ref_dp   <- total_dp - alt_dp
    }

    list(
        chrom      = chrom,
        pos        = pos,
        ref        = ref,
        alt        = alt,
        altFreq    = alt_freq,
        totalDepth = total_dp,
        refCount   = ref_dp,
        altCount   = alt_dp
    )
}


# ============================================================
# Internal: merge multiple samples -> WithinHostExperiment
# ============================================================

#' @keywords internal
.mergeToWHE <- function(parsed_list, colData, callerName, metadata = list()) {
    sample_ids <- as.character(colData$sample_id)
    n_samples  <- length(parsed_list)

    ## 1. Collect all unique positions across samples
    all_positions <- unique(do.call(rbind, lapply(parsed_list, function(p) {
        if (length(p$chrom) == 0L) {
            return(data.frame(chrom = character(), pos = integer(),
                              ref = character(), alt = character(),
                              stringsAsFactors = FALSE))
        }
        data.frame(chrom = p$chrom, pos = p$pos, ref = p$ref, alt = p$alt,
                   stringsAsFactors = FALSE)
    })))

    if (nrow(all_positions) == 0L) {
        ## No variants found - return empty WHE
        gr <- GRanges()
        empty_num <- matrix(nrow = 0, ncol = n_samples,
                            dimnames = list(NULL, sample_ids))
        empty_int <- matrix(integer(0), nrow = 0, ncol = n_samples,
                            dimnames = list(NULL, sample_ids))
        return(.makeWHE(
            rowRanges = gr,
            assays = list(altFreq = empty_num, totalDepth = empty_int,
                          refCount = empty_int, altCount = empty_int),
            colData = colData, callerName = callerName, metadata = metadata
        ))
    }

    ## Sort by chrom + pos
    all_positions <- all_positions[order(all_positions$chrom,
                                        all_positions$pos), ]
    ## Remove exact duplicates (same chrom+pos+ref+alt)
    dup_key <- paste(all_positions$chrom, all_positions$pos,
                     all_positions$ref, all_positions$alt, sep = ":")
    all_positions <- all_positions[!duplicated(dup_key), ]
    n_sites <- nrow(all_positions)

    ## 2. Build GRanges
    gr <- GRanges(
        seqnames = all_positions$chrom,
        ranges   = IRanges(all_positions$pos, width = 1L)
    )
    S4Vectors::mcols(gr)$ref <- all_positions$ref
    S4Vectors::mcols(gr)$alt <- all_positions$alt

    ## 3. Build assay matrices
    freq_mat  <- matrix(NA_real_, nrow = n_sites, ncol = n_samples,
                        dimnames = list(NULL, sample_ids))
    depth_mat <- matrix(NA_integer_, nrow = n_sites, ncol = n_samples,
                        dimnames = list(NULL, sample_ids))
    ref_mat   <- matrix(NA_integer_, nrow = n_sites, ncol = n_samples,
                        dimnames = list(NULL, sample_ids))
    alt_mat   <- matrix(NA_integer_, nrow = n_sites, ncol = n_samples,
                        dimnames = list(NULL, sample_ids))

    master_key <- paste(all_positions$chrom, all_positions$pos,
                        all_positions$ref, all_positions$alt, sep = ":")

    for (i in seq_len(n_samples)) {
        p <- parsed_list[[i]]
        if (length(p$chrom) == 0L) next
        sample_key <- paste(p$chrom, p$pos, p$ref, p$alt, sep = ":")
        idx   <- match(sample_key, master_key)
        valid <- !is.na(idx)
        freq_mat[idx[valid], i]  <- p$altFreq[valid]
        depth_mat[idx[valid], i] <- p$totalDepth[valid]
        ref_mat[idx[valid], i]   <- p$refCount[valid]
        alt_mat[idx[valid], i]   <- p$altCount[valid]
    }

    ## 4. Construct WHE
    .makeWHE(
        rowRanges = gr,
        assays = list(altFreq    = freq_mat,
                      totalDepth = depth_mat,
                      refCount   = ref_mat,
                      altCount   = alt_mat),
        colData    = colData,
        callerName = callerName,
        metadata   = metadata
    )
}
