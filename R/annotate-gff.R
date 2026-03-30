# R/annotate-gff.R
# Annotate iSNV sites from GFF3 gene models

#' @include methods.R
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomicRanges GRanges start findOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<- queryHits
#' @importFrom utils read.delim
NULL

#' Annotate variant sites from a GFF3 gene model
#'
#' Reads a GFF3 file and adds gene/feature annotations to
#' \code{mcols(rowRanges)} of a \code{\link{WithinHostExperiment}}.
#' This automates the annotation step required by
#' \code{\link{dndsWithinHost}}, which expects \code{REF_AA} and
#' \code{ALT_AA} columns.
#'
#' The function overlaps variant positions with GFF features of the
#' specified type (default: \code{"gene"}) and adds a
#' \code{GFF_FEATURE} column containing the gene name or ID.
#' For coding-region annotation (amino acid changes), a reference
#' FASTA is additionally required; this function provides the
#' positional overlap only.
#'
#' @param whe A \code{\link{WithinHostExperiment}} object.
#' @param gff Character scalar. Path to a GFF3 file.
#' @param featureType Character scalar. Which GFF feature type to
#'   extract (default: \code{"gene"}). Common choices:
#'   \code{"gene"}, \code{"CDS"}, \code{"mRNA"}.
#' @param nameCol Character scalar. GFF3 attribute to use as the
#'   feature name (default: \code{"Name"}). Falls back to
#'   \code{"gene"}, then \code{"ID"} if not found.
#'
#' @return A \code{\link{WithinHostExperiment}} with a
#'   \code{GFF_FEATURE} column added to \code{mcols(rowRanges)}.
#'
#' @details
#' The GFF3 file is parsed with base R (no external dependency).
#' Only lines matching \code{featureType} in column 3 are used.
#' Attribute parsing follows the GFF3 specification
#' (key=value pairs separated by semicolons).
#'
#' For full codon-aware annotation (predicting amino acid changes),
#' use dedicated tools such as SnpEff or ANNOVAR upstream, then
#' import the annotated VCF via \code{\link{readWithinHost}}.
#'
#' @export
#' @examples
#' ## Create a minimal GFF3 for demonstration
#' gff_file <- tempfile(fileext = ".gff3")
#' writeLines(c(
#'   "##gff-version 3",
#'   "test_segment\t.\tgene\t1\t500\t.\t+\t.\tName=geneA",
#'   "test_segment\t.\tgene\t501\t1000\t.\t+\t.\tName=geneB"
#' ), gff_file)
#' vcf <- system.file("extdata", "test_donor.vcf",
#'     package = "WithinHostExperiment")
#' whe <- readWithinHost(vcf,
#'     colData = S4Vectors::DataFrame(sample_id = "donor"),
#'     caller = "ivar")
#' whe <- annotateFromGFF(whe, gff_file)
#' S4Vectors::mcols(SummarizedExperiment::rowRanges(whe))$GFF_FEATURE
annotateFromGFF <- function(whe, gff,
                            featureType = "gene",
                            nameCol = "Name") {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    if (!file.exists(gff))
        stop("GFF file not found: ", gff)
    if (!is.character(featureType) || length(featureType) != 1L)
        stop("'featureType' must be a single character string.")

    ## Parse GFF3 (base R, no rtracklayer dependency)
    lines <- readLines(gff)
    lines <- lines[!startsWith(lines, "#") & nchar(lines) > 0]
    if (length(lines) == 0L)
        stop("No feature lines found in GFF file.")

    fields <- strsplit(lines, "\t")
    col_counts <- vapply(fields, length, integer(1))
    fields <- fields[col_counts >= 9L]
    if (length(fields) == 0L)
        stop("No valid 9-column GFF lines found.")

    gff_df <- data.frame(
        seqid  = vapply(fields, `[`, character(1), 1L),
        type   = vapply(fields, `[`, character(1), 3L),
        start  = as.integer(vapply(fields, `[`, character(1), 4L)),
        end    = as.integer(vapply(fields, `[`, character(1), 5L)),
        strand = vapply(fields, `[`, character(1), 7L),
        attrs  = vapply(fields, `[`, character(1), 9L),
        stringsAsFactors = FALSE
    )

    ## Filter to requested feature type
    gff_df <- gff_df[gff_df$type == featureType, , drop = FALSE]
    if (nrow(gff_df) == 0L)
        stop("No features of type '", featureType,
             "' found in GFF file. Available types: ",
             paste(unique(vapply(fields, `[`, character(1), 3L)),
                   collapse = ", "))

    ## Extract feature name from attributes
    .parse_attr <- function(attr_str, key) {
        pairs <- strsplit(attr_str, ";")[[1L]]
        for (p in pairs) {
            kv <- strsplit(trimws(p), "=")[[1L]]
            if (length(kv) == 2L && kv[1L] == key) return(kv[2L])
        }
        NA_character_
    }

    ## Try nameCol, then fallbacks
    gff_df$name <- vapply(gff_df$attrs, .parse_attr, character(1),
                          key = nameCol)
    if (all(is.na(gff_df$name)) && nameCol != "gene") {
        gff_df$name <- vapply(gff_df$attrs, .parse_attr, character(1),
                              key = "gene")
    }
    if (all(is.na(gff_df$name))) {
        gff_df$name <- vapply(gff_df$attrs, .parse_attr, character(1),
                              key = "ID")
    }

    ## Build GRanges for GFF features
    gff_gr <- GRanges(
        seqnames = gff_df$seqid,
        ranges   = IRanges(gff_df$start, gff_df$end)
    )

    ## Overlap with variant positions
    rr <- rowRanges(whe)
    hits <- findOverlaps(rr, gff_gr)

    ## Assign feature names
    feature_vec <- rep(NA_character_, length(rr))
    feature_vec[queryHits(hits)] <- gff_df$name[
        S4Vectors::subjectHits(hits)]

    mcols(rr)$GFF_FEATURE <- feature_vec
    ## Update rowRanges in place (preserve class)
    whe_new <- whe
    SummarizedExperiment::rowRanges(whe_new) <- rr
    whe_new
}
