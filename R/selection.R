# R/selection.R
# Within-host dN/dS and selection coefficients

#' @include methods.R annotate-codon.R
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
#' Estimates the ratio of nonsynonymous to synonymous divergence among
#' the iSNVs of each sample with the counting method of Nei and Gojobori
#' (1986). The numbers of synonymous (\eqn{S}) and nonsynonymous
#' (\eqn{N}) sites are counted over the complete codons of every CDS in
#' \code{gff}, using the reference sequence in \code{refFasta} and the
#' standard genetic code. With \eqn{S_d} synonymous and \eqn{N_d}
#' nonsynonymous iSNVs present in a sample,
#' \deqn{p_S = S_d / S, \quad p_N = N_d / N,}
#' each proportion is corrected for multiple hits with Jukes and Cantor
#' (1969), \eqn{d = -\frac{3}{4}\log(1 - \frac{4}{3}p)}, and
#' \eqn{d_N/d_S} is their ratio. Under neutrality the expected ratio is
#' 1; values below 1 suggest purifying selection and values above 1
#' diversifying selection, but within-host ratios computed from few
#' iSNVs are highly variable and should be read with their counts.
#'
#' Each iSNV present in a sample (non-missing, non-zero frequency, and
#' passing QC when \code{usePassedOnly = TRUE}) is counted once,
#' irrespective of its frequency. Amino acid classes come from
#' \code{\link{annotateCodonChange}}, which should be run with the same
#' \code{gff} and \code{refFasta}. Variants outside a CDS, or without an
#' amino acid annotation, are not counted. Nonsense changes count as
#' nonsynonymous. Codons shared by overlapping CDS features of the same
#' gene are counted once.
#'
#' @param whe A \code{\link{WithinHostExperiment}} annotated with
#'   \code{\link{annotateCodonChange}}.
#' @param gff Character scalar. Path to the GFF3 file with CDS features
#'   used for the annotation.
#' @param refFasta Character scalar. Path to the reference genome FASTA
#'   used for the annotation.
#' @param gffFeatureCol Character scalar. Column name in
#'   \code{mcols(rowRanges)} containing the gene annotation
#'   (default: \code{"GFF_FEATURE"}).
#' @param refAACol Character scalar. Column name in
#'   \code{mcols(rowRanges)} for the reference amino acid
#'   (default: \code{"REF_AA"}).
#' @param altAACol Character scalar. Column name in
#'   \code{mcols(rowRanges)} for the alternative amino acid
#'   (default: \code{"ALT_AA"}).
#' @param usePassedOnly Logical. If \code{TRUE} (default), only use
#'   variants where \code{qcPass == TRUE}.
#' @param nameCol Character scalar. GFF3 attribute holding the gene name
#'   (default: \code{"gene"}), as in \code{\link{annotateCodonChange}}.
#'
#' @return A \code{DataFrame} with one row per sample and gene carrying
#'   at least one counted iSNV (one row with \code{gene = NA} when a
#'   sample has none):
#'   \describe{
#'     \item{sample_id}{Sample identifier.}
#'     \item{nS, nN}{Synonymous and nonsynonymous iSNVs in the sample.}
#'     \item{S_sites, N_sites}{Synonymous and nonsynonymous sites over
#'       all CDS features.}
#'     \item{pS, pN}{Proportions \eqn{S_d/S} and \eqn{N_d/N}.}
#'     \item{dNdS}{Jukes-Cantor corrected \eqn{d_N/d_S} (\code{NA} when
#'       \eqn{d_S} is zero or a proportion reaches 0.75).}
#'     \item{gene, gene_nS, gene_nN, gene_S_sites, gene_N_sites,
#'       gene_dNdS}{The same quantities for each gene.}
#'   }
#'
#' @references
#' Nei M, Gojobori T (1986). Simple methods for estimating the
#' numbers of synonymous and nonsynonymous nucleotide substitutions.
#' \emph{Mol Biol Evol} 3:418-426.
#' \doi{10.1093/oxfordjournals.molbev.a040410}
#'
#' Jukes TH, Cantor CR (1969). Evolution of protein molecules. In Munro
#' HN (ed.), \emph{Mammalian Protein Metabolism}, pp. 21-132. Academic
#' Press, New York.
#'
#' @seealso \code{\link{annotateCodonChange}}
#'
#' @export
#' @examples
#' ref <- tempfile(fileext = ".fa")
#' writeLines(c(">seg1", "ATGCTGAAAGGGTAA"), ref)
#' gff <- tempfile(fileext = ".gff3")
#' writeLines(c("##gff-version 3",
#'     "seg1\t.\tCDS\t1\t15\t.\t+\t0\tgene=geneA"), gff)
#' gr <- GenomicRanges::GRanges("seg1",
#'     IRanges::IRanges(c(1, 6), width = 1))
#' S4Vectors::mcols(gr)$ref <- c("A", "G")
#' S4Vectors::mcols(gr)$alt <- c("G", "A")
#' whe <- WithinHostExperiment(
#'     assays = list(altFreq = matrix(c(0.1, 0.2), ncol = 1)),
#'     rowRanges = gr,
#'     colData = S4Vectors::DataFrame(sample_id = "S1"))
#' whe <- annotateCodonChange(whe, gff, ref)
#' dndsWithinHost(whe, gff = gff, refFasta = ref)
dndsWithinHost <- function(whe, gff, refFasta,
                           gffFeatureCol = "GFF_FEATURE",
                           refAACol = "REF_AA",
                           altAACol = "ALT_AA",
                           usePassedOnly = TRUE,
                           nameCol = "gene") {
    if (!is(whe, "WithinHostExperiment"))
        stop("'whe' must be a WithinHostExperiment object.")
    if (missing(gff) || missing(refFasta))
        stop("'gff' and 'refFasta' are required: dN/dS is normalised by ",
             "the synonymous and nonsynonymous sites of the coding ",
             "sequences.")
    if (!is.character(gff) || length(gff) != 1L || !file.exists(gff))
        stop("GFF file not found: ", gff)
    if (!is.character(refFasta) || length(refFasta) != 1L ||
        !file.exists(refFasta))
        stop("Reference FASTA not found: ", refFasta)

    rr <- rowRanges(whe)
    mc <- mcols(rr)
    for (column in c(refAACol, altAACol, gffFeatureCol)) {
        if (!column %in% colnames(mc))
            stop("Column '", column, "' not found in mcols(rowRanges). ",
                 "Run annotateCodonChange() first. Available columns: ",
                 paste(colnames(mc), collapse = ", "))
    }

    sites <- .cds_site_counts(.read_gff_cds(gff, nameCol),
                              .read_fasta_simple(refFasta))
    if (nrow(sites) == 0L)
        stop("No complete codons found for the CDS features in 'gff' ",
             "on the sequences in 'refFasta'.")
    S_total <- sum(sites$S_sites)
    N_total <- sum(sites$N_sites)

    ref_aa <- as.character(mc[[refAACol]])
    alt_aa <- as.character(mc[[altAACol]])
    gene_vals <- as.character(mc[[gffFeatureCol]])
    coding <- !is.na(ref_aa) & !is.na(alt_aa) & !is.na(gene_vals) &
        gene_vals %in% sites$gene
    is_syn <- coding & ref_aa == alt_aa
    is_non <- coding & ref_aa != alt_aa

    freq_mat <- assay(whe, "altFreq")
    qc_mat <- if (usePassedOnly && "qcPass" %in% assayNames(whe)) {
        assay(whe, "qcPass")
    } else {
        NULL
    }
    sample_ids <- colnames(freq_mat)
    if (is.null(sample_ids))
        sample_ids <- paste0("S", seq_len(ncol(freq_mat)))

    jukes_cantor <- function(p) {
        ifelse(p < 0.75, -0.75 * log(1 - 4 * p / 3), NA_real_)
    }
    ratio <- function(nS, nN, S, N) {
        dS <- jukes_cantor(nS / S)
        dN <- jukes_cantor(nN / N)
        if (is.na(dS) || is.na(dN) || dS == 0) NA_real_ else dN / dS
    }

    results <- lapply(seq_len(ncol(freq_mat)), function(j) {
        freq_j <- freq_mat[, j]
        if (!is.null(qc_mat)) freq_j[!qc_mat[, j]] <- NA
        present <- !is.na(freq_j) & freq_j > 0

        nS <- sum(present & is_syn)
        nN <- sum(present & is_non)
        base_row <- data.frame(
            sample_id = sample_ids[j], nS = nS, nN = nN,
            S_sites = S_total, N_sites = N_total,
            pS = nS / S_total, pN = nN / N_total,
            dNdS = ratio(nS, nN, S_total, N_total),
            stringsAsFactors = FALSE)

        genes <- sort(unique(gene_vals[present & coding]))
        if (length(genes) == 0L) {
            return(cbind(base_row, data.frame(
                gene = NA_character_, gene_nS = NA_integer_,
                gene_nN = NA_integer_, gene_S_sites = NA_real_,
                gene_N_sites = NA_real_, gene_dNdS = NA_real_,
                stringsAsFactors = FALSE)))
        }
        gene_rows <- lapply(genes, function(g) {
            in_gene <- present & coding & gene_vals == g
            gS <- sites$S_sites[match(g, sites$gene)]
            gN <- sites$N_sites[match(g, sites$gene)]
            gnS <- sum(in_gene & is_syn)
            gnN <- sum(in_gene & is_non)
            data.frame(gene = g, gene_nS = gnS, gene_nN = gnN,
                       gene_S_sites = gS, gene_N_sites = gN,
                       gene_dNdS = ratio(gnS, gnN, gS, gN),
                       stringsAsFactors = FALSE)
        })
        cbind(base_row[rep(1L, length(genes)), , drop = FALSE],
              do.call(rbind, gene_rows))
    })

    out <- do.call(rbind, results)
    rownames(out) <- NULL
    DataFrame(out)
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
#' Feder AF et al. (2016). More effective drugs lead to harder
#' selective sweeps in the evolution of drug resistance in HIV-1.
#' \emph{eLife} 5:e10670. \doi{10.7554/eLife.10670}
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
    delta_logit <- diff(logit_p)
    delta_t <- diff(times) / generation_time

    if (any(delta_t <= 0))
        stop("'times' must contain strictly increasing values.")

    delta_logit / delta_t
}
