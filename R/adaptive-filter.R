# R/adaptive-filter.R
# Technology-aware ISNVFilter presets

#' Create a technology-appropriate ISNVFilter
#'
#' Returns an \code{\link{ISNVFilter}} with default thresholds tuned
#' for the error profile and read characteristics of the specified
#' sequencing technology.
#'
#' @param technology Character scalar. One of \code{"illumina"},
#'   \code{"ont"}, \code{"pacbio_hifi"}, or \code{"generic"}.
#'
#' @return An \code{\link{ISNVFilter}} object.
#'
#' @details
#' Each technology has distinct error modes that motivate different
#' default thresholds:
#'
#' \describe{
#'   \item{\strong{Illumina} (\code{minFreq = 0.03, minDepth = 100})}{
#'     Illumina short reads have low per-base error rates (~0.1--0.5\%)
#'     but systematic context-dependent errors. A 3\% frequency
#'     threshold safely exceeds the noise floor while retaining true
#'     low-frequency variants. Grubaugh et al. (2019) measured iSNVs
#'     above 3\% accurately from at least 1,000 RNA copies sequenced
#'     to at least 400x, so 100x is a permissive floor; low-frequency
#'     false positives inflate diversity estimates markedly (McCrone &
#'     Lauring 2016).
#'   }
#'   \item{\strong{Oxford Nanopore (ONT)} (\code{minFreq = 0.10,
#'     minDepth = 60, maxStrandBias = Inf})}{
#'     ONT reads have higher per-base error rates (~1--5\% with R10
#'     chemistry), requiring a more conservative frequency threshold
#'     of 10\% to avoid false positives; Bull et al. (2020) found that
#'     ONT sequencing does not accurately detect variants at low
#'     read-count frequencies, and reported over 99\% sensitivity and
#'     precision for single-nucleotide variants above about 60-fold
#'     depth, which sets the depth threshold. Strand bias filtering is
#'     disabled (\code{Inf}) rather than given an untested threshold;
#'     supply \code{maxStrandBias} explicitly when the data support one.
#'   }
#'   \item{\strong{PacBio HiFi} (\code{minFreq = 0.05, minDepth = 50})}{
#'     HiFi consensus reads reach about 99.8\% read accuracy via
#'     circular consensus sequencing (Wenger et al. 2019), allowing a
#'     frequency threshold below ONT but above Illumina. Depth of 50x
#'     reflects the lower throughput but high per-read accuracy.
#'   }
#'   \item{\strong{Generic}}{
#'     Falls back to the standard \code{\link{ISNVFilter}} defaults
#'     (\code{minFreq = 0.03, minDepth = 100}), suitable for Illumina
#'     or when the technology is unknown.
#'   }
#' }
#'
#' @references
#' McCrone JT, Lauring AS (2016). Measurements of intrahost viral
#' diversity are extremely sensitive to systematic errors in variant
#' calling. \emph{J Virol} 90:6884-6895.
#' \doi{10.1128/JVI.00667-16}
#'
#' Grubaugh ND et al. (2019). An amplicon-based sequencing framework
#' for accurately measuring intrahost virus diversity using PrimalSeq
#' and iVar. \emph{Genome Biology} 20:8.
#' \doi{10.1186/s13059-018-1618-7}
#'
#' Bull RA et al. (2020). Analytical validity of nanopore sequencing
#' for rapid SARS-CoV-2 genome analysis. \emph{Nature Communications}
#' 11:6272. \doi{10.1038/s41467-020-20075-6}
#'
#' Wenger AM et al. (2019). Accurate circular consensus long-read
#' sequencing improves variant detection and assembly of a human
#' genome. \emph{Nature Biotechnology} 37:1155--1162.
#' \doi{10.1038/s41587-019-0217-9}
#'
#' @seealso \code{\link{ISNVFilter}} for the underlying constructor,
#'   \code{\link{flagISNV}} for applying filters.
#'
#' @export
#' @examples
#' # Illumina defaults
#' adaptiveFilter("illumina")
#'
#' # ONT with strand bias disabled
#' filt <- adaptiveFilter("ont")
#' slot(filt, "maxStrandBias")  # Inf
#'
#' # PacBio HiFi
#' adaptiveFilter("pacbio_hifi")
#'
#' # Generic (same as ISNVFilter())
#' adaptiveFilter("generic")
adaptiveFilter <- function(technology = c("illumina", "ont",
                                          "pacbio_hifi", "generic")) {
    technology <- match.arg(technology)

    switch(technology,
        illumina = ISNVFilter(
            minFreq       = 0.03,
            minDepth      = 100L
        ),
        ont = ISNVFilter(
            minFreq       = 0.10,
            minDepth      = 60L,
            maxStrandBias = Inf
        ),
        pacbio_hifi = ISNVFilter(
            minFreq       = 0.05,
            minDepth      = 50L
        ),
        generic = ISNVFilter()
    )
}
