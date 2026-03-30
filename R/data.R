# R/data.R
# Documentation for example datasets

#' Example WithinHostExperiment dataset
#'
#' A small synthetic \code{\link{WithinHostExperiment}} with 3 samples
#' (donor, recipient, independent) containing iSNV data on a 1000bp
#' viral genome segment. The donor and recipient form a transmission
#' pair (pair_id = "pair_1").
#'
#' @format A \code{\link{WithinHostExperiment}} with 20 variant sites
#'   and 3 samples. Assays: altFreq, totalDepth, refCount, altCount,
#'   qcPass.
#'
#' @source Synthetically generated for package examples and testing.
#'
#' @docType data
#' @name example_whe
#' @usage data(example_whe)
#'
#' @examples
#' data(example_whe)
#' example_whe
#' transmissionPairs(example_whe)
NULL
