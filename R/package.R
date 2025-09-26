#' admixspatial: Spatial Visualization of ADMIXTURE Results
#'
#' A comprehensive R package for visualizing ADMIXTURE results with spatial context.
#' The package supports both library usage and command-line interface, providing
#' tools for creating publication-ready visualizations of population genetic structure.
#'
#' @details
#' The package provides functionality for:
#' \itemize{
#'   \item Reading and validating ADMIXTURE input files (.Q, .fam, popmap, coordinates)
#'   \item Creating admixture bar plots with customizable color palettes
#'   \item Generating spatial maps with automatic geographic data integration
#'   \item Cross-validation error plots for model selection
#'   \item Command-line interface for batch processing
#' }
#'
#' @section Library Usage:
#' The package can be used programmatically:
#' \code{library(admixspatial)}
#' \code{visualize_admixture(q_file, fam_file, popmap_file, coords_file)}
#'
#' @section CLI Usage:
#' The package includes a command-line script:
#' \code{Rscript inst/scripts/admixture_visualization.R --conf_file config.yaml}
#'
#' @author Sebastian Guzman Rodriguez
#' @name admixspatial-package
#' @aliases admixspatial
#' @keywords package
NULL