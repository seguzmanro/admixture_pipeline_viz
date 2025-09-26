#' Command-Line Interface Functions
#'
#' This module contains functions for handling command-line arguments
#' and executing the pipeline from the command line.

#' Get the directory of the currently running R script
#'
#' @return The directory of the R script
#' @export
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    # When running in package context, return package directory
    return(system.file("", package = "admixspatial"))
  }
  return(dirname(sub("--file=", "", file_arg[1])))
}

#' Main CLI function
#'
#' @export
main <- function() {
  parser <- argparse::ArgumentParser(
    description = "Spatial Visualization of ADMIXTURE Results"
  )

  # Config file
  parser$add_argument("--conf_file", default = NULL,
                      help = "Path to config.yaml file (overrides CLI args)")

  # Input modes
  parser$add_argument("--q_file", default = NULL,
                      help = "Path to single Q file (single mode)")
  parser$add_argument("--input_dir", default = "",
                      help = "Directory with multiple Q files (batch mode)")
  parser$add_argument("--log_dir", default = NULL,
                      help = "Optional: Directory with .out files for CV error plot")

  # Common inputs
  parser$add_argument("--k_value", type = "integer", default = NULL,
                      help = "Optional: Process only this K value")
  parser$add_argument("--fam", default = "",
                      help = "PLINK .fam file")
  parser$add_argument("--popmap", default = "",
                      help = "Population mapping CSV")
  parser$add_argument("--coords", default = "",
                      help = "Population coordinates CSV")
  parser$add_argument("--output_dir", default = "results",
                      help = "Output directory")
  parser$add_argument("--parallel_workers", type = "integer", default = 4,
                      help = "Number of parallel workers")
  parser$add_argument("--dpi", type = "integer", default = 300,
                      help = "DPI resolution for PNG plots (maps)")

  # Spatial options
  parser$add_argument("--bbox", type = "character", default = NULL,
                      help = "Custom bbox as 'xmin,ymin,xmax,ymax' string")
  parser$add_argument("--padding", type = "double", default = NULL,
                      help = "Padding for auto bbox (default: 0.15)")

  # Labels
  parser$add_argument("--labels", action="store_true",
                      help = "Add population labels to maps")

  args <- parser$parse_args()

  conf_loaded <- FALSE
  if (!is.null(args$conf_file) && file.exists(args$conf_file)) {
    config <- yaml::read_yaml(args$conf_file)
    # Override args with config
    for (name in names(config)) {
      if (name %in% names(args)) {
        args[[name]] <- config[[name]]
      }
    }
    conf_loaded <- TRUE
    message("Configuration loaded from: ", args$conf_file)
  }

  # Parse bbox if character
  if (!is.null(args$bbox) && is.character(args$bbox)) {
    bbox_vals <- as.numeric(strsplit(args$bbox, ",")[[1]])
    if (length(bbox_vals) != 4 || any(is.na(bbox_vals))) {
      stop("Invalid bbox format. Use 'xmin,ymin,xmax,ymax'")
    }
    args$bbox <- list(xmin = bbox_vals[1], ymin = bbox_vals[2],
                      xmax = bbox_vals[3], ymax = bbox_vals[4])
  }

  # Validate mode
  if (is.null(args$q_file) && (is.null(args$input_dir) || args$input_dir == "")) {
    stop("Must provide either --q_file or --input_dir")
  }
  if (!is.null(args$q_file) && !is.null(args$input_dir) && args$input_dir != "") {
    warning("Both --q_file and --input_dir provided; using single mode (--q_file)")
  }

  # Configure parallel processing
  future::plan(future::multisession, workers = args$parallel_workers)

  # Execute pipeline
  tryCatch({
    execute_pipeline(args)
  }, error = function(e) {
    message("Pipeline failed: ", e$message)
    quit(status = 1)
  })
}