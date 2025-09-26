#' Main Visualization Functions
#'
#' This module contains the main workflow functions for creating
#' ADMIXTURE visualizations.

#' Main pipeline execution flow
#'
#' This is the core function that orchestrates the entire ADMIXTURE visualization
#' pipeline. It processes input files, generates visualizations, and saves outputs
#' to the specified directory.
#'
#' @param args List of arguments containing:
#'   \itemize{
#'     \item{q_file}{Path to single Q file (character)}
#'     \item{input_dir}{Directory with multiple Q files (character)}
#'     \item{fam}{Path to PLINK .fam file (character)}
#'     \item{popmap}{Path to population mapping CSV (character)}
#'     \item{coords}{Path to coordinates CSV (character)}
#'     \item{output_dir}{Output directory (character)}
#'     \item{k_value}{Specific K value to process (integer, optional)}
#'     \item{log_dir}{Directory with .out files for CV plot (character, optional)}
#'     \item{dpi}{Resolution for PNG plots (integer)}
#'     \item{labels}{Add population labels to maps (logical)}
#'     \item{bbox}{Custom bounding box (list, optional)}
#'     \item{padding}{Padding for adaptive bbox (numeric)}
#'     \item{parallel_workers}{Number of parallel workers (integer)}
#'     \item{cache_dir}{Directory for caching spatial data (character, optional)}
#'   }
#' @return NULL (invisibly)
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item{Validates and loads input data files}
#'   \item{Processes Q files (single or batch mode)}
#'   \item{Generates cross-validation error plots (if log files provided)}
#'   \item{Creates spatial basemap with elevation, water features, and boundaries}
#'   \item{Generates admixture bar plots for each K value}
#'   \item{Creates spatial admixture maps with pie charts}
#'   \item{Saves session information for reproducibility}
#' }
#' @examples
#' \dontrun{
#' args <- list(
#'   q_file = "admixture_results/pop.3.Q",
#'   fam = "data.fam",
#'   popmap = "popmap.csv",
#'   coords = "coordinates.csv",
#'   output_dir = "results",
#'   dpi = 300,
#'   labels = TRUE
#' )
#' execute_pipeline(args)
#' }
#' @export
execute_pipeline <- function(args) {
  # Create output directory - resolve extdata path if needed
  output_dir <- args$output_dir
  if (output_dir == "extdata/Pdepint/plot_results") {
    # Create output in the installed extdata directory
    extdata_dir <- system.file("extdata", package = "admixspatial")
    output_dir <- file.path(extdata_dir, "Pdepint", "plot_results")
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Load and validate data - resolve null paths to installed extdata
  fam_path <- args$fam
  popmap_path <- args$popmap
  coords_path <- args$coords

  if (is.null(fam_path)) {
    fam_path <- system.file("extdata/Pdepint/Pdepint_M095_noSele.fam", package = "admixspatial")
  }
  if (is.null(popmap_path)) {
    popmap_path <- system.file("extdata/Pdepint/Pops_Pdepint.csv", package = "admixspatial")
  }
  if (is.null(coords_path)) {
    coords_path <- system.file("extdata/Pdepint/Pdepint_Lat_Long_Coords.csv", package = "admixspatial")
  }

  data <- read_input_data(fam_path, popmap_path, coords_path)

  # Determine mode and process Q files
  if (!is.null(args$q_file)) {
    message("Processing single Q file mode")
    is_single <- TRUE
    q_path <- args$q_file
    if (!file.exists(q_path)) {
      stop("Single Q file not found: ", q_path)
    }
  } else {
    message("Processing batch mode from directory")
    is_single <- FALSE
    if (is.null(args$input_dir)) {
      # Use installed extdata directory
      q_path <- system.file("extdata/Pdepint", package = "admixspatial")
    } else {
      q_path <- args$input_dir
    }
    if (!dir.exists(q_path)) {
      stop("Input directory not found: ", q_path)
    }
  }

  Q_matrices <- process_q_files(q_path, data$fam, is_single)

  # Filter by k_value if specified
  k_names_to_process <- names(Q_matrices)
  if (!is.null(args$k_value)) {
    actual_ks <- sapply(Q_matrices, ncol)
    matching_names <- names(Q_matrices)[actual_ks == args$k_value]
    if (length(matching_names) == 0) {
      stop("No Q file(s) with K = ", args$k_value, " found.")
    }
    k_names_to_process <- matching_names
    message("Filtered to K = ", args$k_value)
  }

  # Process CV errors if log_dir provided
  cv_plotted <- FALSE
  log_dir <- args$log_dir
  if (is.null(log_dir)) {
    # Use installed extdata directory for log files
    log_dir <- system.file("extdata/Pdepint", package = "admixspatial")
  }
  if (dir.exists(log_dir)) {
    log_files <- list.files(log_dir, pattern = "\\.out$", full.names = TRUE)
    if (length(log_files) > 0) {
      cv_err <- purrr::map_dfr(log_files, extract_cv_error)
      cv_err <- cv_err[!is.na(cv_err$K), ]  # Remove invalid
      if (nrow(cv_err) >= 2) {
        cv_plot <- format_model_selection_plot(cv_err, title="ADMIXTURE Cross-Validation Error", method='min')
        save_genomics_plot(
          cv_plot,
          file.path(args$output_dir, "cv_error_plot.pdf"),
          width = 8, height = 6
        )
        message("Generated CV error plot")
        cv_plotted <- TRUE
      } else {
        message("Insufficient valid .out files for CV plot (need at least 2 K values)")
      }
    } else {
      message("No .out files found in log_dir")
    }
  }

  # Spatial setup - use provided cache_dir or default to package library
  if (is.null(args$cache_dir)) {
    package_dir <- find.package("admixspatial")
    cache_path <- file.path(package_dir, "cache_maps")
    message("Using default cache directory (package library): ", cache_path)
  } else {
    cache_path <- args$cache_dir
    message("Using custom cache directory: ", cache_path)
    # Ensure custom cache directory exists
    if (!dir.exists(cache_path)) {
      dir.create(cache_path, recursive = TRUE)
    }
  }
  spatial_data <- load_basemap_data(
    data$popcoords,
    lon_col = 'lon',
    lat_col = 'lat',
    elev_zoom = 5,
    cache_dir = cache_path
  )

  # Create bounding box (custom or adaptive)
  pop_sf <- create_population_sf(data$popcoords)
  if (!is.null(args$bbox) && is.list(args$bbox) && length(args$bbox) == 4) {
    study_bbox <- sf::st_bbox(unlist(args$bbox), crs = sf::st_crs(4326))
    message("Using custom bounding box: xmin=", args$bbox$xmin, ", ymin=", args$bbox$ymin,
            ", xmax=", args$bbox$xmax, ", ymax=", args$bbox$ymax)
  } else {
    padding <- ifelse(is.null(args$padding), 0.15, args$padding)
    adapt_bbox <- adaptive_bbox(pop_sf, padding = padding)
    study_bbox <- sf::st_bbox(c(adapt_bbox$xmin, adapt_bbox$xmax, adapt_bbox$ymax, adapt_bbox$ymin), crs = sf::st_crs(4326))
    message("Using adaptive bounding box with padding: ", padding)
  }

  # Generate base map
  base_map <- create_base_map(spatial_data, study_bbox)

  # Generate visualizations for each selected K
  for (k_name in k_names_to_process) {
    k_val <- ncol(Q_matrices[[k_name]])

    if (k_val > 1) {
      # Admixture barplot
      bar_plot <- format_admixture_barplot(
        Q_matrices[[k_name]],
        data$popmap,
        k_val
      )

      save_genomics_plot(
        bar_plot,
        file.path(args$output_dir, paste0("admixture_barplot_K", k_val, ".pdf")),
        width = 10,
        height = 6
      )

      message("Generated barplot for K = ", k_val)

      # Admixture map
      spatial_components <- generate_spatial_map(
        Q_matrices[[k_name]],
        data$popmap,
        data$popcoords,
        k_value = k_val
      )

      full_map <- base_map +
        tmap::tm_shape(spatial_components$pop_sf) +
        tmap::tm_symbols(
          shape = 'pop',
          shapes = spatial_components$pie_grobs,
          size = 0.7,
          border.lwd = 0.1,
          col_alpha = 0.5
        ) +
        tmap::tm_title(
          paste("Genetic Structure - K =", k_val),
          position = c("center", "top")
        )

      # Add population labels if requested
      if (!is.null(args$labels) && args$labels) {
        full_map <- full_map + tmap::tm_text("pop", size = 0.7, ymod = 1.7, col = "black", fontface = "bold")
      }

      tmap::tmap_save(
        full_map,
        file.path(args$output_dir, paste0("admixture_map_K", k_val, ".png")),
        width = 10, height = 8, dpi = args$dpi
      )

      message("Generated map for K = ", k_val)
    }
  }

  # Save session info for reproducibility
  writeLines(
    capture.output(sessionInfo()),
    file.path(args$output_dir, "session_info.txt")
  )

  message("Pipeline completed successfully. Outputs in: ", args$output_dir)
}

#' High-level wrapper function for library usage
#'
#' This is the main function for users who want to use the package programmatically.
#' It provides a comprehensive interface for generating ADMIXTURE visualizations
#' with full control over all parameters.
#'
#' @param q_file Path to single Q file (character, optional). If provided, processes
#'   only this file. Mutually exclusive with input_dir.
#' @param input_dir Directory containing multiple Q files (character, optional).
#'   If provided, processes all .Q files in the directory. Mutually exclusive with q_file.
#' @param fam_file Path to PLINK .fam file (character, required). Contains sample
#'   information and must match the Q file dimensions.
#' @param popmap_file Path to population mapping CSV (character, required).
#'   Two-column file with sample IDs and population assignments.
#' @param coords_file Path to coordinates CSV (character, required).
#'   Three-column file with population names, longitude, and latitude.
#' @param output_dir Output directory (character, default: "results").
#'   Directory where all output files will be saved.
#' @param k_value Specific K value to process (integer, optional).
#'   If provided, only processes Q files with this K value. If NULL, processes all K values.
#' @param log_dir Directory with .out files for CV plot (character, optional).
#'   If provided, generates cross-validation error plots from ADMIXTURE log files.
#' @param dpi Resolution for PNG plots (integer, default: 300).
#'   Higher values create sharper but larger image files.
#' @param labels Add population labels to maps (logical, default: FALSE).
#'   If TRUE, adds population names as text labels on the spatial maps.
#' @param bbox Custom bounding box as c(xmin, ymin, xmax, ymax) (numeric vector, optional).
#'   If provided, uses this bounding box instead of adaptive calculation.
#' @param padding Padding for adaptive bbox (numeric, default: 0.15).
#'   Fraction of the data range to add as padding around population coordinates.
#' @param parallel_workers Number of parallel workers (integer, default: 4).
#'   Number of CPU cores to use for parallel processing. Increase for faster processing
#'   on multi-core systems.
#' @param cache_dir Directory for caching spatial data (character, default: package library).
#'   Absolute path where downloaded spatial data will be cached. If NULL, uses the package
#'   library location for consistent caching across working directories.
#' @return NULL (invisibly). All outputs are saved to files in the output directory.
#' @details
#' This function provides a complete interface for generating ADMIXTURE visualizations.
#' It automatically:
#' \itemize{
#'   \item{Validates all input files and data consistency}
#'   \item{Downloads and caches spatial data (elevation, water features, boundaries)}
#'   \item{Generates admixture bar plots for each K value}
#'   \item{Creates spatial maps with pie charts showing population admixture}
#'   \item{Produces cross-validation error plots (if log files provided)}
#'   \item{Saves session information for reproducibility}
#' }
#' @examples
#' \dontrun{
#' # Single Q file with custom settings
#' visualize_admixture(
#'   q_file = "admixture_results/population.3.Q",
#'   fam_file = "data/population.fam",
#'   popmap_file = "data/population_popmap.csv",
#'   coords_file = "data/population_coordinates.csv",
#'   output_dir = "admixture_results",
#'   k_value = 3,
#'   labels = TRUE,
#'   dpi = 300,
#'   parallel_workers = 8
#' )
#'
#' # Batch processing with CV error plots
#' visualize_admixture(
#'   input_dir = "admixture_results",
#'   log_dir = "admixture_results",
#'   fam_file = "data/population.fam",
#'   popmap_file = "data/population_popmap.csv",
#'   coords_file = "data/population_coordinates.csv",
#'   output_dir = "batch_results",
#'   parallel_workers = 4
#' )
#' }
#' @export
visualize_admixture <- function(q_file = NULL,
                                input_dir = NULL,
                                fam_file,
                                popmap_file,
                                coords_file,
                                output_dir = "results",
                                k_value = NULL,
                                log_dir = NULL,
                                dpi = 300,
                                labels = FALSE,
                                bbox = NULL,
                                padding = 0.15,
                                parallel_workers = 4,
                                cache_dir = NULL) {

  # Validate inputs
  if (is.null(q_file) && (is.null(input_dir) || !dir.exists(input_dir))) {
    stop("Must provide either q_file or input_dir")
  }

  # Configure parallel processing
  future::plan(future::multisession, workers = parallel_workers)

  # Create argument list
  args <- list(
    q_file = q_file,
    input_dir = input_dir,
    fam = fam_file,
    popmap = popmap_file,
    coords = coords_file,
    output_dir = output_dir,
    k_value = k_value,
    log_dir = log_dir,
    dpi = dpi,
    labels = labels,
    bbox = bbox,
    padding = padding,
    parallel_workers = parallel_workers,
    cache_dir = cache_dir
  )

  # Execute pipeline
  tryCatch({
    execute_pipeline(args)
  }, error = function(e) {
    stop("Pipeline failed: ", e$message)
  })
}

#' Quick visualization function with minimal parameters
#'
#' This is a convenience function for users who want to quickly generate
#' ADMIXTURE visualizations with sensible defaults. It enables population
#' labels by default and uses optimized settings for most use cases.
#'
#' @param q_file Path to Q file (character, required). Path to a single ADMIXTURE
#'   .Q file containing admixture proportions.
#' @param fam_file Path to .fam file (character, required). PLINK .fam file
#'   containing sample information.
#' @param popmap_file Path to popmap CSV (character, required). Two-column CSV
#'   file mapping sample IDs to population names.
#' @param coords_file Path to coordinates CSV (character, required). Three-column
#'   CSV file with population coordinates (pop, lon, lat).
#' @param output_dir Output directory (character, default: "admixture_results").
#'   Directory where visualization files will be saved.
#' @param cache_dir Directory for caching spatial data (character, default: package library).
#'   Absolute path where downloaded spatial data will be cached. If NULL, uses the package
#'   library location for consistent caching across working directories.
#' @return NULL (invisibly). All outputs are saved to files in the output directory.
#' @details
#' This function is designed for quick visualization with minimal configuration.
#' It automatically:
#' \itemize{
#'   \item{Enables population labels on maps}
#'   \item{Uses optimal DPI (300) for publication-quality figures}
#'   \item{Processes all K values found in the Q file}
#'   \item{Generates both bar plots and spatial maps}
#'   \item{Caches spatial data for faster subsequent runs}
#' }
#' @examples
#' \dontrun{
#' # Quick visualization with default settings
#' quick_visualize(
#'   q_file = "admixture_results/population.3.Q",
#'   fam_file = "data/population.fam",
#'   popmap_file = "data/population_popmap.csv",
#'   coords_file = "data/population_coordinates.csv"
#' )
#'
#' # Custom output directory
#' quick_visualize(
#'   q_file = "results/pop.4.Q",
#'   fam_file = "data.fam",
#'   popmap_file = "popmap.csv",
#'   coords_file = "coordinates.csv",
#'   output_dir = "my_results"
#' )
#' }
#' @export
quick_visualize <- function(q_file, fam_file, popmap_file, coords_file, output_dir = "admixture_results", cache_dir = NULL) {
  visualize_admixture(
    q_file = q_file,
    fam_file = fam_file,
    popmap_file = popmap_file,
    coords_file = coords_file,
    output_dir = output_dir,
    labels = TRUE,
    cache_dir = cache_dir
  )
}