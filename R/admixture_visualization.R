#!/usr/bin/env Rscript
# Modular Admixture Spatial Visualization Pipeline
# Author: [Sebastian Guzman] (modified for modularity)
# Description: Generates visualization maps from individual or batch Q files, with optional CV error plots.
# K value is determined from Q file content (number of columns).
# Uses config.yaml for configuration; CLI overrides config.

get_script_dir <- function() {
  #' Get the directory of the currently running R script
  #'
  #' @return The directory of the R script.
  #' @note This function is designed to work when the script is run with Rscript.
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    stop("Cannot determine script path. Please run with Rscript.")
  }
  return(dirname(sub("--file=", "", file_arg[1])))
}

# Determine the script's directory
script_dir <- get_script_dir()

# Load Core Packages ----------------------------------------------------------

suppressPackageStartupMessages({
  library(argparse)
  library(yaml)
  library(sf)
  library(terra)
  library(tmap)
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(future.apply)
  library(purrr)
  library(tools)
})

# Load Custom Modules ---------------------------------------------------------
source(file.path(script_dir, "spatial_helpers.R"))
source(file.path(script_dir, "plot_formats.R"))

# Data Processing Functions ---------------------------------------------------
read_input_data <- function(fam_path, popmap_path, coords_path) {
  #' Read and validate input files with improved error handling
  #' 
  #' @return List of validated data frames
  
  tryCatch({
    # Read files
    fam <- read.table(fam_path, stringsAsFactors = FALSE)[, 2]
    popmap <- read.csv(popmap_path, header = TRUE)[, 1:2]
    popcoords <- read.csv(coords_path, header = TRUE)[, 1:3]
    
    # Standardize column names
    colnames(popmap) <- c("indv", "pop")
    colnames(popcoords) <- c("pop", "lon", "lat")
    
    # Validate structure
    if (nrow(popmap) != length(fam)) {
      stop("FAM file contains ", length(fam), " samples but popmap has ", 
           nrow(popmap), " entries")
    }
    
    if (!all(fam %in% popmap$indv)) {
      missing <- setdiff(fam, popmap$indv)
      stop(length(missing), " samples in FAM not found in popmap. First missing: ", 
           missing[1])
    }

    # Add coordinate validation
    if (any(is.na(popcoords$lon))) {
      stop("Missing longitude values for populations: ", 
           paste(popcoords$pop[is.na(popcoords$lon)], collapse = ", "))
    }
    if (any(is.na(popcoords$lat))) {
      stop("Missing latitude values for populations: ", 
           paste(popcoords$pop[is.na(popcoords$lat)], collapse = ", "))
    }
    
    return(list(
      fam = fam,
      popmap = popmap,
      popcoords = popcoords
    ))
  }, error = function(e) {
    stop("Input data error: ", e$message)
  })
}

extract_cv_error <- function(log_file) {
  #' Extract CV error from ADMIXTURE log file
  #' 
  #' @param log_file Path to .out file
  #' @return Data frame with K and CV_error
  
  content <- readLines(log_file)
  cv_line <- grep("^CV error", content, value = TRUE)
  if (length(cv_line) == 0) {
    return(data.frame(K = NA, CV_error = NA))
  }
  k <- as.integer(sub(".*K=(\\d+).*", "\\1", cv_line))
  error <- as.numeric(sub("CV error \\(.*\\): (.*)", "\\1", cv_line))
  data.frame(K = k, CV_error = error)
}

read_q_file <- function(q_path, sample_ids) {
  #' Read Q file, detect K from columns, validate row count
  #' 
  #' @param q_path Path to .Q file
  #' @param sample_ids Vector of sample IDs from FAM
  #' @return Q matrix or NULL if K==1
  
  if (!file.exists(q_path)) {
    stop("Q file not found: ", q_path)
  }
  
  q_mat <- as.matrix(read.table(q_path, header = FALSE))
  k_val <- ncol(q_mat)
  
  if (k_val == 1) {
    return(NULL)
  }
  
  if (nrow(q_mat) != length(sample_ids)) {
    stop("Q file has ", nrow(q_mat), " rows but FAM has ", length(sample_ids), " samples")
  }
  
  rownames(q_mat) <- sample_ids
  colnames(q_mat) <- paste0("Cluster", 1:k_val)
  return(q_mat)
}

process_q_files <- function(input_dir_or_qfile, sample_ids, is_single = FALSE) {
  #' Process Q files (single or batch)
  #' 
  #' @param input_dir_or_qfile Directory or single file path
  #' @param sample_ids Sample IDs
  #' @param is_single Logical: single file mode?
  #' @return List of Q matrices
  
  if (is_single) {
    q_mat <- read_q_file(input_dir_or_qfile, sample_ids)
    if (is.null(q_mat)) {
      stop("Single Q file has K=1, which is not supported.")
    }
    return(list(single = q_mat))
  } else {
    q_files <- list.files(input_dir_or_qfile, pattern = "\\.Q$", full.names = TRUE)
    if (length(q_files) == 0) {
      stop("No .Q files found in directory: ", input_dir_or_qfile)
    }
    
    q_matrices <- future_lapply(q_files, read_q_file, sample_ids = sample_ids)
    names(q_matrices) <- sub("\\.Q$", "", basename(q_files))
    q_matrices <- compact(q_matrices)
    
    if (length(q_matrices) == 0) {
      stop("No valid Q files with K > 1 found.")
    }
    return(q_matrices)
  }
}

# Main Workflow ---------------------------------------------------------------
execute_pipeline <- function(args) {
  #' Main pipeline execution flow
  #' 
  #' @return NULL
  
  # Create output directory
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load and validate data
  data <- read_input_data(args$fam, args$popmap, args$coords)
  
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
    if (is.null(args$input_dir) || !dir.exists(args$input_dir)) {
      stop("Input directory required for batch mode and must exist.")
    }
    q_path <- args$input_dir
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
  if (!is.null(args$log_dir) && dir.exists(args$log_dir)) {
    log_files <- list.files(args$log_dir, pattern = "\\.out$", full.names = TRUE)
    if (length(log_files) > 0) {
      cv_err <- map_dfr(log_files, extract_cv_error)
      cv_err <- cv_err[!is.na(cv_err$K), ]  # Remove invalid
      if (nrow(cv_err) >= 2) {
        cv_plot <- format_cv_plot(cv_err)
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
  
  # Spatial setup
  cache_path <- file.path(script_dir, "cache_maps")
  message("Using cache directory: ", cache_path)
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
    study_bbox <- st_bbox(unlist(args$bbox), crs = st_crs(4326))
    message("Using custom bounding box: xmin=", args$bbox$xmin, ", ymin=", args$bbox$ymin, 
            ", xmax=", args$bbox$xmax, ", ymax=", args$bbox$ymax)
  } else {
    padding <- ifelse(is.null(args$padding), 0.15, args$padding)
    adapt_bbox <- adaptive_bbox(pop_sf, padding = padding)
    study_bbox <- st_bbox(c(adapt_bbox$xmin, adapt_bbox$xmax, adapt_bbox$ymax, adapt_bbox$ymin), crs = st_crs(4326))
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
        tm_shape(spatial_components$pop_sf) +
        tm_symbols(
          shape = 'pop',
          shapes = spatial_components$pie_grobs,
          size = 0.7,
          border.lwd = 0.1,
          col_alpha = 0.5
        ) +
        tm_title(
          paste("Genetic Structure - K =", k_val),
          position = c("center", "top")
        )
      
      # Add population labels if requested
      if (!is.null(args$labels) && args$labels) {
        full_map <- full_map + tm_text("pop", size = 0.7, ymod = 1.7, col = "black", fontface = "bold")
      }
      
      tmap_save(
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

# Command-line Interface ------------------------------------------------------
main <- function() {
  parser <- argparse::ArgumentParser(
    description = "Modular Spatial Visualization of ADMIXTURE Results"
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

# Execute ---------------------------------------------------------------------
if (!interactive()) {
  main()
}
