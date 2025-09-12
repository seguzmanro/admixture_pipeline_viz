#!/usr/bin/env Rscript
# Admixture Spatial Visualization Pipeline (Refactored)
# Author: [Sebastian Guzman]
# Description: Integrates ADMIXTURE results with spatial data using modular helpers

# Load Core Packages ----------------------------------------------------------
#.libPaths('~/R/tmap_v4_lib/')
suppressPackageStartupMessages({
  library(argparse)
  library(dotenv)
  library(sf)
  library(terra)
  library(tmap)
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(future.apply)
  library(purrr)
})

# Load Custom Modules ---------------------------------------------------------
source("spatial_helpers.R")
source("plot_formats.R")

# Configuration ---------------------------------------------------------------
initialize_environment <- function() {
  #' Load environment variables and validate spatial setup
  #' 
  #' @return NULL
  
  project_home <- '~/Documents/Sebastian/UFRGS/Doctorado_2023/Pipelines/admixture_pipeline'
  dotenv::load_dot_env(file.path(project_home, ".paxil_env"))
  
}

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

process_admixture_results <- function(input_dir, sample_ids) {
  #' Process ADMIXTURE output files with parallelization
  #' 
  #' @return List containing CV errors and Q matrices

  # Extract CV errors (unchanged)
  extract_cv_error <- function(log_file) {
    content <- readLines(log_file)
    cv_line <- grep("^CV error", content, value = TRUE)
    k <- as.integer(sub(".*K=(\\d+).*", "\\1", cv_line))
    error <- as.numeric(sub("CV error \\(.*\\): (.*)", "\\1", cv_line))
    data.frame(K = k, CV_error = error)
  }
  
  # Process log files using absolute paths
  log_files <- list.files(input_dir, pattern = "*\\.out$", full.names = TRUE)  # Changed pattern
  if (length(log_files) == 0) {
    stop("No log files found in: ", input_dir)
  }
  cv_err <- map_dfr(log_files, extract_cv_error)
  
  # Process Q files with absolute paths
  q_files <- list.files(input_dir, pattern = "\\.Q$", full.names = TRUE)
  
  # Verify files actually exist
  missing_files <- q_files[!file.exists(q_files)]
  if (length(missing_files) > 0) {
    stop("Q files not found: ", paste(basename(missing_files), collapse = ", "))
  }
  
  admix_k <- future_lapply(q_files, function(file) {
    k_val <- as.integer(sub(".*\\.(\\d+)\\.Q$", "\\1", basename(file)))
    if (k_val == 1) return(NULL)
    
    q_mat <- as.matrix(read.table(file))
    rownames(q_mat) <- sample_ids
    colnames(q_mat) <- paste0("Cluster", 1:k_val)
    q_mat
  })
  
  # Name results using file basenames
  names(admix_k) <- sub("\\.Q$", "", basename(q_files))
  compact(admix_k)

  return(list(CV_errors = cv_err, Q_matrices = admix_k))
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
  
  # Process ADMIXTURE results
  admix_results <- process_admixture_results(
    args$input_dir,
    data$fam
  )
  
  # Determine which K values to process based on the --k_value argument
  k_names_to_process <- names(admix_results$Q_matrices)
  
  if (!is.null(args$k_value)) {
    message("User specified K = ", args$k_value, ". Filtering results.")
    
    # Find the name in the list corresponding to the specified K
    # This pattern looks for ".K$" at the end of the string to be precise
    target_k_name <- grep(paste0("\\.", args$k_value, "$"), k_names_to_process, value = TRUE)
    
    if (length(target_k_name) == 0) {
      stop("Specified K value (", args$k_value, ") not found in the input directory.")
    }
    
    # Overwrite the list of K's to process with only the target one
    k_names_to_process <- target_k_name
  }

  if (is.null(args$k_value)) {
    # Generate CV plot using plot_formats
    cv_plot <- format_cv_plot(admix_results$CV_errors)
    save_genomics_plot(
      cv_plot, 
      file.path(args$output_dir, "cv_error_plot.pdf"),
      width = 8, height = 6
    )
  }
  # Generate admixture barplots (compoplots) for each K
  for (k_name in k_names_to_process) {
    k_value <- as.integer(sub(".*\\.(\\d+)", "\\1", k_name))
    
    if (k_value > 1) { # Skip K=1
      # Create the barplot
      bar_plot <- format_admixture_barplot(
        admix_results$Q_matrices[[k_name]],
        data$popmap,
        k_value
      )
      
      # Save as PDF
      save_genomics_plot(
        bar_plot,
        file.path(args$output_dir, paste0("admixture_barplot_K", k_value, ".pdf")),
        width = 10, 
        height = 6
      )
      
      message("Generated barplot for K = ", k_value)
    }
  }
  # Load spatial data using spatial helpers
  spatial_data <- load_basemap_data(data$popcoords, lon_col = 'lon', lat_col = 'lat', elev_zoom = 5)

  # Create adaptive bounding box
  pop_sf <- create_population_sf(data$popcoords)
  adapt_bbox <- adaptive_bbox(pop_sf, padding=0.15)
  study_bbox <- st_bbox(c(adapt_bbox$xmin, adapt_bbox$xmax, adapt_bbox$ymax, adapt_bbox$ymin), crs = st_crs(4326))

  # Generate admixture maps
  base_map <- create_base_map(spatial_data, study_bbox)

  # Process each K value sequentially instead of in parallel
  for (k_name in k_names_to_process) {
    k_val <- as.integer(sub(".*\\.(\\d+)", "\\1", k_name))
  
    # Generate spatial components
    spatial_components <- generate_spatial_map(
      admix_results$Q_matrices[[k_name]],
      data$popmap,
      data$popcoords,
      k_value = k_val
    )
    
    # Create the full map
    full_map <- base_map +
      tm_shape(spatial_components$pop_sf) +
      tm_symbols(shape='pop',
        shapes = spatial_components$pie_grobs,
        size = 0.7,
        border.lwd = 0.1,
        col_alpha = 0.5
      ) +
      tm_title(
        paste("Genetic Structure - K =", k_val),
        position = c("center", "top")
      )
    
    # Save high-resolution output
    tmap_save(
      full_map, 
      file.path(args$output_dir, paste0("admixture_map_K", k_val, ".png")),
      width = 10, height = 8, dpi = args$dpi
    )
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
  # Initialize environment
  initialize_environment()

  parser <- argparse::ArgumentParser(
    description = "Spatial Visualization of ADMIXTURE Results"
  )
  # Add arguments for input files and directories
  parser$add_argument("--input_dir", default = Sys.getenv("EXAMPLE_ADMIX_DIR"), 
                      help = "ADMIXTURE output directory")
  parser$add_argument("--k_value", type = "integer", default = NULL,
                      help = "Optional: Process only this K value (default: all K)")
  parser$add_argument("--fam", default = Sys.getenv("EXAMPLE_FAM"), 
                      help = "PLINK .fam file")
  parser$add_argument("--popmap", default = Sys.getenv("EXAMPLE_POPMAP"), 
                      help = "Population mapping CSV")
  parser$add_argument("--coords", default = Sys.getenv("EXAMPLE_COORDS"), 
                      help = "Population coordinates CSV")
  parser$add_argument("--output_dir", default = Sys.getenv("EXAMPLE_OUTPUT_DIR"), 
                      help = "Output directory")
  parser$add_argument("--parallel_workers", type = "integer", default = 4,
                      help = "Number of parallel workers")
  parser$add_argument("--dpi", type = "integer", default = 300,
                      help = "dpi resolution for PNG plots (mainly the generated maps)")
  
  args <- parser$parse_args()
  
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