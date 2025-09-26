#' Data Processing Functions for ADMIXTURE Results
#'
#' This module contains functions for reading, validating, and processing
#' ADMIXTURE input files and related data.

#' Read and validate input files with improved error handling
#'
#' This function reads and validates the three main input files required for
#' ADMIXTURE visualization: PLINK .fam file, population mapping file, and
#' population coordinates file. It performs comprehensive validation to ensure
#' data consistency.
#'
#' @param fam_path Path to PLINK .fam file (character).
#'   Standard PLINK format with sample information.
#' @param popmap_path Path to population mapping CSV file (character).
#'   Two-column CSV file with sample IDs and population assignments.
#' @param coords_path Path to population coordinates CSV file (character).
#'   Three-column CSV file with population names, longitude, and latitude.
#' @return List of validated data frames containing:
#'   \itemize{
#'     \item{fam}{Vector of sample IDs from .fam file}
#'     \item{popmap}{Data frame with sample-to-population mapping}
#'     \item{popcoords}{Data frame with population coordinates}
#'   }
#' @details
#' The function performs the following validations:
#' \itemize{
#'   \item{Checks that all input files exist and are readable}
#'   \item{Validates that sample counts match between .fam and popmap files}
#'   \item{Ensures all samples in .fam file have population assignments}
#'   \item{Verifies coordinate data is complete (no missing lon/lat values)}
#'   \item{Standardizes column names for consistent processing}
#' }
#' @examples
#' \dontrun{
#' # Read and validate input files
#' data <- read_input_data(
#'   fam_path = "population.fam",
#'   popmap_path = "population_popmap.csv",
#'   coords_path = "population_coordinates.csv"
#' )
#'
#' # Access validated data
#' sample_ids <- data$fam
#' popmap <- data$popmap
#' coordinates <- data$popcoords
#' }
#' @export
read_input_data <- function(fam_path, popmap_path, coords_path) {
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

#' Extract CV error from ADMIXTURE log file
#'
#' @param log_file Path to .out file
#' @return Data frame with K and CV_error
#' @export
extract_cv_error <- function(log_file) {
  content <- readLines(log_file)
  cv_line <- grep("^CV error", content, value = TRUE)
  if (length(cv_line) == 0) {
    return(data.frame(K = NA, CV_error = NA))
  }
  k <- as.integer(sub(".*K=(\\d+).*", "\\1", cv_line))
  error <- as.numeric(sub("CV error \\(.*\\): (.*)", "\\1", cv_line))
  data.frame(K = k, CV_error = error)
}

#' Read Q file, detect K from columns, validate row count
#'
#' @param q_path Path to .Q file
#' @param sample_ids Vector of sample IDs from FAM
#' @return Q matrix or NULL if K==1
#' @export
read_q_file <- function(q_path, sample_ids) {
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

#' Process Q files (single or batch)
#'
#' @param input_dir_or_qfile Directory or single file path
#' @param sample_ids Sample IDs
#' @param is_single Logical: single file mode?
#' @return List of Q matrices
#' @export
process_q_files <- function(input_dir_or_qfile, sample_ids, is_single = FALSE) {
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

    q_matrices <- future.apply::future_lapply(q_files, read_q_file, sample_ids = sample_ids)
    names(q_matrices) <- sub("\\.Q$", "", basename(q_files))
    q_matrices <- purrr::compact(q_matrices)

    if (length(q_matrices) == 0) {
      stop("No valid Q files with K > 1 found.")
    }
    return(q_matrices)
  }
}