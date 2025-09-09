#' Spatial Data Processing Utilities for Population Genetics Visualization
#' 
#' This module provides helper functions for loading and processing spatial data
#' used in genomic visualization pipelines. Designed for integration with
#' admixture mapping workflows.

library(sf)
library(terra)
library(plyr)
library(dplyr)

#' Load and process country boundary data
#'
#' @param countries Vector of country codes to load (URY, ARG, BRA)
#' @param simplify_tolerance Simplification tolerance in degrees (0 = no simplification)
#' @return Combined and processed sf object of country boundaries
load_country_data <- function(countries = c("URY", "ARG", "BRA"), 
                             simplify_tolerance = 0.01) {
  # Validate input
  if (!all(countries %in% c("URY", "ARG", "BRA"))) {
    stop("Supported countries: URY, ARG, BRA")
  }
  
  # Environment variable mapping
  env_vars <- c(
    URY = Sys.getenv("URY_ADM"),
    ARG = Sys.getenv("ARG_ADM"),
    BRA = Sys.getenv("BRA_ADM")
  )
  
  # Load and process each country
  country_data <- lapply(countries, function(code) {
    path <- env_vars[[code]]
    if (!file.exists(path)) {
      warning(paste("File not found for", code, ":", path))
      return(NULL)
    }
    
    suppressMessages(
      read_sf(path) %>% 
        select(admin = matches("ADM0|NAME_0")) %>% 
        mutate(country_code = code) %>%
        st_make_valid()
    )
  })
  
  # Combine and simplify
  combined <- do.call(rbind, compact(country_data)) %>%
    st_simplify(preserveTopology = TRUE, dTolerance = simplify_tolerance)
  
  return(combined)
}

#' Load and filter water body data
#'
#' @param countries Vector of country codes to load (URY, BRA)
#' @param filter_expr Expression to filter water bodies (e.g., "PERENNIAL")
#' @param tolerance Simplification tolerance in degrees
#' @return Processed sf object of water bodies
load_water_data <- function(countries = c("URY", "BRA"), 
                           filter_expr = "Perennial|Permanent",
                           tolerance = 0.005) {
  # Validate input
  valid_countries <- c("URY", "BRA")
  invalid <- setdiff(countries, valid_countries)
  if (length(invalid) > 0) {
    stop("Supported water data countries: ", paste(valid_countries, collapse = ", "))
  }
  
  # Environment variable mapping
  env_vars <- c(
    URY = Sys.getenv("URY_WAT"),
    BRA = Sys.getenv("BRA_WAT")
  )
  
  # Load and process each country's water data
  water_data <- lapply(countries, function(code) {
    path <- env_vars[[code]]
    if (!file.exists(path)) {
      warning(paste("Water data not found for", code, ":", path))
      return(NULL)
    }
    
    suppressMessages(
      read_sf(path) %>%
        filter(grepl(filter_expr, HYC_DESCRI, ignore.case = TRUE)) %>%
        mutate(country_code = code) %>%
        st_make_valid() %>%
        st_simplify(preserveTopology = TRUE, dTolerance = tolerance)
    )
  })
  
  return(do.call(rbind, compact(water_data)))
}

#' Load and merge elevation data
#'
#' @param countries Vector of country codes (URY, BRA)
#' @return Merged SpatRaster object
load_elevation_data <- function(countries = c("URY", "BRA")) {
  # Validate input
  valid_countries <- c("URY", "BRA")
  invalid <- setdiff(countries, valid_countries)
  if (length(invalid) > 0) {
    stop("Supported elevation data countries: ", paste(valid_countries, collapse = ", "))
  }
  
  # Environment variable mapping
  env_vars <- c(
    URY = Sys.getenv("URY_ALT"),
    BRA = Sys.getenv("BRA_ALT")
  )
  
  # Load and merge rasters
  elev_rasters <- lapply(countries, function(code) {
    path <- env_vars[[code]]
    
    # Check for existence - now supports both .tif and .grd
    if (!file.exists(path)) {
      # Try alternative extensions
      tif_path <- sub("\\.grd$", ".tif", path, ignore.case = TRUE)
      if (file.exists(tif_path)) {
        path <- tif_path
      } else {
        warning(paste("Elevation data not found for", code, ":", path))
        return(NULL)
      }
    }
    
    # Load with terra
    tryCatch({
      r <- rast(path)
      return(r)
    }, error = function(e) {
      warning(paste("Failed to load elevation for", code, ":", e$message))
      return(NULL)
    })
  })

  # Remove NULL elements and merge
  elev_rasters <- compact(elev_rasters)
  if (length(elev_rasters) == 0) return(NULL)
  
  merged <- do.call(merge, elev_rasters)
  names(merged) <- "elevation"
  
  return(merged)
}

#' Crop spatial objects to bounding box
#'
#' @param spatial_obj Spatial object (sf or SpatRaster)
#' @param bbox Bounding box (xmin, xmax, ymin, ymax)
#' @return Cropped spatial object
crop_to_bbox <- function(spatial_obj, bbox) {
  if (inherits(spatial_obj, "sf")) {
    return(st_crop(spatial_obj, bbox))
  } else if (inherits(spatial_obj, "SpatRaster")) {
    return(crop(spatial_obj, ext(bbox)))
  } else {
    stop("Unsupported spatial object type")
  }
}

#' Validate spatial environment variables
#'
#' @param required_vars Vector of required environment variables
#' @return TRUE if all exist, throws error otherwise
validate_spatial_env <- function(required_vars = c(
  "URY_ADM", "ARG_ADM", "BRA_ADM",
  "URY_WAT", "BRA_WAT",
  "URY_ALT", "BRA_ALT"
)) {
  missing_vars <- setdiff(required_vars, names(Sys.getenv()))
  if (length(missing_vars) > 0) {
    stop("Missing required environment variables:\n",
         paste("-", missing_vars, collapse = "\n"))
  }
  return(TRUE)
}

#' Generate standardized color palette for spatial elements
#'
#' @param element Spatial element type ("country", "water", "elevation")
#' @return List with fill and color parameters
get_spatial_palette <- function(element = c("country", "water", "elevation")) {
  element <- match.arg(element)
  
  palettes <- list(
    country = list(fill = NA, color = "black", lwd = 0.8),
    water = list(fill = "#a6cee3", color = NA, alpha = 0.6),
    elevation = list(palette = "Greys", alpha = 0.7, legend = FALSE)
  )
  
  return(palettes[[element]])
}

#' Create population coordinate sf object
#'
#' @param coords_df Data frame with population coordinates
#' @param crs Coordinate reference system (default: EPSG:4326)
#' @return sf object with point geometries
create_population_sf <- function(coords_df, crs = 4326) {
  required_cols <- c("pop", "lon", "lat")
  if (!all(required_cols %in% names(coords_df))) {
    stop("coords_df must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  coords_df %>%
    st_as_sf(coords = c("lon", "lat"), crs = crs) %>%
    distinct(pop, .keep_all = TRUE)  # Ensure unique populations
}

#' Calculate adaptive bounding box from population coordinates
#'
#' @param pop_sf Population sf object
#' @param padding Percentage padding around points (default: 15%)
#' @return Bounding box (xmin, ymin, xmax, ymax)
adaptive_bbox <- function(pop_sf, padding = 0.15) {
  bbox <- st_bbox(pop_sf)
  x_range <- bbox["xmax"] - bbox["xmin"]
  y_range <- bbox["ymax"] - bbox["ymin"]
  
  list(
    xmin = bbox["xmin"] - padding * x_range,
    xmax = bbox["xmax"] + padding * x_range,
    ymin = bbox["ymin"] - padding * y_range,
    ymax = bbox["ymax"] + padding * y_range
  )
}

create_base_map <- function(spatial_data, bbox) {
  #' Create reusable base map template
  #' 
  #' @return tmap object
  
  tm_shape(spatial_data$altitude, bbox = bbox, raster.downsample = FALSE) + 
    tm_raster(style = "cont", palette = "Greys", legend.show = FALSE) +
  tm_shape(spatial_data$countries) + 
    tm_borders(lwd = 1.2, col = "black") +
  tm_shape(spatial_data$water) + 
    tm_polygons(col = "dodgerblue", alpha = 0.3, border.col = NA) +
  tm_layout(
    legend.position = c("right", "bottom")
  ) +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom"))
}

generate_spatial_map <- function(admix_df, popmap, popcoords, spatial_data, k_value, bbox) {
  #' Generate spatial admixture-type visualization components for mapping. Specifically, creates the pie-chart grobs representing population admixture-type proportions for each population and for a determined k_value.
  #' 
  #' @return List with population SF object and pie chart grobs
  
  # Convert admixture matrix to data frame with sample IDs
  admix_df <- as.data.frame(admix_df)
  admix_df$sample <- rownames(admix_df)
  
  # Merge with popmap to get population assignments
  admix_pop <- admix_df %>%
    left_join(popmap, by = c("sample" = "indv"))
  
  # Calculate population-level admixture proportions
  pop_admix <- admix_pop %>%
    group_by(pop) %>%
    summarise(across(starts_with("Cluster"), mean)) %>%
    left_join(popcoords, by = "pop")
  
  # Check for missing coordinates
  missing_coords <- pop_admix$pop[is.na(pop_admix$lon) | is.na(pop_admix$lat)]
  if (length(missing_coords) > 0) {
    stop("Populations missing coordinates: ", 
         paste(missing_coords, collapse = ", "))
  }
  
  # Create population sf object
  pop_sf <- create_population_sf(pop_admix)
  row.names(pop_sf) <- pop_sf$pop
  
  # Generate pie charts
  pie_grobs <- setNames(
    lapply(split(pop_admix, pop_admix$pop), function(pop_data) {
      proportions <- as.numeric(pop_data[, grep("Cluster", names(pop_data))])
      names(proportions) <- colnames(pop_data)[grep("Cluster", names(pop_data))]
      ggplotGrob(population_pie_chart(proportions, colors = genetic_palette(k_value)))
    }),
    pop_admix$pop
  )
  
  return(list(pop_sf = pop_sf, pie_grobs = pie_grobs))
}