#' Spatial Data Processing Functions
#'
#' This module contains functions for processing spatial data and creating
#' geographic visualizations of ADMIXTURE results.

#' Get OSM water data with caching
#'
#' @param bbox_poly The bounding box polygon
#' @param osm_cache_dir Cache directory for OSM data
#' @return List containing rivers, lakes_poly, and lakes_multi
#' @export
get_osm_water_data <- function(bbox_poly, osm_cache_dir) {
  if (!dir.exists(osm_cache_dir)) {
    dir.create(osm_cache_dir, recursive = TRUE)
  }

  # Create a unique key for this bounding box
  bbox_key <- digest::digest(bbox_poly)
  cache_file <- file.path(osm_cache_dir, paste0(bbox_key, ".rds"))

  # Check if cached data exists
  if (file.exists(cache_file)) {
    message("Loading OSM water data from cache")
    return(readRDS(cache_file))
  }

  message("Downloading OSM water data (this may take a moment)")

  # Query for RIVERS
  rivers_query <- osmdata::opq(bbox = bbox_poly) |>
    osmdata::add_osm_feature(key = 'water', value = 'river') |>
    osmdata::osmdata_sf()

  # Query for LAKES
  lakes_query <- osmdata::opq(bbox = bbox_poly) |>
    osmdata::add_osm_feature(key = 'natural', value = 'water') |>
    osmdata::add_osm_feature(key = 'water', value = 'lake') |>
    osmdata::osmdata_sf()

  # Extract the features
  water_data <- list(
    rivers = rivers_query$osm_lines,
    lakes_poly = lakes_query$osm_polygons,
    lakes_multi = lakes_query$osm_multipolygons
  )

  # Save to cache
  saveRDS(water_data, cache_file)
  message("OSM water data saved to cache")

  return(water_data)
}

#' Load basemap data (boundaries, lakes, rivers, elevation) around given points
#'
#' This function downloads and processes spatial data required for creating
#' geographic visualizations of population genetic structure. It combines data
#' from multiple sources including elevation, political boundaries, water features,
#' and geographic reference data.
#'
#' @param pop_coords An sf object (POINT geometry) or data.frame with lon/lat columns.
#'   Population coordinates in WGS84 (EPSG:4326) coordinate system.
#' @param lon_col Character, column name for longitude if pop_coords is a data.frame
#'   (default: "lon").
#' @param lat_col Character, column name for latitude if pop_coords is a data.frame
#'   (default: "lat").
#' @param padding Numeric, degrees to expand the bbox around the population coordinates
#'   (default: 2). Controls how much area around the population points is included.
#' @param simplify_tolerance Numeric, tolerance for geometry simplification in degrees
#'   (default: 0.01). Higher values create simpler geometries but may reduce detail.
#' @param elev_zoom Integer, zoom level for elevation tiles (1–14; higher = finer, heavier).
#'   Controls the resolution of elevation data (default: 5).
#' @param cache_dir Directory for caching spatial data (character, required).
#'   Absolute path where downloaded spatial data will be cached for faster subsequent runs.
#' @return A list with sf/SpatRaster objects containing:
#'   \itemize{
#'     \item{countries}{Country/state boundaries as sf polygons}
#'     \item{lakes}{Natural Earth lake polygons}
#'     \item{rivers}{Natural Earth river polylines}
#'     \item{osm_water}{List with OSM water features (rivers, lakes)}
#'     \item{altitude}{Elevation raster (SpatRaster)}
#'     \item{bbox}{Bounding box polygon for the study area}
#'   }
#' @details
#' The function automatically downloads and caches spatial data from multiple sources:
#' \itemize{
#'   \item{Elevation data from SRTM (Shuttle Radar Topography Mission)}
#'   \item{Country boundaries from Natural Earth}
#'   \item{Water features from Natural Earth and OpenStreetMap}
#' }
#' Data is cached locally to speed up subsequent analyses of the same geographic region.
#' @examples
#' \dontrun{
#' # Load spatial data for population coordinates
#' coords <- data.frame(
#'   pop = c("POP1", "POP2", "POP3"),
#'   lon = c(-56.1, -55.8, -55.5),
#'   lat = c(-32.4, -32.1, -31.8)
#' )
#'
#' spatial_data <- load_basemap_data(
#'   pop_coords = coords,
#'   cache_dir = "cache_maps",
#'   elev_zoom = 8,  # Higher resolution elevation
#'   padding = 1.5   # More padding around points
#' )
#' }
#' @export
load_basemap_data <- function(pop_coords,
                              lon_col = "lon",
                              lat_col = "lat",
                              padding = 2,
                              simplify_tolerance = 0.01,
                              elev_zoom = 5,
                              cache_dir = NULL) {

  if (is.null(cache_dir)) {
    stop("An absolute path must be provided for 'cache_dir'.")
  }

  # Ensure sf object
  if (!inherits(pop_coords, "sf")) {
    if (!all(c(lon_col, lat_col) %in% names(pop_coords))) {
      stop("If pop_coords is a data.frame, it must have columns: ", lon_col, ", ", lat_col)
    }
    pop_coords <- sf::st_as_sf(pop_coords, coords = c(lon_col, lat_col), crs = 4326)
  }

  # Compute bbox + padding
  bbox <- sf::st_bbox(pop_coords)
  bbox_poly <- bbox + c(-padding,-padding,padding,padding)

  # Make sure cache dir exists
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  # Countries
  countries <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf") |>
    sf::st_make_valid() |>
    sf::st_crop(bbox_poly) |>
    sf::st_simplify(dTolerance = simplify_tolerance, preserveTopology = TRUE)

  # Natural Earth Lakes
  ne_lakes <- rnaturalearth::ne_download(scale = 10, type = "lakes", category = "physical",
                           returnclass = "sf") |>
    sf::st_make_valid() |>
    sf::st_crop(bbox_poly)

  # Natural Earth Rivers
  ne_rivers <- rnaturalearth::ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical",
                            returnclass = "sf") |>
    sf::st_make_valid() |>
    sf::st_crop(bbox_poly)

  # OSM water data
  osm_water_data <- get_osm_water_data(bbox_poly, osm_cache_dir = file.path(cache_dir, "osm_water"))

  # OSM Rivers (typically represented as lines)
  osm_rivers_lines <- osm_water_data$rivers |>
    dplyr::select(geometry) |>
    sf::st_as_sf() |>
    sf::st_make_valid()

  # OSM Lakes (typically represented as polygons)
  osm_lakes_polygons <- osm_water_data$lakes_poly |>
    dplyr::select(geometry) |>
    sf::st_as_sf() |>
    sf::st_make_valid()
  if (!is.null(osm_lakes_polygons)) {
    osm_lakes_polygons$area <- sf::st_area(osm_lakes_polygons)
    osm_large_lakes <- osm_lakes_polygons |> dplyr::filter(area > units::set_units(1, km^2))
  }

  # Also check for multipolygons which might contain additional water features
  osm_lakes_multipolygons <- osm_water_data$lakes_multi |>
    dplyr::select(geometry) |>
    sf::st_as_sf() |>
    sf::st_make_valid()

  water_features <- list(
    osm_rivers = osm_rivers_lines,
    osm_lakes_poly = osm_large_lakes,
    osm_lakes_multi = osm_lakes_multipolygons
  )


  # Elevation (SpatRaster; cache by bbox_poly and zoom level)
  cached_elev <- file.path(cache_dir, paste0(
    "elev", paste0(round(bbox_poly, 4), collapse = ""), '_', elev_zoom, ".tif"
  ))
  if (file.exists(cached_elev)) {
    elev <- terra::rast(cached_elev)
  } else {
    elev <- elevatr::get_elev_raster(
      locations = pop_coords,
      z = elev_zoom,
      expand = 5,
      clip = "bbox"
    ) |> terra::crop(bbox_poly) |> terra::rast()

    elev[elev[] < 0] <- NA
    names(elev) <- "elevation"

    terra::writeRaster(elev, cached_elev, overwrite = TRUE)
  }

  return(list(
    countries = countries,
    lakes = ne_lakes,
    rivers = ne_rivers,
    osm_water = water_features,
    altitude = elev,
    bbox = bbox_poly
  ))
}

#' Create population coordinate sf object
#'
#' @param coords_df Data frame with population coordinates
#' @param crs Coordinate reference system (default: EPSG:4326)
#' @return sf object with point geometries
#' @export
create_population_sf <- function(coords_df, crs = 4326) {
  required_cols <- c("pop", "lon", "lat")
  if (!all(required_cols %in% names(coords_df))) {
    stop("coords_df must contain columns: ", paste(required_cols, collapse = ", "))
  }

  coords_df |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = crs) |>
    dplyr::distinct(pop, .keep_all = TRUE)  # Ensure unique populations
}

#' Calculate adaptive bounding box from population coordinates
#'
#' @param pop_sf Population sf object
#' @param padding Percentage padding around points (default: 15%)
#' @return Bounding box (xmin, ymin, xmax, ymax)
#' @export
adaptive_bbox <- function(pop_sf, padding = 0.12) {
  bbox <- sf::st_bbox(pop_sf)
  x_range <- bbox["xmax"] - bbox["xmin"]
  y_range <- bbox["ymax"] - bbox["ymin"]

  list(
    xmin = bbox["xmin"] - padding * x_range,
    xmax = bbox["xmax"] + padding * x_range,
    ymin = bbox["ymin"] - padding * y_range,
    ymax = bbox["ymax"] + padding * y_range
  )
}

#' Create reusable base map template
#'
#' @param spatial_data Spatial data list from load_basemap_data
#' @param bbox Bounding box
#' @return tmap object
#' @export
create_base_map <- function(spatial_data, bbox) {
  tm <- tmap::tm_shape(spatial_data$altitude, bbox = bbox, raster.downsample = FALSE) +
      tmap::tm_raster(col='elevation', col.scale = tmap::tm_scale_continuous(values='brewer.greys', value.na='lightblue'), col_alpha = 0.9, col.legend=tmap::tm_legend(show=FALSE)) +
    tmap::tm_shape(spatial_data$countries) +
      tmap::tm_borders(lwd = 1.2)

  if (!is.null(spatial_data$osm_water$osm_rivers)) {
    tm <- tm + tmap::tm_shape(spatial_data$osm_water$osm_rivers) +
      tmap::tm_lines(col = "darkblue", lwd = 1.2, col_alpha = 0.7)
  }

  if (!is.null(spatial_data$osm_water$osm_lakes_poly)) {
    tm <- tm + tmap::tm_shape(spatial_data$osm_water$osm_lakes_poly) +
      tmap::tm_polygons(col = "darkblue", fill = "lightblue", lwd = 1.2, fill_alpha = 0.7)
  }

  if (!is.null(spatial_data$osm_water$osm_lakes_multi)) {
    tm <- tm + tmap::tm_shape(spatial_data$osm_water$osm_lakes_multi) +
      tmap::tm_polygons(col = "darkblue", fill = "lightblue", lwd = 1.2, fill_alpha = 0.7)
  }

  tm <- tm + tmap::tm_shape(spatial_data$rivers) +
      tmap::tm_lines(col="darkblue", lwd=1.2) +
    tmap::tm_shape(spatial_data$lakes) +
      tmap::tm_polygons(col="darkblue", fill = "lightblue", fill_alpha = 0.9) +
    tmap::tm_layout(
      legend.show = FALSE
    ) +
    tmap::tm_compass(position = c("right", "top")) +
    tmap::tm_scalebar(position = c("left", "bottom"))

  return(tm)
}

#' Generate spatial admixture-type visualization components for mapping
#'
#' @param admix_df Admixture matrix data frame
#' @param popmap Population map data frame
#' @param popcoords Population coordinates data frame
#' @param k_value Number of clusters (K)
#' @return List with population SF object and pie chart grobs
#' @export
generate_spatial_map <- function(admix_df, popmap, popcoords, k_value) {
  # Convert admixture matrix to data frame with sample IDs
  admix_df <- as.data.frame(admix_df)
  admix_df$sample <- rownames(admix_df)

  # Merge with popmap to get population assignments
  admix_pop <- admix_df |>
    dplyr::left_join(popmap, by = c("sample" = "indv"))

  # Calculate population-level admixture proportions
  pop_admix <- admix_pop |>
    dplyr::group_by(pop) |>
    dplyr::summarise(dplyr::across(dplyr::starts_with("Cluster"), mean)) |>
    dplyr::left_join(popcoords, by = "pop")

  pop_admix <- as.data.frame(pop_admix)

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
      ggplot2::ggplotGrob(population_pie_chart(proportions, colors = genetic_palette(k_value)))
    }),
    pop_admix$pop
  )

  return(list(pop_sf = pop_sf, pie_grobs = pie_grobs))
}