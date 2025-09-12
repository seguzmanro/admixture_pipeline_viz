# spatial_helpers.R

library(sf)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(elevatr)
library(terra)

#' Load basemap data (boundaries, lakes, rivers, elevation) around given points
#'
#' @param pop_coords An sf object (POINT geometry) or data.frame with lon/lat columns.
#' @param lon_col Character, column name for longitude if pop_coords is a data.frame.
#' @param lat_col Character, column name for latitude if pop_coords is a data.frame.
#' @param padding Numeric, degrees to expand the bbox around the population coordinates.
#' @param simplify_tolerance Numeric, tolerance for geometry simplification (degrees).
#' @param elev_zoom Integer, zoom level for elevation tiles (1–14; higher = finer, heavier).
#'
#' @return A list with sf/SpatRaster objects: countries, lakes, rivers, elevation, bbox
#'
load_basemap_data <- function(pop_coords,
                              lon_col = "lon",
                              lat_col = "lat",
                              padding = 2,
                              simplify_tolerance = 0.01,
                              elev_zoom = 5,
                              cache_dir = "cache_maps") {
  # Ensure sf object
  if (!inherits(pop_coords, "sf")) {
    if (!all(c(lon_col, lat_col) %in% names(pop_coords))) {
      stop("If pop_coords is a data.frame, it must have columns: ", lon_col, ", ", lat_col)
    }
    pop_coords <- st_as_sf(pop_coords, coords = c(lon_col, lat_col), crs = 4326)
  }
  
  # Compute bbox + padding
  bbox <- st_bbox(pop_coords)
  bbox_poly <- bbox + c(-padding,-padding,padding,padding)
  
  # Make sure cache dir exists
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  # Countries
  countries <- ne_countries(scale = 10, returnclass = "sf") %>%
    st_make_valid() %>%
    st_crop(bbox_poly) %>%
    st_simplify(dTolerance = simplify_tolerance, preserveTopology = TRUE)
  
  #sf_use_s2(FALSE)
  # Lakes
  lakes <- ne_download(scale = 10, type = "lakes", category = "physical",
                       returnclass = "sf") %>%
    st_make_valid() %>%
    st_crop(bbox_poly)
    #st_make_valid()
  
  # Rivers
  rivers <- ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical",
                        returnclass = "sf") %>%
    st_make_valid() %>%
    st_crop(bbox_poly)
    #st_make_valid()
  
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
    ) %>% crop(bbox_poly) %>% terra::rast()
    
    elev[elev[] < 0] <- NA
    names(elev) <- "elevation"
    
    terra::writeRaster(elev, cached_elev, overwrite = TRUE)
  }
  
  return(list(
    countries = countries,
    lakes = lakes,
    rivers = rivers,
    altitude = elev,
    bbox = bbox_poly
  ))
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
adaptive_bbox <- function(pop_sf, padding = 0.12) {
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
    tm_raster(col='elevation', col.scale = tm_scale_continuous(values='brewer.greys', value.na='lightblue'), col_alpha = 0.9, col.legend=tm_legend(show=FALSE)) +
  tm_shape(spatial_data$countries) + 
    tm_borders(lwd = 1.2) +
  tm_shape(spatial_data$lakes) + 
    tm_polygons(col="darkblue", fill = "lightblue", fill_alpha = 0.9) +
  tm_shape(spatial_data$rivers) + 
    tm_lines(col="darkblue", lwd=1.2) +
  tm_layout(
    legend.show = FALSE
  ) +
  tm_compass(position = c("right", "top")) +
  tm_scalebar(position = c("left", "bottom"))
}

generate_spatial_map <- function(admix_df, popmap, popcoords, k_value) {
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
      ggplotGrob(population_pie_chart(proportions, colors = genetic_palette(k_value)))
    }),
    pop_admix$pop
  )
  
  return(list(pop_sf = pop_sf, pie_grobs = pie_grobs))
}