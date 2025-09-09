#' Environment Configuration Template for Spatial Genomics Pipeline
#' 
#' This file provides a template for setting up the required environment variables.
#' Users should create a copy named `.paxil_env` in their project home directory
#' and populate it with their local paths. The pipeline will automatically load
#' these variables when run.
#' 
#' Instructions:
#' 1. Copy this file to your project home directory as `.paxil_env`
#' 2. Replace the placeholder paths with your actual file paths
#' 3. Set the PAXIL_HOME environment variable to your project directory:
#'    - Linux/macOS: add `export PAXIL_HOME=/path/to/project` to ~/.bashrc
#'    - Windows: set as system environment variable
#' 4. The pipeline will load these variables automatically

# Project home directory (set as system environment variable)
# PAXIL_HOME will be set at system level, not here

# =====================
# COUNTRY BOUNDARY DATA
# =====================

# Argentina administrative boundaries
ARG_ADM="/path/to/argentina_adm_shapefile.shp"

# Brazil administrative boundaries
BRA_ADM="/path/to/brazil_adm_shapefile.shp"

# Uruguay administrative boundaries
URY_ADM="/path/to/uruguay_adm_shapefile.shp"

# =================
# WATER BODY DATA
# =================

# Uruguay water bodies
URY_WAT="/path/to/uruguay_water_shapefile.shp"

# Brazil water bodies
BRA_WAT="/path/to/brazil_water_shapefile.shp"

# ====================
# ELEVATION RASTER DATA
# ====================

# Brazil altitude raster
BRA_ALT="/path/to/brazil_elevation.tif"

# Uruguay altitude raster
URY_ALT="/path/to/uruguay_elevation.tif"

# ========================
# OPTIONAL - CUSTOM EXTENTS
# ========================
# Uncomment and modify if you want to override default map boundaries
# STUDY_REGION_XMIN=-58.5
# STUDY_REGION_XMAX=-50
# STUDY_REGION_YMIN=-35.5
# STUDY_REGION_YMAX=-29.5

# ========================
# EXAMPLE DATA PATHS
# ========================
# For demonstration/testing purposes only
# EXAMPLE_FAM="/path/to/example_data/plink.fam"
# EXAMPLE_POPMAP="/path/to/example_data/popmap.csv"
# EXAMPLE_COORDS="/path/to/example_data/coordinates.csv"
# EXAMPLE_ADMIX_DIR="/path/to/example_data/admixture_results/"