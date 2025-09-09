# Petunia Admixture Spatial Visualization Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)](https://www.r-project.org/)

## Overview

This pipeline processes ADMIXTURE results and creates spatial visualizations of population structure in *Petunia* species from Brazilian and Uruguayan pampas. It integrates genetic admixture proportions with geographic data to produce:

- Cross-validation error plots for K selection
- Individual ancestry bar plots
- Population-level admixture maps overlaid on topography
- Geographic distribution of genetic clusters

![Example Output](examples/Pdepint/plot_results/admixture_map_K3.png)
*Spatial visualization of genetic structure (K=3) in Petunia  populations*

## Features

- **Spatial-Genomic Integration**: Combines genetic data with elevation models, water bodies, and country boundaries
- **Publication-Ready Visualizations**: Creates high-quality figures for scientific publications
- **Reproducible Workflow**: Configuration management and session tracking
- **Parallel Processing**: Efficient handling of large genomic datasets
- **Customizable**: Adaptable to other species and study regions

## Installation

### Prerequisites

- R (≥ 4.0.0)
- System dependencies for spatial packages:
  ```bash
  # Ubuntu/Debian
  sudo apt-get install libudunits2-dev libgdal-dev libgeos-dev libproj-dev
  
  # Fedora/CentOS
  sudo dnf install udunits2-devel gdal-devel geos-devel proj-devel
  ```
### Setup
- Clone the repository:
  ```bash
  git clone https://github.com/yourusername/petunia-admixture.git
  cd petunia-admixture
  ```
- Install required R packages:  
  ```bash
  Rscript -e "install.packages(c('argparse', 'dotenv', 'sf', 'terra', 'tmap', 'ggplot2', 'dplyr', 'tidyr', 'future.apply', 'purrr',   'RColorBrewer'))"
  ```
- Configure environment:
  ```bash
  # Set project home directory
  echo "PAXIL_HOME=$(pwd)" >> ~/.Renviron
  ```
- Copy and edit environment template
  ```bash
  cp config/environment_template.R .paxil_env
  nano .paxil_env  # Update paths to your spatial data
  ```
- Verify configuration:
  ```bash
  Rscript config/verify_environment.R
  ```
## Usage
### Basic Execution
  ```bash
  Rscript R/admixture_visualization.R \
    --input_dir path/to/admixture_results \
    --fam path/to/plink.fam \
    --popmap path/to/popmap.csv \
    --coords path/to/coordinates.csv \
    --output_dir results
    --dpi 300
  ```
### Command Line Options
|Argument|Description|Default|
|---|---|---|
|--input_dir|ADMIXTURE output directory|```$EXAMPLE_ADMIX_DIR```|
|--fam|PLINK .fam file|```$EXAMPLE_FAM```|
|--popmap|CSV mapping samples to populations|```$EXAMPLE_POPMAP```|
|--coords|CSV with population coordinates|```$EXAMPLE_COORDS```|
|--output_dir|Output directory|```$EXAMPLE_OUTPUT_DIR```|
|--parallel_workers|Number of parallel workers|4|
|--dpi|Resolution for the PNG plots (mainly, the maps)|300|

## Input Files
### ADMIXTURE Output Directory:
  ```.Q``` files (admixture proportions)

  ```.out``` files (log files with CV errors)
### PLINK .fam File:
```text
  FAM001 SAMPLE001 0 0 0 -9
  FAM002 SAMPLE002 0 0 0 -9
  ...
```
### Population Map (popmap.csv):
```csv
  indv,pop
  SAMPLE001,POP_A
  SAMPLE002,POP_B
```
### Population Coordinates (coordinates.csv):
```csv
  pop,lon,lat
  POP_A,-56.123,-32.456
  POP_B,-55.987,-33.210
```
## Output
The pipeline generates these outputs in the specified directory:
```
results/
  ├── cv_error_plot.pdf           # CV error plot with optimal K
  ├── admixture_map_K2.png        # Spatial map for K=2
  ├── admixture_map_K3.png        # Spatial map for K=3
  ├── admixture_barplot_K2.pdf    # Admixture proportions barplot for K=2
  ├── admixture_barplot_K2.pdf    # Admixture proportions barplot for K=3
  └── session_info.txt            # Reproducibility information
```

## Example Dataset
Run with included sample data:
```bash
  Rscript R/admixture_visualization.R \
    --input_dir examples/Pdepint \
    --fam examples/Pdepint/Pdepint_M095_noSele.fam \
    --popmap examples/Pdepint/Pops_Pdepint.csv \
    --coords examples/Pdepint/Pdepint_Lat_Long_Coords.csv \
    --output_dir examples/Pdepint/plot_results
```

# Spatial Data Requirements
|Data Type|Description|Source|
|---|---|---|
|Country Boundaries|Administrative boundaries|GADM|
|Water Bodies|Lakes and rivers|Natural Earth|
|Elevation|Digital elevation models|SRTM|

# Customization
## For Other Species/Regions
- Modify bounding box in command line or code

- Update coordinate reference system in R/spatial_helpers.R

- Adjust color palettes in R/plot_formats.R

## Extending Functionality
- Add new visualizations in R/plot_formats.R

- Incorporate additional spatial layers in R/spatial_helpers.R

- Integrate with other genetics tools (PLINK, fastSTRUCTURE)

# Reproducibility
## The pipeline generates:

- Session information

- Package versions

- Execution parameters

- Environment configuration

# Citation
If using this software, please cite:
```bibtex
  @software{petunia_admixture_pipeline,
    author = {Sebastián Guzmán},
    title = {Petunia Admixture Spatial Visualization Pipeline},
    year = {2025},
    url = {https://github.com/seguzmanro/petunia-admixture}
  }
```

# License
MIT License - see LICENSE for details

# Support
For questions or issues:

1. Check Troubleshooting Guide

2. Open a GitHub Issue

---

<b>Developed by:</b> Sebastian Guzman

<b>Affiliation:</b> Departamento de Genetica - Universidade Federal do Rio Grande do Sul

<b>Contact:</b> sebastian.guzman@ufrgs.br