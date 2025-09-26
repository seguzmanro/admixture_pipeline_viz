#!/usr/bin/env Rscript
# Spatial Visualization of ADMIXTURE Results - CLI Interface
# This script provides a command-line interface to the admixspatial package

# Load the package
library(admixspatial)

# Execute main function if not running interactively
if (!interactive()) {
  main()
}