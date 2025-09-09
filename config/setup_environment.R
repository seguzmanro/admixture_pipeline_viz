#' Environment Setup Script
#' 
#' This script helps configure the environment for first-time users.
#' Run interactively to set up your project environment.

# Load required package
if (!requireNamespace("usethis", quietly = TRUE)) {
  install.packages("usethis")
}
library(usethis)

setup_environment <- function() {
  # Get project home directory
  project_home <- getwd()
  
  # Create .paxil_env file if doesn't exist
  env_file <- file.path(project_home, ".paxil_env")
  if (!file.exists(env_file)) {
    # Copy template to .paxil_env
    template_path <- system.file("config/environment_template.R", package = "yourpackage")
    if (template_path == "") {
      template_path <- "config/environment_template.R"
    }
    file.copy(template_path, env_file)
    
    # Open for editing
    usethis::edit_file(env_file)
    
    message("Created .paxil_env file. Please update the paths in this file.")
  } else {
    message(".paxil_env file already exists at: ", env_file)
  }
  
  # Set PAXIL_HOME in .Renviron
  renviron_path <- usethis::scoped_path_r(c("user", "project"), ".Renviron", envvar = "R_ENVIRON_USER")
  if (!file.exists(renviron_path)) file.create(renviron_path)
  
  current_env <- readLines(renviron_path)
  if (!any(grepl("^PAXIL_HOME=", current_env))) {
    write(paste0("PAXIL_HOME='", normalizePath(project_home), "'"), 
          renviron_path, append = TRUE)
    message("Added PAXIL_HOME to ", renviron_path)
  } else {
    message("PAXIL_HOME already set in ", renviron_path)
  }
  
  # Instructions for user
  message("\nSETUP COMPLETE. NEXT STEPS:")
  message("1. Edit the .paxil_env file with your actual paths")
  message("2. Restart R session for changes to take effect")
  message("3. Verify setup with: source('config/verify_environment.R')")
}

# Run setup interactively
if (interactive()) {
  setup_environment()
}