#' Plot Formatting Utilities for Genomic Visualizations
#'
#' This module provides standardized formatting functions for creating
#' publication-ready plots in population genomics workflows.

#' Theme for scientific publications
#'
#' @param base_size Base font size
#' @param base_family Base font family
#' @return Custom ggplot2 theme
#' @export
theme_genomics <- function(base_size = 11, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Core elements
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.3, color = "grey90"),
      panel.background = ggplot2::element_blank(),

      # Text elements
      plot.title = ggplot2::element_text(
        size = grid::unit(1.2, "lines"),
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 10)
      ),
      plot.subtitle = ggplot2::element_text(
        size = grid::unit(1.0, "lines"),
        hjust = 0.5,
        color = "grey30"
      ),
      axis.title = ggplot2::element_text(
        size = grid::unit(0.9, "lines"),
        face = "bold"
      ),
      axis.text = ggplot2::element_text(size = grid::unit(0.8, "lines")),
      legend.title = ggplot2::element_text(size = grid::unit(0.85, "lines"), face = "bold"),
      legend.text = ggplot2::element_text(size = grid::unit(0.8, "lines")),

      # Facet elements
      strip.text = ggplot2::element_text(
        size = grid::unit(0.9, "lines"),
        face = "bold",
        margin = ggplot2::margin(5, 0, 5, 0)
      ),
      strip.background = ggplot2::element_rect(fill = "grey92", color = NA),

      # Positioning
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.6, "lines"),
      plot.margin = ggplot2::margin(15, 15, 15, 15)
    )
}

#' Color palette for genetic components
#'
#' Generates a color palette optimized for visualizing genetic admixture components.
#' The palette combines multiple color schemes to provide good contrast and
#' accessibility for different numbers of genetic clusters.
#'
#' @param n Number of colors needed (integer, 1-40).
#'   Maximum of 40 colors supported. If n > 40, returns first 40 colors.
#' @return Vector of hex color codes (character vector).
#' @details
#' The palette is constructed from multiple RColorBrewer palettes in order of
#' visual contrast:
#' \enumerate{
#'   \item{Set1 (8 colors): Highest contrast, most vibrant}
#'   \item{Set3 (12 colors): Good contrast, varied and distinct}
#'   \item{Set2 (8 colors): Muted but still distinct}
#'   \item{Paired (12 colors): Lower contrast pairs}
#' }
#' This ordering ensures that the most visually distinct colors are used first.
#' @examples
#' # Get colors for 5 genetic clusters
#' colors <- genetic_palette(5)
#' print(colors)
#' # [1] "#377EB8" "#FF7F00" "#4DAF4A" "#E41A1C" "#984EA3"
#'
#' # Use in custom plots
#' barplot(1:10, col = genetic_palette(10))
#'
#' # Check available colors
#' length(genetic_palette(40))  # Maximum supported
#' @export
genetic_palette <- function(n) {
  max_colors <- 40
  if (n > max_colors) {
    warning("Requested ", n, " colors but maximum is ", max_colors)
    n <- max_colors
  }

  # Extended color palette combining multiple RColorBrewer palettes
  palette_40 <- c(
    # Set1 (8 colors) - highest contrast, most vibrant
    "#377EB8", "#FF7F00", "#4DAF4A", "#E41A1C","#984EA3",
    "#FFFF33", "#A65628", "#F781BF",

    # Set3 (12 colors) - good contrast, varied and distinct
    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
    "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
    "#CCEBC5", "#FFED6F",

    # Set2 (8 colors) - muted but still distinct
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
    "#FFD92F", "#E5C494", "#B3B3B3",

    # Paired (12 colors) - lowest contrast (pairs of similar colors)
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
    "#FFFF99", "#B15928"
  )

  return(palette_40[1:n])
}

#' Format admixture bar plot
#'
#' Creates a publication-ready bar plot showing individual admixture proportions
#' for each genetic cluster. Individuals are sorted by population and cluster
#' dominance for better visualization.
#'
#' @param q_matrix Data frame with admixture proportions (numeric matrix).
#'   Rows represent individuals, columns represent genetic clusters.
#' @param popmap Population map data frame (data.frame).
#'   Two-column data frame with sample IDs and population assignments.
#' @param k_value Current K value (integer).
#'   Number of genetic clusters (columns in q_matrix).
#' @param sort_populations Whether to sort by population (logical, default: TRUE).
#'   If TRUE, sorts individuals by population and cluster dominance.
#' @param title Custom title for the plot (character, optional).
#'   If NULL, generates automatic title based on K value.
#' @return Formatted ggplot object ready for display or saving.
#' @details
#' The function creates a stacked bar plot where:
#' \itemize{
#'   \item{Each bar represents an individual}
#'   \item{Bar segments show proportion of ancestry from each genetic cluster}
#'   \item{Individuals are grouped by population (if sort_populations = TRUE)}
#'   \item{Colors are automatically assigned using the genetic_palette function}
#'   \item{Plot uses the genomics theme for publication-ready appearance}
#' }
#' @examples
#' \dontrun{
#' # Create admixture bar plot
#' q_data <- read_q_file("admixture.3.Q", sample_ids)
#' popmap <- read.csv("popmap.csv")
#'
#' bar_plot <- format_admixture_barplot(
#'   q_matrix = q_data,
#'   popmap = popmap,
#'   k_value = 3,
#'   title = "Population Structure - K = 3"
#' )
#'
#' # Save the plot
#' save_genomics_plot(bar_plot, "admixture_barplot.pdf", width = 10, height = 6)
#' }
#' @export
format_admixture_barplot <- function(q_matrix, popmap, k_value, sort_populations = TRUE, title = NULL) {
  # Prepare data from Q matrix
  data <- as.data.frame(q_matrix)
  data$sample <- rownames(data)

  # Merge with population data
  data <- merge(data, popmap, by.x = "sample", by.y = "indv")

  # Convert to long format
  long_data <- data |>
    tidyr::pivot_longer(
      cols = dplyr::starts_with("Cluster"),
      names_to = "cluster",
      values_to = "proportion"
    ) |>
    dplyr::mutate(cluster = factor(cluster, levels = paste0("Cluster", 1:k_value)))

  # Sort by population and cluster proportions
  if (sort_populations) {
    pop_order <- long_data |>
      dplyr::group_by(pop) |>
      dplyr::summarize(dominant_cluster = which.max(tapply(proportion, cluster, mean))) |>
      dplyr::arrange(dominant_cluster, pop) |>
      dplyr::pull(pop)

    long_data$pop <- factor(long_data$pop, levels = pop_order)

    sample_order <- long_data |>
      dplyr::group_by(sample) |>
      dplyr::arrange(desc(proportion)) |>
      dplyr::slice(1) |>
      dplyr::arrange(pop, cluster, desc(proportion)) |>
      dplyr::pull(sample)

    long_data$sample <- factor(long_data$sample, levels = sample_order)
  }

  # Create automatic title if not provided
  if (is.null(title)) {
    plot_title <- paste("Admixture Proportions - K =", k_value)
  } else {
    plot_title <- title
  }

  # Create plot
  ggplot2::ggplot(long_data, ggplot2::aes(x = sample, y = proportion, fill = cluster)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::scale_fill_manual(values = genetic_palette(k_value)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(
      title = plot_title,
      x = "Samples",
      y = "Ancestry Proportion",
      fill = "Genetic\nCluster"
    ) +
    theme_genomics() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = grid::unit(0.6, "lines")),
      panel.spacing = ggplot2::unit(0.1, "lines"),
      legend.position = "right"
    ) +
    if (sort_populations) {
      ggplot2::facet_grid(~ pop, scales = "free_x", space = "free", switch = "x")
    }
}

#' Format model selection plot (CV error, AIC, BIC, etc.)
#'
#' @param data Data frame with K values and model selection criteria
#' @param title Custom title for the plot (optional)
#' @param method Method for selecting optimal K: 'min' (default) or 'diffNgroup'
#' @return Formatted ggplot object
#' @export
format_model_selection_plot <- function(data, title = NULL, method = "min") {
  # Validate input data structure
  if (ncol(data) != 2) {
    stop("Input data must have exactly 2 columns")
  }

  if (!"K" %in% colnames(data)) {
    stop("First column must be named 'K'")
  }

  # Validate method parameter
  if (!method %in% c("min", "diffNgroup")) {
    stop("Method must be either 'min' or 'diffNgroup'")
  }

  # Get the column names
  k_col <- "K"
  y_col <- setdiff(colnames(data), k_col)

  # Determine optimal K based on method
  if (method == "diffNgroup") {
    # Calculate differences between consecutive K values
    if (nrow(data) < 2) {
      stop("diffNgroup method requires at least 2 K values")
    }

    # Calculate absolute differences between consecutive values
    diffs <- abs(diff(data[[y_col]]))
    # Find the K where the difference is maximized
    optimal_idx <- which.max(diffs) + 1
    optimal_k <- data$K[optimal_idx]

    # Create automatic title if not provided
    if (is.null(title)) {
      title <- paste("Model Selection -", y_col, "(diffNgroup)")
    }
  } else {
    # Default minimum behavior
    optimal_k <- data$K[which.min(data[[y_col]])]

    # Create automatic title if not provided
    if (is.null(title)) {
      title <- paste("Model Selection -", y_col)
    }
  }

  # Create the plot
  ggplot2::ggplot(data, ggplot2::aes(x = K, y = !!rlang::sym(y_col))) +
    ggplot2::geom_line(color = "#1F77B4", linewidth = 1) +
    ggplot2::geom_point(color = "#1F77B4", size = 3) +
    ggplot2::geom_vline(xintercept = optimal_k, linetype = "dashed", color = "#D62728") +
    ggplot2::annotate(
      "text",
      x = optimal_k + 0.2,
      y = max(data[[y_col]]) * 0.95,
      label = paste("Optimal K =", optimal_k, "(", method, ")"),
      color = "#D62728",
      hjust = 0,
      size = 4
    ) +
    ggplot2::scale_x_continuous(breaks = unique(data$K)) +
    ggplot2::labs(
      title = title,
      x = "Number of Clusters (K)",
      y = y_col
    ) +
    theme_genomics()
}

#' Create a pie chart for population admixture
#'
#' @param proportions Named vector of admixture proportions
#' @param colors Vector of colors for clusters
#' @return ggplot object of pie chart
#' @export
population_pie_chart <- function(proportions, colors = NULL) {
  if (is.null(colors)) {
    colors <- genetic_palette(length(proportions))
  }

  data <- data.frame(
    cluster = names(proportions),
    proportion = proportions,
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(data, ggplot2::aes(x = "", y = proportion, fill = cluster)) +
    ggplot2::geom_col(width = 1, color = "white", linewidth = 0.2) +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}

#' Format spatial admixture map
#'
#' Creates a publication-ready spatial map showing population admixture proportions
#' as pie charts overlaid on a geographic basemap with elevation, water features,
#' and administrative boundaries.
#'
#' @param basemap Base map tmap object created by create_base_map.
#'   Should include elevation, boundaries, and water features.
#' @param q_matrix Data frame with admixture proportions (numeric matrix).
#'   Rows represent individuals, columns represent genetic clusters.
#' @param popmap Population map data frame (data.frame).
#'   Two-column data frame with sample IDs and population assignments.
#' @param popcoords Population coordinates data frame (data.frame).
#'   Three-column data frame with population names, longitude, and latitude.
#' @param k_value Number of genetic clusters (K) (integer, optional).
#'   If NULL, inferred from number of columns in q_matrix.
#' @param labels Add population labels to map (logical, default: FALSE).
#'   If TRUE, adds population names as text labels on the spatial map.
#' @param title Custom title for the map (character, optional).
#'   If NULL, generates automatic title based on K value.
#' @return tmap object ready for display or saving.
#' @details
#' The function creates a spatial admixture map where:
#' \itemize{
#'   \item{Each pie chart represents a population}
#'   \item{Pie segments show proportion of ancestry from each genetic cluster}
#'   \item{Base map includes elevation, water features, and administrative boundaries}
#'   \item{Population labels can be optionally displayed}
#'   \item{Colors are automatically assigned using the genetic_palette function}
#' }
#' @examples
#' \dontrun{
#' # Load required data
#' basemap <- create_base_map(spatial_data, study_bbox)
#' q_data <- read_q_file("admixture.3.Q", sample_ids)
#' popmap <- read.csv("popmap.csv")
#' coords <- read.csv("coordinates.csv")
#'
#' # Create spatial admixture map
#' spatial_map <- format_admixture_spatial(
#'   basemap = basemap,
#'   q_matrix = q_data,
#'   popmap = popmap,
#'   popcoords = coords,
#'   labels = TRUE,
#'   title = "Population Genetic Structure - K = 3"
#' )
#'
#' # Save the map
#' tmap::tmap_save(spatial_map, "admixture_map.png", width = 10, height = 8, dpi = 300)
#' }
#' @export
format_admixture_spatial <- function(basemap, q_matrix, popmap, popcoords, labels = FALSE, title = NULL) {

  # Determine k_value from Q matrix
  k_value <- ncol(q_matrix)
  
  # Generate spatial components
  spatial_components <- generate_spatial_map(
    admix_df = q_matrix,
    popmap = popmap,
    popcoords = popcoords,
    k_value = k_value
  )

  # Create automatic title if not provided
  if (is.null(title)) {
    plot_title <- paste("Genetic Structure - K =", k_value)
  } else {
    plot_title <- title
  }

  # Create the full map
  full_map <- basemap +
    tmap::tm_shape(spatial_components$pop_sf) +
    tmap::tm_symbols(
      shape = 'pop',
      shapes = spatial_components$pie_grobs,
      size = 0.7,
      border.lwd = 0.1,
      col_alpha = 0.5
    ) +
    tmap::tm_title(
      plot_title,
      position = c("center", "top")
    )

  # Add population labels if requested
  if (labels) {
    full_map <- full_map + tmap::tm_text("pop", size = 0.7, ymod = 1.7, col = "black", fontface = "bold")
  }

  return(full_map)
}

#' Save plot with standardized settings
#'
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution
#' @export
save_genomics_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ext <- tools::file_ext(filename)

  if (ext == "pdf") {
    ggplot2::ggsave(
      filename,
      plot,
      width = width,
      height = height,
      device = cairo_pdf
    )
  } else if (ext %in% c("png", "tiff", "jpeg")) {
    ggplot2::ggsave(
      filename,
      plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  } else if (ext == "svg") {
    ggplot2::ggsave(
      filename,
      plot,
      width = width,
      height = height
    )
  } else {
    stop("Unsupported file format: ", ext)
  }
}