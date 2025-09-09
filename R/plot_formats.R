#' Custom Plot Formatting Utilities for Genomic Visualizations
#' 
#' This module provides standardized formatting functions for creating
#' publication-ready plots in population genomics workflows.

library(ggplot2)
library(RColorBrewer)
library(ggtext)

#' Theme for scientific publications
#'
#' @param base_size Base font size
#' @param base_family Base font family
#' @return Custom ggplot2 theme
theme_genomics <- function(base_size = 11, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Core elements
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
      panel.background = element_blank(),
      
      # Text elements
      plot.title = element_text(
        size = rel(1.2),
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        size = rel(1.0),
        hjust = 0.5,
        color = "grey30"
      ),
      axis.title = element_text(
        size = rel(0.9),
        face = "bold"
      ),
      axis.text = element_text(size = rel(0.8)),
      legend.title = element_text(size = rel(0.85), face = "bold"),
      legend.text = element_text(size = rel(0.8)),
      
      # Facet elements
      strip.text = element_text(
        size = rel(0.9),
        face = "bold",
        margin = margin(5, 0, 5, 0)
      ),
      strip.background = element_rect(fill = "grey92", color = NA),
      
      # Positioning
      legend.position = "right",
      legend.key.size = unit(0.6, "lines"),
      plot.margin = margin(15, 15, 15, 15)
    )
}

#' Color palette for genetic components
#'
#' @param n Number of colors needed
#' @return Vector of color values
genetic_palette <- function(n) {
  max_colors <- 13
  if (n > max_colors) {
    warning("Requested ", n, " colors but maximum is ", max_colors)
    n <- max_colors
  }
  
  # Extended color-blind friendly palette
  base_palette <- c(
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
    "#AEC7E8", "#FFBB78"
  )
  
  return(base_palette[1:n])
}

#' Format admixture bar plot
#'
#' @param data Data frame with admixture proportions
#' @param k_value Current K value
#' @param sort_populations Whether to sort by population (logical)
#' @return Formatted ggplot object
format_admixture_barplot <- function(q_matrix, popmap, k_value, sort_populations = TRUE) {
  # Prepare data from Q matrix
  data <- as.data.frame(q_matrix)
  data$sample <- rownames(data)
  
  # Merge with population data
  data <- merge(data, popmap, by.x = "sample", by.y = "indv")
  
  # Convert to long format
  long_data <- data %>%
    pivot_longer(
      cols = starts_with("Cluster"),
      names_to = "cluster",
      values_to = "proportion"
    ) %>%
    mutate(cluster = factor(cluster, levels = paste0("Cluster", 1:k_value)))
  
  # Sort by population and cluster proportions
  if (sort_populations) {
    pop_order <- long_data %>%
      group_by(pop) %>%
      summarize(dominant_cluster = which.max(tapply(proportion, cluster, mean))) %>%
      arrange(dominant_cluster, pop) %>%
      pull(pop)
    
    long_data$pop <- factor(long_data$pop, levels = pop_order)
    
    sample_order <- long_data %>%
      group_by(sample) %>%
      arrange(desc(proportion)) %>%
      slice(1) %>%
      arrange(pop, cluster, desc(proportion)) %>%
      pull(sample)
    
    long_data$sample <- factor(long_data$sample, levels = sample_order)
  }
  
  # Create plot
  ggplot(long_data, aes(x = sample, y = proportion, fill = cluster)) +
    geom_col(width = 1) +
    scale_fill_manual(values = genetic_palette(k_value)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = paste("Admixture Proportions - K =", k_value),
      x = "Samples",
      y = "Ancestry Proportion",
      fill = "Genetic\nCluster"
    ) +
    theme_genomics() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(0.6)),
      panel.spacing = unit(0.1, "lines"),
      legend.position = "right"
    ) +
    if (sort_populations) {
      facet_grid(~ pop, scales = "free_x", space = "free", switch = "x")
    }
}

#' Format cross-validation error plot
#'
#' @param cv_data Data frame with CV errors
#' @return Formatted ggplot object
format_cv_plot <- function(cv_data) {
  min_k <- cv_data$K[which.min(cv_data$CV_error)]
  
  ggplot(cv_data, aes(x = K, y = CV_error)) +
    geom_line(color = "#1F77B4", linewidth = 1) +
    geom_point(color = "#1F77B4", size = 3) +
    geom_vline(xintercept = min_k, linetype = "dashed", color = "#D62728") +
    annotate(
      # aes(x = min_k + 0.2, y = max(CV_error) * 0.95),
      "text", x = min_k + 0.2, y = max(cv_data$CV_error) * 0.95,
      label = paste("Optimal K =", min_k),
      color = "#D62728",
      hjust = 0,
      size = 4
    ) +
    scale_x_continuous(breaks = unique(cv_data$K)) +
    labs(
      title = "ADMIXTURE Cross-Validation Error",
      x = "Number of Clusters (K)",
      y = "Cross-Validation Error"
    ) +
    theme_genomics()
}

#' Create a pie chart for population admixture
#'
#' @param proportions Named vector of admixture proportions
#' @param colors Vector of colors for clusters
#' @return ggplot object of pie chart
population_pie_chart <- function(proportions, colors = NULL) {
  if (is.null(colors)) {
    colors <- genetic_palette(length(proportions))
  }
  
  data <- data.frame(
    cluster = names(proportions),
    proportion = proportions,
    stringsAsFactors = FALSE
  )
  
  ggplot(data, aes(x = "", y = proportion, fill = cluster)) +
    geom_col(width = 1, color = "white", linewidth = 0.2) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = colors) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0)
    )
}

#' Save plot with standardized settings
#'
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution
save_genomics_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ext <- tools::file_ext(filename)
  
  if (ext == "pdf") {
    ggsave(
      filename,
      plot,
      width = width,
      height = height,
      device = cairo_pdf
    )
  } else if (ext %in% c("png", "tiff", "jpeg")) {
    ggsave(
      filename,
      plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  } else if (ext == "svg") {
    ggsave(
      filename,
      plot,
      width = width,
      height = height
    )
  } else {
    stop("Unsupported file format: ", ext)
  }
}

#' Create a PCA biplot with population coloring
#'
#' @param pca_data Data frame with PCA coordinates
#' @param pop_col Column containing population information
#' @param pc_x Principal component for x-axis (default: PC1)
#' @param pc_y Principal component for y-axis (default: PC2)
#' @param variance Vector with explained variance percentages
#' @return ggplot object
format_pca_biplot <- function(pca_data, pop_col = "pop", 
                             pc_x = "PC1", pc_y = "PC2",
                             variance = NULL) {
  if (!is.null(variance)) {
    x_lab <- paste0(pc_x, " (", round(variance[1], 1), "%)")
    y_lab <- paste0(pc_y, " (", round(variance[2], 1), "%)")
  } else {
    x_lab <- pc_x
    y_lab <- pc_y
  }
  
  ggplot(pca_data, aes_string(x = pc_x, y = pc_y, color = pop_col)) +
    geom_point(alpha = 0.8, size = 2.5) +
    stat_ellipse(level = 0.8, linewidth = 0.5) +
    scale_color_brewer(palette = "Set2") +
    labs(
      x = x_lab,
      y = y_lab,
      color = "Population"
    ) +
    theme_genomics() +
    theme(
      legend.position = "bottom",
      legend.box.spacing = unit(0, "cm")
    )
}

#' Create Manhattan plot for genome-wide association study
#'
#' @param gwas_data Data frame with GWAS results
#' @param chr_col Chromosome column name
#' @param pos_col Position column name
#' @param p_col P-value column name
#' @param sig_threshold Significance threshold
#' @return ggplot object
format_manhattan_plot <- function(gwas_data, 
                                  chr_col = "CHR", 
                                  pos_col = "BP", 
                                  p_col = "P",
                                  sig_threshold = 5e-8) {
  # Prepare data
  plot_data <- gwas_data %>%
    group_by(!!sym(chr_col)) %>%
    arrange(!!sym(chr_col), !!sym(pos_col)) %>%
    mutate(
      position_cum = cumsum(c(0, diff(!!sym(pos_col))) + 
        sum(lag(!!sym(pos_col), default = 0)),
      .before = !!sym(pos_col)
    )) %>%
    ungroup()
  
  # Axis settings
  axis_df <- plot_data %>%
    group_by(!!sym(chr_col)) %>%
    summarize(center = mean(position_cum))
  
  # Create plot
  ggplot(plot_data, aes(x = position_cum, y = -log10(!!sym(p_col)), 
                       color = factor(!!sym(chr_col) %% 2))) +
    geom_point(alpha = 0.7, size = 1) +
    geom_hline(yintercept = -log10(sig_threshold), 
               color = "red", 
               linetype = "dashed") +
    scale_x_continuous(
      breaks = axis_df$center,
      labels = axis_df[[chr_col]]
    ) +
    scale_color_manual(values = c("#1F77B4", "#AEC7E8")) +
    labs(
      x = "Chromosome",
      y = "-log<sub>10</sub>(p-value)"
    ) +
    theme_genomics() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      axis.title.y = element_markdown()
    )
}