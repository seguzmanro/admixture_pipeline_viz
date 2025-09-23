#' Custom Plot Formatting Utilities for Genomic Visualizations
#' 
#' This module provides standardized formatting functions for creating
#' publication-ready plots in population genomics workflows.

library(ggplot2)
library(RColorBrewer)
library(ggtext)
library(cowplot)

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
  max_colors <- 40
  if (n > max_colors) {
    warning("Requested ", n, " colors but maximum is ", max_colors)
    n <- max_colors
  }

  # Extended color palette combining multiple RColorBrewer palettes
  # Similar to Python's approach with tab20b + tab20c
  # Ordered by contrast: high contrast first, lower contrast last
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
#' @param q_matrix Data frame with admixture proportions
#' @param popmap Population map data frame
#' @param k_value Current K value
#' @param sort_populations Whether to sort by population (logical)
#' @param title Custom title for the plot (optional)
#' @return Formatted ggplot object
format_admixture_barplot <- function(q_matrix, popmap, k_value, sort_populations = TRUE, title = NULL) {
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
  
  # Create automatic title if not provided
  if (is.null(title)) {
    plot_title <- paste("Admixture Proportions - K =", k_value)
  } else {
    plot_title <- title
  }

  # Create plot
  ggplot(long_data, aes(x = sample, y = proportion, fill = cluster)) +
    geom_col(width = 1) +
    scale_fill_manual(values = genetic_palette(k_value)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = plot_title,
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

#' Format model selection plot (CV error, AIC, BIC, etc.)
#'
#' @param data Data frame with K values and model selection criteria
#' @param title Custom title for the plot (optional)
#' @param method Method for selecting optimal K: 'min' (default) or 'diffNgroup'
#' @return Formatted ggplot object
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
    # diffNgroup finds the K where the difference is most significant
    if (nrow(data) < 2) {
      stop("diffNgroup method requires at least 2 K values")
    }

    # Calculate absolute differences between consecutive values
    diffs <- abs(diff(data[[y_col]]))
    # Find the K where the difference is maximized
    # Add 1 because diff reduces the length by 1
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
  ggplot(data, aes(x = K, y = !!sym(y_col))) +
    geom_line(color = "#1F77B4", linewidth = 1) +
    geom_point(color = "#1F77B4", size = 3) +
    geom_vline(xintercept = optimal_k, linetype = "dashed", color = "#D62728") +
    annotate(
      "text",
      x = optimal_k + 0.2,
      y = max(data[[y_col]]) * 0.95,
      label = paste("Optimal K =", optimal_k, "(", method, ")"),
      color = "#D62728",
      hjust = 0,
      size = 4
    ) +
    scale_x_continuous(breaks = unique(data$K)) +
    labs(
      title = title,
      x = "Number of Clusters (K)",
      y = y_col
    ) +
    theme_genomics()
}

#' Format cross-validation error plot (deprecated - use format_model_selection_plot)
#'
#' @param cv_data Data frame with CV errors
#' @param method Method for selecting optimal K: 'min' (default) or 'diffNgroup'
#' @return Formatted ggplot object
#' @keywords internal
format_cv_plot <- function(cv_data, method = "min") {
  warning("format_cv_plot is deprecated. Use format_model_selection_plot instead.")
  format_model_selection_plot(cv_data, "ADMIXTURE Cross-Validation Error", method)
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
#' Format DAPC scatter plot
#'
#' Creates a standardized ggplot version of the DAPC scatter plot with cumulative variance inset.
#'
#' @param dapc DAPC object from adegenet package
#' @param nGroups Number of genetic clusters/groups (ignored if popmap is provided)
#' @param popmap Population map data frame (optional). If provided, nGroups will be set to the number of unique populations.
#' @return ggplot object with inset
#'
format_dapc_scatterplot <- function(dapc, nGroups, popmap = NULL) {
  # Extract coordinates and group assignments
  ld_scores <- dapc$ind.coord
  grp <- dapc$grp

  # Prepare data for plotting
  plot_data <- data.frame(
    LD1 = ld_scores[, 1],
    LD2 = ld_scores[, 2],
    group = grp,
    sample = rownames(ld_scores)
  )

  # Handle popmap if provided
  if (!is.null(popmap)) {
    # Validate popmap structure
    if (!"indv" %in% colnames(popmap) || !"pop" %in% colnames(popmap)) {
      stop("popmap must contain 'indv' and 'pop' columns")
    }

    # Merge with popmap to get population information
    plot_data <- merge(plot_data, popmap, by.x = "sample", by.y = "indv", all.x = TRUE)

    # Check for missing population data
    missing_pops <- sum(is.na(plot_data$pop))
    if (missing_pops > 0) {
      warning(missing_pops, " samples have no population assignment in popmap")
    }

    # Determine number of groups from unique populations
    unique_pops <- unique(na.omit(plot_data$pop))
    actual_nGroups <- length(unique_pops)

    if (actual_nGroups == 0) {
      stop("No valid population assignments found in popmap")
    }

    # Override nGroups with actual number of populations
    nGroups <- actual_nGroups

    # Create group labels from population names
    plot_data$group <- factor(plot_data$pop, levels = unique_pops, labels = unique_pops)

    message("Using ", nGroups, " populations from popmap: ", paste(unique_pops, collapse = ", "))
  } else {
    # Use generic cluster labels if no popmap provided
    plot_data$group <- factor(grp, labels = paste("Cluster", 1:nGroups))
  }
  
  # Use standardized genetic palette
  cluster_colors <- genetic_palette(nGroups)
  
  # Calculate variance explained for LD axes (discriminant axes)
  if (length(dapc$eig) >= 2) {
    var_ld1 <- round(100 * dapc$eig[1] / sum(dapc$eig), 1)
    var_ld2 <- round(100 * dapc$eig[2] / sum(dapc$eig), 1)
  } else {
    var_ld1 <- 100
    var_ld2 <- 0
  }
  
  # Create main scatter plot
  p <- ggplot(plot_data, aes(x = LD1, y = LD2, color = group)) +
    geom_point(pch = 16, size = 1.5, alpha = 0.8) +
    scale_color_manual(values = cluster_colors) +
    labs(
      x = paste("LD1 (", var_ld1, "%)", sep = ""),
      y = paste("LD2 (", var_ld2, "%)", sep = ""),
      color = "Cluster"
    ) +
    theme_genomics() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.8, "lines"),
      legend.text = element_text(size = rel(0.8))
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  # Prepare inset data for cumulative PCA variance
  pca_eig <- dapc$pca.eig
  cum_var <- 100 * cumsum(pca_eig) / sum(pca_eig)
  inset_data <- data.frame(
    axis = 1:length(cum_var),
    cum_var = cum_var
  )
  
  # Create inset plot
  inset_p <- ggplot(inset_data, aes(x = axis, y = cum_var)) +
    geom_hline(yintercept = 100, color = "lightgrey", linetype = "dashed", alpha = 0.5) +
    geom_step(color = "black", linewidth = 1.5) +
    geom_point(color = "black", size = 0.5) +
    scale_x_continuous(breaks = pretty(1:length(cum_var), n = 5),
                       expand = expansion(mult = c(0.01, 0.05))) +
    scale_y_continuous(breaks = pretty(c(0, cum_var), n = 5), limits = c(0, 100),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = "PCA axis", y = "Cumulated variance (%)") +
    theme_genomics(base_size = 8) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey80"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(size = rel(0.8), color = "black"),
      axis.title = element_text(size = rel(1.0), color = "black"),
      legend.position = "none"
    )
  
  # Combine with inset using cowplot (position similar to original: bottom left, but adjusted for ggplot)
  p_with_inset <- ggdraw(p) +
    draw_plot(inset_p, x = 0.02, y = 0.02, width = 0.28, height = 0.28)
  
  return(p_with_inset)
}