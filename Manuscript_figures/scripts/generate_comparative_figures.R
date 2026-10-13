script_arguments <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
validation_root <- if (length(script_arguments)) {
  dirname(normalizePath(sub("^--file=", "", script_arguments[[1L]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
repeat {
  if (file.exists(file.path(validation_root, "R", "project_paths.R"))) break
  parent <- dirname(validation_root)
  if (identical(parent, validation_root)) stop("Could not locate the manuscript project.")
  validation_root <- parent
}
source(file.path(validation_root, "R", "project_paths.R"))

required_packages <- c("ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}

simulation_name <- Sys.getenv("DCS_COMPARISON_RESULTS_NAME", "production")
real_data_name <- Sys.getenv(
  "DCS_REAL_FULL_RESULTS_NAME", "production_30000_unique_locations"
)
simulation_root <- file.path(
  validation_results_dir, "Simulations", "model_comparison", "results", simulation_name
)
real_data_root <- file.path(
  validation_results_dir, "Real_data_analysis", "model_comparison", "results", real_data_name
)
figure_dir <- Sys.getenv(
  "DCS_MANUSCRIPT_FIGURE_DIR", file.path(validation_root, "Luc_i_a", "fig")
)
asset_dir <- Sys.getenv(
  "DCS_MANUSCRIPT_ASSET_DIR", manuscript_assets_dir
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

theme_manuscript <- function() {
  ggplot2::theme_classic(base_size = 11, base_family = "sans") +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(1.1, "lines"),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
    )
}

save_figure <- function(plot, stem, width, height) {
  for (directory in c(figure_dir, asset_dir)) {
    ggplot2::ggsave(
      file.path(directory, paste0(stem, ".pdf")), plot,
      device = grDevices::cairo_pdf, width = width, height = height,
      units = "in", bg = "white"
    )
    ggplot2::ggsave(
      file.path(directory, paste0(stem, ".png")), plot,
      width = width, height = height, units = "in", dpi = 450, bg = "white"
    )
  }
}

clustering_file <- file.path(simulation_root, "clustering_replication_metrics.csv")
if (!file.exists(clustering_file)) stop("Missing: ", clustering_file)
clustering <- utils::read.csv(clustering_file, check.names = FALSE)
clustering <- clustering[
  clustering$model %in% c(
    "dcs_sppm", "ppm_common_shock", "sppm_no_common_shock"
  ) & is.finite(clustering$adjusted_rand_index),
]

model_labels <- c(
  dcs_sppm = "DCS-sPPM",
  ppm_common_shock = "PPM + common shock",
  sppm_no_common_shock = "sPPM without common shock"
)
geography_labels <- c(
  regular_strips = "Regular strips",
  irregular_regions = "Irregular regions",
  random_allocation = "Random allocation"
)
separation_labels <- c(
  very_similar = "Similar",
  intermediate = "Intermediate",
  very_different = "Different"
)

ari_summary <- stats::aggregate(
  adjusted_rand_index ~ model + geography + separation,
  clustering, mean
)
ari_summary$model <- factor(
  model_labels[ari_summary$model], levels = unname(model_labels)
)
ari_summary$geography <- factor(
  geography_labels[ari_summary$geography], levels = unname(geography_labels)
)
ari_summary$separation <- factor(
  separation_labels[ari_summary$separation], levels = unname(separation_labels)
)

simulation_plot <- ggplot2::ggplot(
  ari_summary,
  ggplot2::aes(
    x = separation, y = adjusted_rand_index,
    colour = model, linetype = model, shape = model, group = model
  )
) +
  ggplot2::geom_line(linewidth = 0.65) +
  ggplot2::geom_point(size = 2.4, stroke = 0.65) +
  ggplot2::facet_wrap(~ geography, nrow = 1L) +
  ggplot2::scale_colour_manual(values = c("#1B4F72", "#666666", "#B45F06")) +
  ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  ggplot2::scale_shape_manual(values = c(16, 17, 15)) +
  ggplot2::scale_y_continuous(
    limits = c(0.35, 1.01), breaks = seq(0.4, 1, by = 0.2),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(x = "Regime separation", y = "Mean adjusted Rand index") +
  theme_manuscript() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 9.5, angle = 45, hjust = 1, vjust = 1),
    strip.text.x = ggplot2::element_blank()
  )

save_figure(
  simulation_plot, "comparative_simulation_partition_recovery", 8.6, 3.85
)

surface_file <- file.path(real_data_root, "prediction_surface_difference_vs_dcs.csv")
parameter_file <- file.path(real_data_root, "posterior_parameter_comparison.csv")
regime_file <- file.path(real_data_root, "regime_parameter_comparison.csv")
for (path in c(surface_file, parameter_file, regime_file)) {
  if (!file.exists(path)) stop("Missing: ", path)
}

surface <- utils::read.csv(surface_file, check.names = FALSE)
surface <- surface[surface$process %in% c("u", "w", "y"), ]
surface <- rbind(
  data.frame(
    model = "dcs_sppm", model_label = "DCS-sPPM",
    process = c("u", "w", "y"), mean_surface_correlation = 1,
    mean_surface_rmse = 0, mean_surface_mae = 0,
    posterior_sd_rmse = 0, mean_sd_ratio = 1
  ),
  surface
)
surface$model <- factor(
  surface$model,
  levels = c(
    "dcs_sppm", "ppm_common_shock", "k1_common_shock",
    "sppm_no_common_shock"
  ),
  labels = c(
    "DCS-sPPM", "PPM + common shock", "K=1 + common shock",
    "sPPM without common shock"
  )
)
surface$process <- factor(
  surface$process, levels = c("u", "w", "y"),
  labels = c("Common shock\nu", "Latent state\nw", "Response\ny")
)

surface_plot <- ggplot2::ggplot(
  surface[is.finite(surface$mean_surface_correlation), ],
  ggplot2::aes(
    x = process, y = mean_surface_correlation,
    colour = model, shape = model, group = model
  )
) +
  ggplot2::geom_hline(yintercept = 1, colour = "grey75", linewidth = 0.35) +
  ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.28), linewidth = 0.55) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.28), size = 2.5) +
  ggplot2::scale_colour_manual(
    values = c("#1B4F72", "#666666", "#B45F06", "#7A3E9D")
  ) +
  ggplot2::scale_shape_manual(values = c(16, 17, 15, 18)) +
  ggplot2::scale_y_continuous(
    limits = c(-0.03, 1.02), breaks = seq(0, 1, by = 0.2),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    x = NULL, y = "Correlation with DCS-sPPM posterior mean"
  ) +
  theme_manuscript() +
  ggplot2::theme(legend.position = "none")

parameters <- utils::read.csv(parameter_file, check.names = FALSE)
regimes <- utils::read.csv(regime_file, check.names = FALSE)
variance_values <- rbind(
  data.frame(
    model = parameters$model[parameters$quantity == "measurement_variance"],
    component = "Measurement variance",
    value = parameters$mean[parameters$quantity == "measurement_variance"]
  ),
  data.frame(
    model = regimes$model,
    component = "Mean local innovation variance",
    value = regimes$tau2_mean
  ),
  data.frame(
    model = parameters$model[parameters$quantity == "common_shock_variance"],
    component = "Common-shock variance",
    value = parameters$mean[parameters$quantity == "common_shock_variance"]
  )
)
reference <- variance_values[variance_values$model == "dcs_sppm", c("component", "value")]
names(reference)[2L] <- "reference"
variance_values <- merge(variance_values, reference, by = "component")
variance_values$ratio <- variance_values$value / variance_values$reference
variance_values$model <- factor(
  variance_values$model,
  levels = c(
    "dcs_sppm", "ppm_common_shock", "k1_common_shock",
    "sppm_no_common_shock"
  ),
  labels = c(
    "DCS-sPPM", "PPM + common shock", "K=1 + common shock",
    "sPPM without common shock"
  )
)
variance_values$component <- factor(
  variance_values$component,
  levels = c(
    "Measurement variance", "Mean local innovation variance",
    "Common-shock variance"
  )
)

variance_plot <- ggplot2::ggplot(
  variance_values,
  ggplot2::aes(
    x = component, y = ratio, colour = model, shape = model, group = model
  )
) +
  ggplot2::geom_hline(yintercept = 1, colour = "grey45", linewidth = 0.4) +
  ggplot2::geom_point(
    position = ggplot2::position_dodge(width = 0.35), size = 2.7,
    stroke = 0.7
  ) +
  ggplot2::scale_colour_manual(
    values = c("#1B4F72", "#666666", "#B45F06", "#7A3E9D")
  ) +
  ggplot2::scale_shape_manual(values = c(16, 17, 15, 18)) +
  ggplot2::scale_y_log10(
    breaks = c(0.5, 1, 2, 5, 10),
    labels = c("0.5", "1", "2", "5", "10")
  ) +
  ggplot2::labs(x = NULL, y = "Posterior mean relative to DCS-sPPM") +
  theme_manuscript() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

real_data_plot <- patchwork::wrap_plots(
  surface_plot, variance_plot, widths = c(1, 1.15)
)

save_figure(
  real_data_plot, "comparative_real_data_structural_sensitivity", 9.2, 4.3
)

cat("Created comparative manuscript figures in:", figure_dir, "\n")
