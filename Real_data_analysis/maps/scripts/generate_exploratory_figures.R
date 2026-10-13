required_packages <- c("ggplot2", "maps", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install the required packages: ", paste(missing_packages, collapse = ", "))
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument)) {
  script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]), winslash = "/")
  project_dir <- dirname(script_path)
} else {
  project_dir <- normalizePath(getwd(), winslash = "/")
}
while (!file.exists(file.path(project_dir, "R", "project_paths.R"))) {
  parent <- dirname(project_dir)
  if (identical(parent, project_dir)) stop("Could not locate the manuscript project.")
  project_dir <- parent
}

fit_path <- file.path(
  project_dir,
  "Real_data_analysis", "model_comparison", "results",
  "production_30000_unique_locations", "dcs_sppm", "pilot_fit.rds"
)
figure_dir <- file.path(project_dir, "Luc_i_a", "fig")
result_dir_override <- Sys.getenv("DCS_GRID_RESULT_DIR", "")
if (nzchar(result_dir_override)) {
  fit_path <- file.path(result_dir_override, "pilot_fit.rds")
}
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(project_dir, "R", "manuscript_figure_style.R"))

if (!file.exists(fit_path)) {
  stop("Final analysis fit not found: ", fit_path)
}

fit <- readRDS(fit_path)
required_fit_objects <- c(
  "station_information", "analysis_dates", "Y_imputed_mean", "data_scaling"
)
missing_fit_objects <- setdiff(required_fit_objects, names(fit))
if (length(missing_fit_objects)) {
  stop("The final fit is missing: ", paste(missing_fit_objects, collapse = ", "))
}

stations <- fit$station_information
dates <- as.Date(fit$analysis_dates)
response_standardized <- fit$Y_imputed_mean

if (!identical(dim(response_standardized), c(length(dates), nrow(stations)))) {
  stop("The response matrix does not align with the retained dates and stations.")
}

response <- response_standardized * fit$data_scaling$response_scale +
  fit$data_scaling$response_center

lag1_correlation <- function(values) {
  stats::cor(values[-length(values)], values[-1], use = "complete.obs")
}

stations$annual_mean_pm10 <- colMeans(response, na.rm = TRUE)
stations$lag1_acf <- apply(response, 2L, lag1_correlation)

daily_summary <- data.frame(
  Date = dates,
  Mean = rowMeans(response, na.rm = TRUE),
  Q25 = apply(response, 1L, stats::quantile, probs = 0.25, na.rm = TRUE),
  Median = apply(response, 1L, stats::median, na.rm = TRUE),
  Q75 = apply(response, 1L, stats::quantile, probs = 0.75, na.rm = TRUE)
)

study_bbox <- unlist(fit$pilot_settings$study_bbox)
if (!all(c("xmin", "ymin", "xmax", "ymax") %in% names(study_bbox))) {
  study_bbox <- c(xmin = -11, ymin = 34, xmax = 46, ymax = 72)
}

world <- ggplot2::map_data("world")

base_map <- function() {
  ggplot2::ggplot() +
    manuscript_country_layer(world) +
    ggplot2::coord_quickmap(
      xlim = unname(study_bbox[c("xmin", "xmax")]),
      ylim = unname(study_bbox[c("ymin", "ymax")]),
      expand = FALSE
    ) +
    manuscript_map_theme()
}

panel_network <- base_map() +
  ggplot2::geom_point(
    data = stations,
    ggplot2::aes(x = Longitude, y = Latitude, colour = "Monitoring site"),
    shape = 16,
    size = 1.15,
    alpha = 0.9
  ) +
  ggplot2::scale_colour_manual(
    values = c("Monitoring site" = "grey15"),
    name = NULL
  )

panel_mean <- base_map() +
  ggplot2::geom_point(
    data = stations,
    ggplot2::aes(x = Longitude, y = Latitude, colour = annual_mean_pm10),
    shape = 16,
    size = 1.25,
    alpha = 0.9
  ) +
  ggplot2::scale_colour_viridis_c(
    option = "C",
    end = 0.94,
    name = "Annual mean PM10\n(micrograms m^-3)",
    guide = manuscript_colourbar()
  )

panel_acf <- base_map() +
  ggplot2::geom_point(
    data = stations,
    ggplot2::aes(x = Longitude, y = Latitude, colour = lag1_acf),
    shape = 16,
    size = 1.25,
    alpha = 0.9
  ) +
  ggplot2::scale_colour_viridis_c(
    option = "D",
    direction = -1,
    name = "Lag-1 ACF of\ndaily PM10",
    guide = manuscript_colourbar()
  )

save_manuscript_pair(
  panel_network, "europe_pm10_network", figure_dir,
  width = 7.4, height = 5.8
)
save_manuscript_pair(
  panel_mean, "europe_pm10_mean", figure_dir,
  width = 7.4, height = 5.8
)
save_manuscript_pair(
  panel_acf, "europe_pm10_lag1_acf", figure_dir,
  width = 7.4, height = 5.8
)

temporal_candidate <- ggplot2::ggplot(
  daily_summary,
  ggplot2::aes(x = Date)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = Q25, ymax = Q75),
    fill = "grey75",
    alpha = 0.65
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = Mean),
    colour = "grey10",
    linewidth = 0.45
  ) +
  ggplot2::labs(
    x = NULL,
    y = "PM10 (micrograms m^-3)"
  ) +
  ggplot2::scale_x_date(
    date_breaks = "2 months",
    labels = english_month_labels,
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::theme_minimal(base_size = 11, base_family = "sans") +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.25),
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
    axis.title.y = ggplot2::element_text(size = 10),
    axis.text = ggplot2::element_text(colour = "grey25")
  )

save_manuscript_pair(
  temporal_candidate, "europe_pm10_temporal_summary", figure_dir,
  width = 7.5, height = 3.2
)

audit_lines <- c(
  "Data source and construction audit",
  "==================================",
  paste0("Final fit: ", normalizePath(fit_path, winslash = "/")),
  paste0("Stations: fit$station_information (n = ", nrow(stations), ")"),
  paste0("Dates: fit$analysis_dates (", min(dates), " to ", max(dates), ")"),
  "Daily PM10: fit$Y_imputed_mean transformed back to physical units using",
  "fit$data_scaling$response_center and fit$data_scaling$response_scale.",
  "This matrix is the final ordered and gap-imputed response passed to the production sampler.",
  "Annual mean: column mean of the final daily PM10 matrix.",
  "Lag-1 ACF: Pearson correlation between consecutive values of each final daily station series;",
  "no additional seasonal adjustment or residualization is introduced.",
  "Daily candidate: cross-station mean and interquartile range at each analysis date.",
  "Final analysis workflow: Real_data_analysis/model_comparison/scripts/run_four_models.ps1.",
  "Shared map convention: R/manuscript_figure_style.R."
)
audit_dir <- file.path(
  dirname(fit_path)
)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(audit_lines, file.path(audit_dir, "europe_pm10_exploratory_audit.txt"))

message("Created three exploratory map panels and the temporal summary in: ",
        figure_dir)
