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
source(file.path(validation_root, "R", "manuscript_figure_style.R"))
source(file.path(validation_root, "R", "model_comparison_utils.R"))
`%||%` <- function(left, right) if (is.null(left)) right else left

required <- c("Rcpp", "ggplot2", "data.table", "maps", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Install: ", paste(missing, collapse = ", "))

result_name <- Sys.getenv("DCS_GRID_RESULTS_NAME", "production_30000_unique_locations")
result_dir_override <- Sys.getenv("DCS_GRID_RESULT_DIR", "")
result_dir <- if (nzchar(result_dir_override)) {
  normalizePath(result_dir_override, winslash = "/", mustWork = TRUE)
} else {
  file.path(validation_results_dir, "Real_data_analysis", "model_comparison", "results", result_name, "dcs_sppm")
}
fit_file <- file.path(result_dir, "pilot_fit.rds")
vi_file <- file.path(result_dir, "maps", "clustering", "vi_partition.rds")
if (!file.exists(fit_file)) stop("Missing production fit: ", fit_file)
if (!file.exists(vi_file)) stop("Missing VI representative partition: ", vi_file)

fit <- readRDS(fit_file)
model <- validate_dcs_model_type(
  fit$model %||% fit$model_type %||% fit$comparison_settings$model %||%
    "dcs_sppm"
)
uses_partition <- dcs_model_uses(model, "uses_partition")
uses_spatial_partition <- dcs_model_uses(model, "uses_spatial_partition")
vi_object <- readRDS(vi_file)
vi_partition <- as.integer(vi_object$partition)
vi_partition <- match(vi_partition, unique(vi_partition))
allocation <- as.matrix(fit$cluster)
if (length(vi_partition) != ncol(allocation)) {
  stop("The VI partition does not align with the retained allocation chain.")
}
if (is.null(fit$prediction$grid) || is.null(fit$prediction$locations)) {
  stop("The production fit does not contain the retained prediction grid.")
}

stations <- data.table::as.data.table(fit$station_information)
grid <- data.table::as.data.table(fit$prediction$grid)
prediction_locations <- as.matrix(fit$prediction$locations)
response_scale <- fit$data_scaling$response_scale

latitude_reference <- mean(stations$Latitude) * pi / 180
observed_locations <- cbind(
  x_1000km = stations$Longitude * cos(latitude_reference) * 111.32 / 1000,
  y_1000km = stations$Latitude * 110.57 / 1000
)
similarity_center <- colMeans(observed_locations)
similarity_scale <- apply(observed_locations, 2L, stats::sd)
observed_similarity <- sweep(
  sweep(observed_locations, 2L, similarity_center, "-"),
  2L, similarity_scale, "/"
)
prediction_similarity <- sweep(
  sweep(prediction_locations, 2L, similarity_center, "-"),
  2L, similarity_scale, "/"
)

output_dir <- file.path(result_dir, "predictive_temporal_regime")
figure_dir <- if (model == "dcs_sppm") {
  file.path(validation_root, "Luc_i_a", "fig")
} else {
  file.path(result_dir, "fig")
}
cache_dir <- file.path(dirname(validation_root), "tmp", "sourcecpp_predictive_regime")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(output_dir, "predictive_temporal_regime_summary.rds")
reuse_summary <- tolower(Sys.getenv("DCS_REUSE_REGIME_SUMMARY", "true")) == "true"
draw_stride <- as.integer(Sys.getenv("DCS_REGIME_DRAW_STRIDE", "10"))
threads <- as.integer(Sys.getenv("DCS_REGIME_THREADS", "8"))
if (!is.finite(draw_stride) || draw_stride < 1L) stop("Invalid draw stride.")
if (!is.finite(threads) || threads < 1L) stop("Invalid thread count.")
draw_indices <- seq.int(1L, nrow(allocation), by = draw_stride)

if (reuse_summary && file.exists(summary_file)) {
  regime <- readRDS(summary_file)
  if (!identical(as.integer(regime$draw_indices), as.integer(draw_indices)) ||
      !identical(regime$model_type, model)) {
    message("Stored regime settings differ from the request; recomputing.")
    reuse_summary <- FALSE
  }
}

if (!reuse_summary || !file.exists(summary_file)) {
  if (tolower(Sys.getenv("DCS_REQUIRE_CACHED_MANUSCRIPT", "false")) == "true") {
    stop("A matching saved predictive-regime summary is required for manuscript-only rendering.")
  }
  Rcpp::sourceCpp(
    file.path(validation_root, "src", "predictive_regime_postprocess.cpp"),
    cacheDir = cache_dir,
    rebuild = FALSE,
    showOutput = FALSE
  )
  regime <- predictive_regime_summary_cpp(
    allocation = allocation,
    phi_star = fit$phistar,
    tau2_star = fit$tau2star,
    mass = as.numeric(fit$Mdp),
    observed_similarity = observed_similarity,
    prediction_similarity = prediction_similarity,
    vi_partition = vi_partition,
    draw_indices = draw_indices,
    k0 = 1,
    v0 = 4,
    a_phi = 1,
    b_phi = 1,
    a_tau = 2,
    b_tau = 0.35,
    threads = threads,
    model_type = model
  )
  saveRDS(regime, summary_file, compress = FALSE)
}

site_phi_sum <- numeric(ncol(allocation))
site_phi_second <- numeric(ncol(allocation))
site_tau_sum <- numeric(ncol(allocation))
site_tau_second <- numeric(ncol(allocation))
for (draw in draw_indices) {
  labels <- allocation[draw, ]
  draw_phi <- as.numeric(fit$phistar[[draw]])[labels]
  draw_tau <- sqrt(as.numeric(fit$tau2star[[draw]])[labels])
  site_phi_sum <- site_phi_sum + draw_phi
  site_phi_second <- site_phi_second + draw_phi^2
  site_tau_sum <- site_tau_sum + draw_tau
  site_tau_second <- site_tau_second + draw_tau^2
}
number_of_draws <- length(draw_indices)
stations[, `:=`(
  phi_mean = site_phi_sum / number_of_draws,
  phi_sd = sqrt(pmax(0, site_phi_second / number_of_draws -
                       (site_phi_sum / number_of_draws)^2)),
  tau_mean_original = response_scale * site_tau_sum / number_of_draws,
  tau_sd_original = response_scale * sqrt(pmax(
    0,
    site_tau_second / number_of_draws - (site_tau_sum / number_of_draws)^2
  )),
  vi_cluster = vi_partition
)]

grid[, `:=`(
  phi_mean = as.numeric(regime$phi_mean),
  phi_sd = as.numeric(regime$phi_sd),
  tau_mean_original = response_scale * as.numeric(regime$tau_mean),
  tau_sd_original = response_scale * as.numeric(regime$tau_sd),
  tau2_mean_standardized = as.numeric(regime$tau2_mean),
  dominant_vi_regime = as.integer(regime$dominant_regime),
  maximum_regime_probability = as.numeric(regime$maximum_probability),
  normalized_regime_entropy = as.numeric(regime$normalized_entropy),
  new_regime_probability = as.numeric(regime$new_regime_probability)
)]
vi_phi_summary <- stations[, .(
  dominant_regime_phi = mean(phi_mean)
), by = vi_cluster]
vi_phi_lookup <- numeric(max(vi_partition))
vi_phi_lookup[vi_phi_summary$vi_cluster] <- vi_phi_summary$dominant_regime_phi
stations <- vi_phi_summary[stations, on = "vi_cluster"]
grid[, dominant_regime_phi := ifelse(
  dominant_vi_regime <= length(vi_phi_lookup),
  vi_phi_lookup[pmin(dominant_vi_regime, length(vi_phi_lookup))],
  NA_real_
)]

online_phi <- as.numeric(fit$prediction$regime$phi$mean)
online_tau2 <- as.numeric(fit$prediction$regime$tau2$mean)
safe_cor <- function(first, second) {
  if (stats::sd(first) == 0 || stats::sd(second) == 0) return(NA_real_)
  stats::cor(first, second)
}
validation <- data.frame(
  quantity = c("phi mean", "tau2 mean"),
  correlation_with_online_sampler = c(
    safe_cor(grid$phi_mean, online_phi),
    safe_cor(grid$tau2_mean_standardized, online_tau2)
  ),
  rmse_against_online_sampler = c(
    sqrt(mean((grid$phi_mean - online_phi)^2)),
    sqrt(mean((grid$tau2_mean_standardized - online_tau2)^2))
  ),
  online_retained_draws = fit$prediction$retained_draws,
  posthoc_draws = number_of_draws
)

data.table::fwrite(grid, file.path(output_dir, "predictive_temporal_regime_by_grid.csv"))
data.table::fwrite(stations, file.path(output_dir, "observed_temporal_regime_by_station.csv"))
data.table::fwrite(validation, file.path(output_dir, "online_prediction_validation.csv"))
saveRDS(
  list(
    grid = as.data.frame(grid),
    stations = as.data.frame(stations),
    regime_probability = regime$regime_probability,
    vi_partition = vi_partition,
    draw_indices = draw_indices,
    validation = validation,
    projection = regime$projection,
    new_regime_category = ncol(regime$regime_probability)
  ),
  file.path(output_dir, "predictive_temporal_regime_archive.rds"),
  compress = FALSE
)

world <- ggplot2::map_data("world")
study_bbox <- unlist(fit$pilot_settings$study_bbox)
if (!all(c("xmin", "ymin", "xmax", "ymax") %in% names(study_bbox))) {
  study_bbox <- c(xmin = -11, ymin = 34, xmax = 46, ymax = 72)
}

surface_map <- function(grid_variable, station_variable, legend, limits,
                        transform = "identity",
                        mark_network_only = FALSE) {
  plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(long, lat, group = group),
      fill = manuscript_palette$land,
      colour = NA
    ) +
    ggplot2::geom_raster(
      data = grid,
      ggplot2::aes(Longitude, Latitude, fill = .data[[grid_variable]])
    ) +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(long, lat, group = group),
      fill = NA,
      colour = manuscript_palette$border,
      linewidth = 0.25
    )
  if (mark_network_only) {
    plot <- plot + ggplot2::geom_point(
      data = stations,
      ggplot2::aes(Longitude, Latitude),
      colour = "grey15", alpha = 0.35, size = 0.45
    )
  } else {
    plot <- plot + ggplot2::geom_point(
      data = stations,
      ggplot2::aes(Longitude, Latitude, fill = .data[[station_variable]]),
      shape = 21, colour = "white", stroke = 0.10, size = 0.72
    )
  }
  plot +
    ggplot2::scale_fill_viridis_c(
      name = legend,
      limits = limits,
      trans = transform,
      oob = scales::squish,
      na.value = "grey75",
      guide = manuscript_colourbar()
    ) +
    ggplot2::coord_quickmap(
      xlim = unname(study_bbox[c("xmin", "xmax")]),
      ylim = unname(study_bbox[c("ymin", "ymax")]),
      expand = FALSE
    ) +
    manuscript_map_theme()
}

phi_limits <- c(0, 1)
tau_limits <- as.numeric(stats::quantile(
  c(grid$tau_mean_original, stations$tau_mean_original), c(0.01, 0.99)
))
if (!is.finite(diff(tau_limits)) || diff(tau_limits) <= 0) {
  center <- tau_limits[1L]
  width <- max(abs(center) * 0.01, sqrt(.Machine$double.eps))
  tau_limits <- center + c(-width, width)
}
certainty_limits <- c(0, 1)

phi_plot <- surface_map(
  "phi_mean", "phi_mean",
  "Predictive mean\npersistence", phi_limits
)
tau_plot <- surface_map(
  "tau_mean_original", "tau_mean_original",
  "Predictive mean\ninnovation SD\n(micrograms m^-3)",
  tau_limits, transform = "log10"
)
certainty_plot <- surface_map(
  "maximum_regime_probability", "phi_mean",
  "Maximum probability", certainty_limits,
  mark_network_only = TRUE
)
entropy_plot <- surface_map(
  "normalized_regime_entropy", "phi_mean",
  "Normalized entropy", c(0, 1),
  mark_network_only = TRUE
)
dominant_plot <- surface_map(
  "dominant_regime_phi", "dominant_regime_phi",
  "Dominant-regime\nmean persistence", c(0, 1)
)

save_manuscript_pair(
  phi_plot, "europe_predictive_phi", figure_dir,
  width = 7.4, height = 5.8
)
save_manuscript_pair(
  tau_plot, "europe_predictive_tau", figure_dir,
  width = 7.4, height = 5.8
)
save_manuscript_pair(
  certainty_plot, "europe_predictive_regime_uncertainty", figure_dir,
  width = 7.4, height = 5.8
)
save_manuscript_pair(
  entropy_plot, "europe_predictive_regime_entropy", figure_dir,
  width = 7.4, height = 5.8
)
save_manuscript_pair(
  dominant_plot, "europe_predictive_regime_map", figure_dir,
  width = 7.4, height = 5.8
)

method_lines <- c(
  "Predictive temporal-regime post-processing",
  "===========================================",
  paste0("Production fit: ", normalizePath(fit_file, winslash = "/")),
  paste0("Structural model: ", model, " (", dcs_model_label(model), ")"),
  paste0("VI partition: ", normalizePath(vi_file, winslash = "/")),
  paste0("Prediction locations: ", nrow(grid)),
  paste0("Posterior draws used post hoc: ", number_of_draws,
         " systematic draws (stride ", draw_stride, ")"),
  "No MCMC rerun is required.",
  if (uses_spatial_partition) {
    paste(
      "For every retained draw and grid location, existing-cluster and singleton",
      "weights use the spatial cohesion-3 ratio with standardized coordinates,",
      "k0=1, v0=4, identity L0, and the sampled mass parameter."
    )
  } else if (uses_partition) {
    paste(
      "Predictive allocation uses ordinary PPM weights: occupied-cluster sizes and",
      "the sampled mass parameter, with no spatial-similarity contribution."
    )
  } else {
    "The predictive temporal regime is the single global regime by construction."
  },
  if (uses_partition) {
    paste(
      "The singleton contribution integrates phi and tau moments under Beta(1,1)",
      "and IG(2,0.35), respectively."
    )
  } else {
    "No singleton regime is available under the fixed K=1 specification."
  },
  if (uses_partition) {
    paste(
      "Label invariance: draw-specific cluster probability is projected onto fixed",
      "VI groups in proportion to the observed VI composition of that cluster;",
      "the singleton remains a separate category. No raw label is averaged."
    )
  } else {
    "Label projection is unnecessary because the temporal regime is fixed globally."
  },
  "Maximum regime probability is the largest posterior probability across the VI",
  "groups and the singleton category. Normalized entropy divides Shannon entropy by",
  "log(number of VI groups + 1).",
  paste0("Projection description: ", regime$projection)
)
writeLines(method_lines, file.path(output_dir, "METHOD.txt"))

print(validation)
message("Numerical summaries: ", output_dir)
message("Manuscript figures: ", figure_dir)
