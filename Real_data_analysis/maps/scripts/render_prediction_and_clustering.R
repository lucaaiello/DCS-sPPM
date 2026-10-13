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
source(file.path(validation_root, "R", "online_prediction_utils.R"))
source(file.path(validation_root, "R", "manuscript_figure_style.R"))
source(file.path(validation_root, "R", "model_comparison_utils.R"))
`%||%` <- function(left, right) if (is.null(left)) right else left

required <- c(
  "data.table", "ggplot2", "patchwork", "scales", "salso", "RANN", "maps"
)
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
if (!file.exists(fit_file)) stop("Missing fit: ", fit_file)
fit <- readRDS(fit_file)
model <- validate_dcs_model_type(fit$model %||% fit$model_type %||% "dcs_sppm")
uses_common_shock <- dcs_model_uses(model, "uses_common_shock")
uses_partition <- dcs_model_uses(model, "uses_partition")
prediction <- fit$prediction
if (is.null(prediction)) stop("The fit does not contain online predictions.")
grid <- data.table::as.data.table(prediction$grid)
dates <- as.Date(fit$analysis_dates)
world <- ggplot2::map_data("world")
if (nrow(grid) * length(dates) != length(prediction$y$mean)) {
  stop("Prediction grid and surface dimensions are inconsistent.")
}

map_dir <- file.path(result_dir, "maps")
predictive_dir <- file.path(map_dir, "predictive")
clustering_dir <- file.path(map_dir, "clustering")
figure_dir <- if (model == "dcs_sppm") {
  file.path(validation_root, "Luc_i_a", "fig")
} else {
  file.path(result_dir, "fig")
}
dir.create(predictive_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clustering_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

response_scale <- fit$data_scaling$response_scale
response_center <- fit$data_scaling$response_center
processes <- if (uses_common_shock) c("u", "w", "y") else c("w", "y")
z_values <- c(q005 = stats::qnorm(0.005), q025 = stats::qnorm(0.025),
              q975 = stats::qnorm(0.975), q995 = stats::qnorm(0.995))

scaled_summary <- function(process) {
  summary <- prediction[[process]]
  center <- if (process == "y") response_center else 0
  list(
    mean = summary$mean * response_scale + center,
    sd = summary$sd * response_scale
  )
}
surface_summaries <- setNames(lapply(processes, scaled_summary), processes)

robust_limits <- tolower(Sys.getenv("DCS_MAP_ROBUST_LIMITS", "true")) == "true"
range_limits <- function(values, nonnegative = FALSE) {
  values <- values[is.finite(values)]
  probabilities <- if (robust_limits) c(0.005, 0.995) else c(0, 1)
  limits <- as.numeric(stats::quantile(values, probabilities))
  if (nonnegative) limits[1L] <- 0
  if (!is.finite(diff(limits)) || diff(limits) <= 0) limits <- range(values)
  limits
}

plot_limits <- list()
render_predictive <- tolower(Sys.getenv(
  "DCS_RENDER_PREDICTIVE_MAPS",
  "true"
)) == "true"
for (process in processes) {
  summary <- surface_summaries[[process]]
  lower99 <- summary$mean + z_values["q005"] * summary$sd
  upper99 <- summary$mean + z_values["q995"] * summary$sd
  plot_limits[[process]] <- list(
    mean = range_limits(summary$mean),
    sd = range_limits(summary$sd, TRUE),
    interval = range_limits(c(lower99, upper99))
  )
}
limit_table <- data.table::rbindlist(lapply(processes, function(process) {
  data.table::rbindlist(lapply(names(plot_limits[[process]]), function(quantity) {
    data.table::data.table(
      process = process,
      quantity = quantity,
      lower = plot_limits[[process]][[quantity]][1L],
      upper = plot_limits[[process]][[quantity]][2L],
      robust = robust_limits
    )
  }))
}))
data.table::fwrite(limit_table, file.path(map_dir, "predictive_colour_limits.csv"))
writeLines(
  c(
    paste0("model=", model),
    paste0("uses_common_shock=", uses_common_shock),
    paste0("uses_partition=", uses_partition),
    if (!uses_common_shock) {
      "u maps omitted because the common shock is structurally zero."
    } else {
      "u maps applicable."
    },
    if (!uses_partition) {
      "VI clustering is the single global temporal regime by construction."
    } else {
      "VI clustering applicable."
    }
  ),
  file.path(map_dir, "MODEL_APPLICABILITY.txt")
)

surface_plot <- function(values, limits, legend_title) {
  plotting <- data.frame(grid, value = as.numeric(values))
  ggplot2::ggplot(plotting, ggplot2::aes(Longitude, Latitude, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(long, lat, group = group),
      inherit.aes = FALSE,
      fill = NA,
      colour = manuscript_palette$border,
      linewidth = 0.25,
      alpha = 0.9
    ) +
    ggplot2::scale_fill_viridis_c(
      limits = limits,
      oob = scales::squish,
      name = legend_title,
      guide = manuscript_colourbar()
    ) +
    ggplot2::coord_quickmap(
      xlim = c(-11, 46),
      ylim = c(34, 72),
      expand = FALSE
    ) +
    manuscript_map_theme(base_size = 10)
}

surface_panel_set <- function(process, time) {
  summary <- surface_summaries[[process]]
  mean_values <- summary$mean[time, ]
  sd_values <- summary$sd[time, ]
  list(
    surface_plot(mean_values, plot_limits[[process]]$mean,
                 paste0(process, " mean")),
    surface_plot(sd_values, plot_limits[[process]]$sd,
                 paste0(process, " SD")),
    surface_plot(mean_values + z_values["q025"] * sd_values,
                 plot_limits[[process]]$interval, process),
    surface_plot(mean_values + z_values["q975"] * sd_values,
                 plot_limits[[process]]$interval, process),
    surface_plot(mean_values + z_values["q005"] * sd_values,
                 plot_limits[[process]]$interval, process),
    surface_plot(mean_values + z_values["q995"] * sd_values,
                 plot_limits[[process]]$interval, process)
  )
}

if (render_predictive) for (process in processes) {
  process_dir <- file.path(predictive_dir, process)
  dir.create(process_dir, recursive = TRUE, showWarnings = FALSE)
  for (time in seq_along(dates)) {
    combined <- patchwork::wrap_plots(
      surface_panel_set(process, time), ncol = 3
    )
    ggplot2::ggsave(
      file.path(process_dir, sprintf("%s_%s.png", process, dates[time])),
      combined,
      width = 14,
      height = 8,
      dpi = 160,
      bg = "white"
    )
  }
}

manuscript_dates <- as.Date(c("2024-04-09", "2024-12-15"))
for (date_index in seq_along(manuscript_dates)) {
  date <- manuscript_dates[date_index]
  time <- match(date, dates)
  if (is.na(time)) stop("Requested manuscript date is absent: ", date)
  plot <- patchwork::wrap_plots(surface_panel_set("y", time)[1:4], ncol = 2)
  save_manuscript_pair(
    plot,
    stem = paste0("y_", date),
    directories = figure_dir,
    width = 10,
    height = 8,
    dpi = 320
  )
}

allocation <- as.matrix(fit$cluster)
vi_runs <- as.integer(Sys.getenv("DCS_VI_RUNS", "32"))
vi_cores <- as.integer(Sys.getenv("DCS_VI_CORES", "8"))
set.seed(as.integer(Sys.getenv("DCS_VI_SEED", "20260924")))
vi_file <- file.path(clustering_dir, "vi_partition.rds")
reuse_vi <- tolower(Sys.getenv("DCS_REUSE_VI", "true")) == "true"
require_cache <- tolower(Sys.getenv("DCS_REQUIRE_CACHED_MANUSCRIPT", "false")) == "true"
if (require_cache && uses_partition && (!reuse_vi || !file.exists(vi_file))) {
  stop("A saved VI partition is required for manuscript-only rendering.")
}
if (!uses_partition) {
  vi_partition <- rep.int(1L, ncol(allocation))
} else if (reuse_vi && file.exists(vi_file)) {
  vi_partition <- as.integer(readRDS(vi_file)$partition)
  if (length(vi_partition) != ncol(allocation)) {
    if (require_cache) stop("The saved VI partition does not match the final monitoring inventory.")
    message("Stored VI partition has the wrong number of stations; recomputing it.")
    vi_partition <- as.integer(salso::salso(
      allocation,
      loss = salso::VI(),
      nRuns = vi_runs,
      nCores = vi_cores
    ))
  }
} else {
  vi_partition <- as.integer(salso::salso(
    allocation,
    loss = salso::VI(),
    nRuns = vi_runs,
    nCores = vi_cores
  ))
}
vi_partition <- match(vi_partition, unique(vi_partition))
n_vi_clusters <- max(vi_partition)

site_phi_mean <- numeric(ncol(allocation))
site_tau2_mean <- numeric(ncol(allocation))
site_vi_cohesion <- numeric(ncol(allocation))
vi_members <- split(seq_len(ncol(allocation)), vi_partition)
cluster_phi_mean <- numeric(n_vi_clusters)
cluster_phi_m2 <- numeric(n_vi_clusters)
cluster_innovation_sd_mean <- numeric(n_vi_clusters)
cluster_innovation_sd_m2 <- numeric(n_vi_clusters)
for (draw in seq_len(nrow(allocation))) {
  labels <- allocation[draw, ]
  draw_phi <- as.numeric(fit$phistar[[draw]])[labels]
  draw_tau2 <- as.numeric(fit$tau2star[[draw]])[labels]
  site_phi_mean <- site_phi_mean + draw_phi
  site_tau2_mean <- site_tau2_mean + draw_tau2
  for (cluster in seq_along(vi_members)) {
    members <- vi_members[[cluster]]
    draw_cluster_phi <- mean(draw_phi[members])
    phi_delta <- draw_cluster_phi - cluster_phi_mean[cluster]
    cluster_phi_mean[cluster] <- cluster_phi_mean[cluster] + phi_delta / draw
    cluster_phi_m2[cluster] <- cluster_phi_m2[cluster] +
      phi_delta * (draw_cluster_phi - cluster_phi_mean[cluster])
    draw_cluster_innovation_sd <- mean(sqrt(draw_tau2[members]))
    tau_delta <- draw_cluster_innovation_sd -
      cluster_innovation_sd_mean[cluster]
    cluster_innovation_sd_mean[cluster] <-
      cluster_innovation_sd_mean[cluster] + tau_delta / draw
    cluster_innovation_sd_m2[cluster] <-
      cluster_innovation_sd_m2[cluster] + tau_delta *
      (draw_cluster_innovation_sd - cluster_innovation_sd_mean[cluster])
    if (length(members) == 1L) {
      site_vi_cohesion[members] <- site_vi_cohesion[members] + 1
    } else {
      counts <- table(labels[members])
      site_vi_cohesion[members] <- site_vi_cohesion[members] +
        (as.numeric(counts[as.character(labels[members])]) - 1) /
        (length(members) - 1)
    }
  }
}
site_phi_mean <- site_phi_mean / nrow(allocation)
site_tau2_mean <- site_tau2_mean / nrow(allocation)
site_vi_cohesion <- site_vi_cohesion / nrow(allocation)
cluster_phi_sd <- sqrt(cluster_phi_m2 / (nrow(allocation) - 1))
cluster_innovation_sd_sd <- sqrt(
  cluster_innovation_sd_m2 / (nrow(allocation) - 1)
)

stations <- data.table::as.data.table(fit$station_information)
stations[, `:=`(
  vi_cluster = vi_partition,
  posterior_mean_phi = site_phi_mean,
  posterior_mean_tau2 = site_tau2_mean,
  vi_within_cluster_coclustering = site_vi_cohesion,
  cluster_mean_phi = cluster_phi_mean[vi_partition],
  cluster_sd_phi = cluster_phi_sd[vi_partition],
  cluster_mean_innovation_sd = cluster_innovation_sd_mean[vi_partition],
  cluster_sd_innovation_sd = cluster_innovation_sd_sd[vi_partition],
  cluster_mean_innovation_sd_original =
    response_scale * cluster_innovation_sd_mean[vi_partition],
  cluster_sd_innovation_sd_original =
    response_scale * cluster_innovation_sd_sd[vi_partition]
)]
cluster_summary <- stations[, .(
  size = .N,
  longitude = mean(Longitude),
  latitude = mean(Latitude),
  cluster_mean_phi = unique(cluster_mean_phi),
  cluster_sd_phi = unique(cluster_sd_phi),
  cluster_mean_innovation_sd = unique(cluster_mean_innovation_sd),
  cluster_sd_innovation_sd = unique(cluster_sd_innovation_sd),
  cluster_mean_innovation_sd_original =
    unique(cluster_mean_innovation_sd_original),
  cluster_sd_innovation_sd_original =
    unique(cluster_sd_innovation_sd_original),
  mean_within_cluster_coclustering = mean(vi_within_cluster_coclustering)
), by = vi_cluster]

coordinates <- as.matrix(stations[, .(Longitude, Latitude)])
nearest <- RANN::nn2(coordinates, coordinates, k = min(7L, nrow(coordinates)))$nn.idx
adjacency <- vector("list", n_vi_clusters)
for (site in seq_len(nrow(stations))) {
  first <- vi_partition[site]
  for (neighbor in nearest[site, -1L]) {
    second <- vi_partition[neighbor]
    if (first != second) {
      adjacency[[first]] <- union(adjacency[[first]], second)
      adjacency[[second]] <- union(adjacency[[second]], first)
    }
  }
}
cluster_order <- order(lengths(adjacency), decreasing = TRUE)
graph_colour <- integer(n_vi_clusters)
for (cluster in cluster_order) {
  unavailable <- graph_colour[adjacency[[cluster]]]
  graph_colour[cluster] <- min(setdiff(seq_len(n_vi_clusters), unavailable))
}
n_graph_colours <- max(graph_colour)
palette <- grDevices::hcl.colors(n_graph_colours, "Dynamic")
stations[, graph_colour := graph_colour[vi_cluster]]
cluster_summary[, graph_colour := graph_colour[vi_cluster]]

data.table::fwrite(stations, file.path(clustering_dir, "vi_partition_by_station.csv"))
data.table::fwrite(cluster_summary, file.path(clustering_dir, "vi_cluster_summary.csv"))
saveRDS(
  list(partition = vi_partition, cluster_summary = cluster_summary,
       loss = "variation of information", salso_runs = vi_runs),
  file.path(clustering_dir, "vi_partition.rds"),
  compress = FALSE
)

cluster_map <- ggplot2::ggplot(
  stations,
  ggplot2::aes(Longitude, Latitude, colour = factor(graph_colour))
) +
  ggplot2::geom_polygon(
    data = world,
    ggplot2::aes(long, lat, group = group),
    inherit.aes = FALSE,
    fill = manuscript_palette$land,
    colour = manuscript_palette$border,
    linewidth = 0.25
  ) +
  ggplot2::geom_point(size = 1.15, alpha = 0.9) +
  ggplot2::scale_colour_manual(values = palette, guide = "none") +
  ggplot2::coord_quickmap(
    xlim = c(-11, 46),
    ylim = c(34, 72),
    expand = FALSE
  ) +
  manuscript_map_theme()

parameter_map <- function(variable, legend, limits = NULL,
                          encode_stability = TRUE, transform = "identity") {
  ggplot2::ggplot(stations, ggplot2::aes(
    Longitude,
    Latitude,
    colour = .data[[variable]],
    alpha = if (encode_stability) vi_within_cluster_coclustering else 1
  )) +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(long, lat, group = group),
      inherit.aes = FALSE,
      fill = manuscript_palette$land,
      colour = manuscript_palette$border,
      linewidth = 0.25
    ) +
    ggplot2::geom_point(size = 1.25) +
    ggplot2::scale_colour_viridis_c(
      name = legend,
      limits = limits,
      trans = transform,
      oob = scales::squish,
      guide = manuscript_colourbar()
    ) +
    ggplot2::scale_alpha_continuous(
      name = "Co-clustering",
      range = c(0.25, 1),
      limits = c(0, 1),
      guide = "none"
    ) +
    ggplot2::coord_quickmap(
      xlim = c(-11, 46),
      ylim = c(34, 72),
      expand = FALSE
    ) +
    manuscript_map_theme()
}

save_manuscript_pair(
  cluster_map, "vi_partition_graph_colours",
  c(clustering_dir, figure_dir), width = 7.4, height = 5.8
)
save_manuscript_pair(
  parameter_map("cluster_mean_phi", "Posterior mean\npersistence", limits = c(0, 1)),
  "vi_partition_phi", c(clustering_dir, figure_dir), width = 7.4, height = 5.8
)
save_manuscript_pair(
  parameter_map(
    "cluster_mean_innovation_sd_original",
    "Innovation SD\n(micrograms m^-3)",
    limits = response_scale * stats::quantile(
      cluster_innovation_sd_mean, c(0.01, 0.99)
    ),
    transform = "log10"
  ),
  "vi_partition_tau2", c(clustering_dir, figure_dir), width = 7.4, height = 5.8
)
save_manuscript_pair(
  parameter_map(
    "cluster_sd_phi", "Posterior SD of\npersistence", encode_stability = FALSE
  ),
  "vi_partition_phi_uncertainty", clustering_dir, width = 7.4, height = 5.8
)
save_manuscript_pair(
  parameter_map(
    "cluster_sd_innovation_sd_original",
    "Posterior SD of\ninnovation SD\n(micrograms m^-3)",
    limits = response_scale * stats::quantile(
      cluster_innovation_sd_sd, c(0.01, 0.99)
    ),
    encode_stability = FALSE,
    transform = "log10"
  ),
  "vi_partition_innovation_sd_uncertainty",
  clustering_dir, width = 7.4, height = 5.8
)
save_manuscript_pair(
  parameter_map("vi_within_cluster_coclustering", "Probability"),
  "vi_partition_stability", c(clustering_dir, figure_dir),
  width = 7.4, height = 5.8
)

representative_plot <- ggplot2::ggplot(
  cluster_summary,
  ggplot2::aes(
    cluster_mean_phi,
    cluster_mean_innovation_sd_original,
    size = size,
    colour = mean_within_cluster_coclustering
  )
) +
  ggplot2::geom_point(alpha = 0.8) +
  ggplot2::scale_colour_viridis_c(
    name = "Co-clustering",
    limits = c(0, 1)
  ) +
  ggplot2::scale_size_area(name = "Cluster size", max_size = 12) +
  ggplot2::labs(
    x = "Posterior mean persistence",
    y = "Innovation SD (micrograms m^-3)"
  ) +
  ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA)
  )
save_manuscript_pair(
  representative_plot, "vi_representative_parameter_summary",
  c(clustering_dir, figure_dir), width = 8, height = 6
)

cat("Rendered", if (render_predictive) length(dates) * length(processes) else 0,
    "predictive PNG files.\n")
cat("VI-optimal partition contains", n_vi_clusters, "clusters and uses",
    n_graph_colours, "graph colours.\n")
