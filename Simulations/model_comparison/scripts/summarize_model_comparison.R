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
source(file.path(validation_root, "R", "model_comparison_utils.R"))

results_name <- Sys.getenv("DCS_COMPARISON_RESULTS_NAME", "production")
root <- file.path(
  validation_results_dir, "Simulations", "model_comparison", "results", results_name
)
dir.create(root, recursive = TRUE, showWarnings = FALSE)

bind_rows_compatible <- function(values) {
  values <- values[!vapply(values, is.null, logical(1L))]
  if (!length(values)) return(data.frame())
  columns <- unique(unlist(lapply(values, names), use.names = FALSE))
  values <- lapply(values, function(value) {
    for (column in setdiff(columns, names(value))) value[[column]] <- NA
    value[, columns, drop = FALSE]
  })
  do.call(rbind, values)
}

bind_component <- function(results, component) {
  bind_rows_compatible(lapply(results, `[[`, component))
}

clustering <- collect_model_comparison_results(file.path(root, "clustering"))
prediction <- collect_model_comparison_results(file.path(root, "prediction"))

read_legacy_full_results <- function(directory, kind) {
  files <- list.files(
    directory, pattern = "^metrics[.]rds$", recursive = TRUE,
    full.names = TRUE
  )
  lapply(files, function(file) {
    result <- readRDS(file)
    for (component in intersect(
      c("partition", "parameters", "prediction", "regime_prediction"),
      names(result)
    )) {
      if (is.null(result[[component]])) next
      result[[component]]$model <- "dcs_sppm"
      result[[component]]$model_label <- "DCS-sPPM"
      result[[component]] <- result[[component]][, c(
        "model", "model_label",
        setdiff(names(result[[component]]), c("model", "model_label"))
      )]
    }
    result
  })
}

if (!any(vapply(clustering, function(result) {
  identical(result$settings$model, "dcs_sppm")
}, logical(1L)))) {
  clustering <- c(
    read_legacy_full_results(
      file.path(simulation_clustering_results_dir, "clustering_recovery_100rep"),
      "clustering"
    ),
    clustering
  )
}
if (!any(vapply(prediction, function(result) {
  identical(result$settings$model, "dcs_sppm")
}, logical(1L)))) {
  prediction <- c(
    read_legacy_full_results(
      file.path(simulation_prediction_results_dir, "prediction_experiment_100rep"),
      "prediction"
    ),
    prediction
  )
}

partition_rows <- bind_component(clustering, "partition")
if (nrow(partition_rows)) {
  partition_summary <- summarize_metric_long(
    partition_rows,
    c("model", "model_label", "geography", "separation"),
    c(
      "adjusted_rand_index", "variation_of_information",
      "pairwise_sensitivity", "pairwise_specificity",
      "coclustering_brier", "exact_partition_recovery",
      "posterior_mean_clusters", "elapsed_seconds"
    )
  )
  utils::write.csv(
    partition_rows, file.path(root, "clustering_replication_metrics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    partition_summary, file.path(root, "simulation_comparison_table.csv"),
    row.names = FALSE
  )
}

parameter_rows <- bind_rows_compatible(list(
  bind_component(clustering, "parameters"),
  bind_component(prediction, "parameters")
))
if (nrow(parameter_rows)) {
  utils::write.csv(
    parameter_rows, file.path(root, "parameter_replication_metrics.csv"),
    row.names = FALSE
  )
}

prediction_rows <- bind_component(prediction, "prediction")
if (nrow(prediction_rows)) {
  prediction_summary <- summarize_metric_long(
    prediction_rows,
    c(
      "model", "model_label", "geography", "separation", "process",
      "stratum"
    ),
    c("rmse", "bias", "coverage", "interval_width", "crps", "elapsed_seconds")
  )
  utils::write.csv(
    prediction_rows, file.path(root, "prediction_replication_metrics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    prediction_summary, file.path(root, "prediction_comparison_table.csv"),
    row.names = FALSE
  )
}

regime_rows <- bind_component(prediction, "regime_prediction")
if (nrow(regime_rows)) {
  utils::write.csv(
    regime_rows, file.path(root, "predictive_regime_metrics.csv"),
    row.names = FALSE
  )
}

cat("Wrote model-comparison summaries to:", root, "\n")
