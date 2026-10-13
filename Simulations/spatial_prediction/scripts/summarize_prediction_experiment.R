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
source(file.path(validation_root, "R", "validation_utils.R"))
source(file.path(validation_root, "R", "pilot_simulation_utils.R"))
source(file.path(validation_root, "R", "full_simulation_utils.R"))
source(file.path(validation_root, "R", "clustering_recovery_utils.R"))

results_name <- Sys.getenv(
  "DCS_PREDICTION_RESULTS_NAME",
  "prediction_experiment_smoke"
)
results_dir <- file.path(simulation_prediction_results_dir, results_name)
manifest_file <- file.path(results_dir, "task_manifest.csv")
if (!file.exists(manifest_file)) stop("Task manifest not found: ", manifest_file)
manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)

metric_files <- list.files(
  results_dir,
  pattern = "^metrics\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)
results <- Filter(Negate(is.null), lapply(metric_files, function(file) {
  tryCatch(readRDS(file), error = function(error) NULL)
}))
if (!length(results)) stop("No completed prediction checkpoints were found.")

prediction_metrics <- do.call(rbind, lapply(results, `[[`, "prediction"))
regime_metrics <- do.call(rbind, lapply(results, `[[`, "regime_prediction"))
partition_metrics <- do.call(rbind, lapply(results, `[[`, "partition"))
parameter_metrics <- do.call(rbind, lapply(results, `[[`, "parameters"))
runtime_metrics <- do.call(rbind, lapply(results, function(result) {
  data.frame(
    scenario = result$settings$scenario,
    geography = result$settings$geography,
    separation = result$settings$separation,
    replicate = result$settings$replicate,
    elapsed_seconds = result$runtime$elapsed_seconds,
    prediction_seconds = result$runtime$prediction_seconds,
    prediction_percent = 100 * result$runtime$prediction_seconds /
      result$runtime$elapsed_seconds,
    prediction_draws = result$runtime$prediction_draws,
    prediction_memory_mb = result$runtime$prediction_memory_mb,
    openmp_threads = result$runtime$openmp_threads,
    stringsAsFactors = FALSE
  )
}))

summarize_prediction <- function(group) {
  data.frame(
    scenario = group$scenario[1],
    geography = group$geography[1],
    separation = group$separation[1],
    process = group$process[1],
    stratum = group$stratum[1],
    replications = nrow(group),
    bias = mean(group$bias),
    rmse = mean(group$rmse),
    rmse_mcse = stats::sd(group$rmse) / sqrt(nrow(group)),
    mae = mean(group$mae),
    coverage = mean(group$coverage),
    coverage_mcse = stats::sd(group$coverage) / sqrt(nrow(group)),
    interval_width = mean(group$interval_width),
    crps = mean(group$crps),
    mean_posterior_sd = mean(group$mean_posterior_sd),
    stringsAsFactors = FALSE
  )
}
prediction_summary <- do.call(rbind, lapply(
  split(
    prediction_metrics,
    list(
      prediction_metrics$scenario,
      prediction_metrics$process,
      prediction_metrics$stratum
    ),
    drop = TRUE
  ),
  summarize_prediction
))
rownames(prediction_summary) <- NULL

regime_summary <- do.call(rbind, lapply(
  split(regime_metrics, list(regime_metrics$scenario, regime_metrics$stratum)),
  function(group) {
    data.frame(
      scenario = group$scenario[1],
      geography = group$geography[1],
      separation = group$separation[1],
      stratum = group$stratum[1],
      replications = nrow(group),
      classification_accuracy = mean(group$classification_accuracy),
      accuracy_mcse = stats::sd(group$classification_accuracy) / sqrt(nrow(group)),
      mean_true_class_probability = mean(group$mean_true_class_probability),
      stringsAsFactors = FALSE
    )
  }
))
rownames(regime_summary) <- NULL

partition_summary <- summarize_full_partitions(partition_metrics)
parameter_summary <- summarize_full_parameters(parameter_metrics)

completed_keys <- unique(paste(
  prediction_metrics$scenario,
  prediction_metrics$replicate,
  sep = "::"
))
manifest$key <- paste(manifest$scenario, manifest$replicate, sep = "::")
completion <- do.call(rbind, lapply(split(manifest, manifest$scenario), function(group) {
  done <- group$key %in% completed_keys
  data.frame(
    scenario = group$scenario[1],
    geography = group$geography[1],
    separation = group$separation[1],
    expected = nrow(group),
    completed = sum(done),
    missing = sum(!done),
    stringsAsFactors = FALSE
  )
}))
rownames(completion) <- NULL

failure_files <- list.files(
  results_dir,
  pattern = "^failure\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)
failures <- if (length(failure_files)) {
  do.call(rbind, lapply(failure_files, function(file) {
    failure <- readRDS(file)
    data.frame(
      scenario = failure$settings$scenario,
      replicate = failure$settings$replicate,
      message = failure$message,
      file = file,
      stringsAsFactors = FALSE
    )
  }))
} else {
  data.frame(
    scenario = character(), replicate = integer(),
    message = character(), file = character()
  )
}

outputs <- list(
  prediction_metrics = prediction_metrics,
  prediction_summary = prediction_summary,
  regime_prediction_metrics = regime_metrics,
  regime_prediction_summary = regime_summary,
  partition_metrics = partition_metrics,
  partition_summary = partition_summary,
  parameter_metrics = parameter_metrics,
  parameter_summary = parameter_summary,
  runtime_metrics = runtime_metrics,
  completion_by_scenario = completion,
  failures = failures
)
for (name in names(outputs)) {
  utils::write.csv(
    outputs[[name]],
    file.path(results_dir, paste0(name, ".csv")),
    row.names = FALSE
  )
}

print(completion)
cat("\nRuntime:\n")
print(runtime_metrics)
cat("\nPrediction summary:\n")
print(prediction_summary)
if (nrow(failures)) cat("\nFailures:", nrow(failures), "\n")
