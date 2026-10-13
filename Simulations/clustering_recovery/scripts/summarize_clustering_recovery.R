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
  "DCS_CLUSTER_RESULTS_NAME",
  "clustering_recovery_100rep"
)
if (!grepl("^[A-Za-z0-9_.-]+$", results_name) ||
    results_name %in% c(".", "..")) {
  stop("`DCS_CLUSTER_RESULTS_NAME` must be a safe directory name.")
}
results_dir <- file.path(simulation_clustering_results_dir, results_name)
manifest_file <- file.path(results_dir, "task_manifest.csv")
if (!file.exists(manifest_file)) {
  stop("Task manifest not found: ", manifest_file)
}

manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
metric_files <- list.files(
  results_dir,
  pattern = "^metrics\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)
compact_results <- lapply(metric_files, function(file) {
  tryCatch(readRDS(file), error = function(error) NULL)
})
compact_results <- Filter(Negate(is.null), compact_results)

if (!length(compact_results)) {
  stop("No completed clustering-recovery checkpoints were found.")
}

partition_metrics <- do.call(
  rbind,
  lapply(compact_results, `[[`, "partition")
)
parameter_metrics <- do.call(
  rbind,
  lapply(compact_results, `[[`, "parameters")
)
partition_summary <- summarize_clustering_recovery_partitions(partition_metrics)
parameter_summary <- summarize_full_parameters(parameter_metrics)

completed_keys <- unique(paste(
  partition_metrics$scenario,
  partition_metrics$replicate,
  sep = "::"
))
manifest$key <- paste(manifest$scenario, manifest$replicate, sep = "::")
manifest$status <- ifelse(manifest$key %in% completed_keys, "completed", "missing")
completion_by_scenario <- do.call(rbind, lapply(
  split(manifest, manifest$scenario),
  function(group) {
    data.frame(
      scenario = group$scenario[1],
      geography = group$geography[1],
      separation = group$separation[1],
      expected = nrow(group),
      completed = sum(group$status == "completed"),
      missing = sum(group$status == "missing"),
      stringsAsFactors = FALSE
    )
  }
))
rownames(completion_by_scenario) <- NULL

failure_files <- list.files(
  results_dir,
  pattern = "^failure\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)
failures <- if (length(failure_files)) {
  do.call(rbind, lapply(failure_files, function(file) {
    failure <- tryCatch(readRDS(file), error = function(error) NULL)
    if (is.null(failure)) {
      return(NULL)
    }
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
    scenario = character(),
    replicate = integer(),
    message = character(),
    file = character(),
    stringsAsFactors = FALSE
  )
}

utils::write.csv(
  partition_metrics,
  file.path(results_dir, "partition_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  partition_summary,
  file.path(results_dir, "partition_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  parameter_metrics,
  file.path(results_dir, "parameter_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  parameter_summary,
  file.path(results_dir, "parameter_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  completion_by_scenario,
  file.path(results_dir, "completion_by_scenario.csv"),
  row.names = FALSE
)
utils::write.csv(
  failures,
  file.path(results_dir, "failures.csv"),
  row.names = FALSE
)

print(completion_by_scenario)
cat("\nPartition recovery summary:\n")
print(partition_summary)
if (nrow(failures)) {
  cat("\nFailures:", nrow(failures), "\n")
}
