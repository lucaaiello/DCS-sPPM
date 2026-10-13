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
source(file.path(validation_root, "R", "simulate_dcs_sppm.R"))
source(file.path(validation_root, "R", "sampler_runner.R"))
source(file.path(validation_root, "R", "validation_utils.R"))
source(file.path(validation_root, "R", "pilot_simulation_utils.R"))
source(file.path(validation_root, "R", "full_simulation_utils.R"))
source(file.path(validation_root, "R", "clustering_recovery_utils.R"))
source(file.path(validation_root, "R", "prediction_experiment_utils.R"))
source(file.path(validation_root, "R", "model_comparison_utils.R"))

flag <- function(name, default) {
  value <- tolower(Sys.getenv(name, if (default) "true" else "false"))
  if (!value %in% c("true", "false", "1", "0")) {
    stop(name, " must be true/false or 1/0.")
  }
  value %in% c("true", "1")
}
csv_setting <- function(name, default) {
  value <- trimws(strsplit(Sys.getenv(name, default), ",", fixed = TRUE)[[1L]])
  value[nzchar(value)]
}

smoke <- flag("DCS_COMPARISON_SMOKE", TRUE)
experiment <- match.arg(
  Sys.getenv("DCS_COMPARISON_EXPERIMENT", "both"),
  c("clustering", "prediction", "both")
)
models <- csv_setting(
  "DCS_COMPARISON_MODELS",
  if (smoke) {
    paste(dcs_model_comparison_specs()$model, collapse = ",")
  } else {
    "ppm_common_shock,k1_common_shock,sppm_no_common_shock"
  }
)
models <- vapply(models, validate_dcs_model_type, character(1L))
replications <- as.integer(Sys.getenv(
  "DCS_COMPARISON_REPLICATIONS", if (smoke) "1" else "100"
))
runs <- as.integer(Sys.getenv(
  "DCS_COMPARISON_RUNS", if (smoke) "200" else "8000"
))
burn <- as.integer(Sys.getenv(
  "DCS_COMPARISON_BURN", if (smoke) "100" else "4000"
))
thin <- as.integer(Sys.getenv("DCS_COMPARISON_THIN", "1"))
prediction_thin <- as.integer(Sys.getenv(
  "DCS_COMPARISON_PREDICTION_THIN", if (smoke) "5" else "10"
))
n_side <- as.integer(Sys.getenv(
  "DCS_COMPARISON_N_SIDE", if (smoke) "6" else "12"
))
clustering_times <- as.integer(Sys.getenv(
  "DCS_COMPARISON_CLUSTERING_TIMES", if (smoke) "30" else "365"
))
prediction_times <- as.integer(Sys.getenv(
  "DCS_COMPARISON_PREDICTION_TIMES", if (smoke) "30" else "100"
))
n_neighbors <- as.integer(Sys.getenv(
  "DCS_COMPARISON_NEIGHBORS", if (smoke) "3" else "10"
))
workers <- as.integer(Sys.getenv("DCS_COMPARISON_WORKERS", "1"))
omp_threads <- as.integer(Sys.getenv("DCS_COMPARISON_OMP_THREADS", "1"))
reuse <- flag("DCS_COMPARISON_REUSE", TRUE)
results_name <- Sys.getenv(
  "DCS_COMPARISON_RESULTS_NAME", if (smoke) "smoke" else "production"
)
if (!grepl("^[A-Za-z0-9_.-]+$", results_name)) {
  stop("DCS_COMPARISON_RESULTS_NAME must be a safe directory name.")
}

if (anyNA(c(
  replications, runs, burn, thin, prediction_thin, n_side,
  clustering_times, prediction_times, n_neighbors, workers, omp_threads
)) || replications < 1L || runs <= burn || burn < 0L || thin < 1L ||
    prediction_thin < 1L || n_side < 3L || (n_side^2) %% 3L != 0L ||
    clustering_times < 2L || prediction_times < 2L || n_neighbors < 1L ||
    workers < 1L || omp_threads < 1L) {
  stop("Invalid comparison settings.")
}

results_root <- file.path(
  validation_results_dir, "Simulations", "model_comparison", "results", results_name
)
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
cache_root <- file.path(dirname(validation_root), "tmp", "sourcecpp_model_comparisons")
dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
fixed_prior <- pilot_fixed_prior()
tasks <- comparison_tasks(replications, models)
scenario_filter <- csv_setting("DCS_COMPARISON_SCENARIOS", "")
if (length(scenario_filter)) {
  unknown_scenarios <- setdiff(
    scenario_filter, clustering_recovery_scenarios()$scenario
  )
  if (length(unknown_scenarios)) {
    stop("Unknown comparison scenario(s): ", paste(unknown_scenarios, collapse = ", "))
  }
  tasks <- tasks[vapply(
    tasks, function(task) task$scenario %in% scenario_filter, logical(1L)
  )]
}

format_duration <- function(seconds) {
  seconds <- max(0, round(seconds))
  sprintf(
    "%02d:%02d:%02d", seconds %/% 3600L,
    (seconds %% 3600L) %/% 60L, seconds %% 60L
  )
}

parallel_apply_with_progress <- function(cluster, tasks, fun, kind) {
  n_tasks <- length(tasks)
  n_workers <- length(cluster)
  results <- vector("list", n_tasks)
  started_at <- Sys.time()
  completed <- 0L
  failed <- 0L

  submit <- function(worker, index) {
    parallel:::sendCall(cluster[[worker]], fun, list(tasks[[index]]), tag = index)
  }

  for (index in seq_len(min(n_workers, n_tasks))) submit(index, index)

  while (completed < n_tasks) {
    received <- parallel:::recvOneResult(cluster)
    results[[received$tag]] <- received$value
    completed <- completed + 1L
    if (identical(received$value$status, "failed")) failed <- failed + 1L

    elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
    remaining <- elapsed / completed * (n_tasks - completed)
    cat(sprintf(
      "\r[%s] Progress: %5.1f%% | completed: %d/%d | elapsed: %s | eta: %s | failed: %d",
      kind, 100 * completed / n_tasks, completed, n_tasks,
      format_duration(elapsed), format_duration(remaining), failed
    ))
    flush.console()

    next_index <- completed + n_workers
    if (next_index <= n_tasks) submit(received$node, next_index)
  }
  cat("\n")
  results
}

run_experiment <- function(kind) {
  sampler_source <- if (kind == "clustering") {
    "sampler_global_beta_single_w.cpp"
  } else {
    "sampler_global_beta_single_w_prediction.cpp"
  }
  cpp_file <- file.path(validation_root, "src", sampler_source)
  output_dir <- file.path(results_root, kind)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  run_one <- function(task) {
    if (kind == "clustering") {
      run_clustering_comparison_task(
        task, output_dir, runs, burn, thin, n_side, clustering_times,
        n_neighbors, fixed_prior, reuse, sampler_source
      )
    } else {
      run_prediction_comparison_task(
        task, output_dir, runs, burn, thin, prediction_thin, n_side,
        prediction_times, n_neighbors, fixed_prior, reuse, sampler_source
      )
    }
  }

  if (workers == 1L) {
    Sys.setenv(OMP_NUM_THREADS = as.character(omp_threads))
    Rcpp::sourceCpp(
      cpp_file, cacheDir = file.path(cache_root, kind), rebuild = FALSE
    )
    status <- lapply(seq_along(tasks), function(index) {
      cat(sprintf("[%s %d/%d] %s\n", kind, index, length(tasks),
                  tasks[[index]]$display_name))
      tryCatch(
        run_one(tasks[[index]]),
        error = function(error) list(
          status = "failed", task = tasks[[index]]$display_name,
          message = conditionMessage(error)
        )
      )
    })
  } else {
    cluster <- parallel::makePSOCKcluster(min(workers, length(tasks)))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterExport(
      cluster,
      c(
        "validation_root", "cache_root", "kind", "cpp_file", "tasks",
        "output_dir", "runs", "burn", "thin", "prediction_thin",
        "n_side", "clustering_times", "prediction_times", "n_neighbors",
        "fixed_prior", "reuse", "sampler_source", "omp_threads"
      ),
      envir = environment()
    )
    parallel::clusterEvalQ(cluster, {
      Sys.setenv(OMP_NUM_THREADS = as.character(omp_threads))
      source(file.path(validation_root, "R", "project_paths.R"))
      source(file.path(validation_root, "R", "simulate_dcs_sppm.R"))
      source(file.path(validation_root, "R", "sampler_runner.R"))
      source(file.path(validation_root, "R", "validation_utils.R"))
      source(file.path(validation_root, "R", "pilot_simulation_utils.R"))
      source(file.path(validation_root, "R", "full_simulation_utils.R"))
      source(file.path(validation_root, "R", "clustering_recovery_utils.R"))
      source(file.path(validation_root, "R", "prediction_experiment_utils.R"))
      source(file.path(validation_root, "R", "model_comparison_utils.R"))
      worker_cache <- file.path(cache_root, kind, paste0("worker_", Sys.getpid()))
      dir.create(worker_cache, recursive = TRUE, showWarnings = FALSE)
      Rcpp::sourceCpp(cpp_file, cacheDir = worker_cache, rebuild = FALSE)
      NULL
    })
    status <- parallel_apply_with_progress(cluster, tasks, function(task) {
      tryCatch({
        if (kind == "clustering") {
          run_clustering_comparison_task(
            task, output_dir, runs, burn, thin, n_side, clustering_times,
            n_neighbors, fixed_prior, reuse, sampler_source
          )
        } else {
          run_prediction_comparison_task(
            task, output_dir, runs, burn, thin, prediction_thin, n_side,
            prediction_times, n_neighbors, fixed_prior, reuse, sampler_source
          )
        }
      }, error = function(error) list(
        status = "failed", task = task$display_name,
        message = conditionMessage(error)
      ))
    }, kind)
    parallel::stopCluster(cluster)
    on.exit(NULL, add = FALSE)
  }

  status_table <- do.call(rbind, lapply(seq_along(status), function(index) {
    item <- status[[index]]
    data.frame(
      model = tasks[[index]]$model,
      scenario = tasks[[index]]$scenario,
      replicate = tasks[[index]]$replicate,
      status = item$status %||% "completed",
      message = item$message %||% "",
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(
    status_table, file.path(output_dir, "run_status.csv"), row.names = FALSE
  )
  if (any(status_table$status == "failed")) {
    warning(sum(status_table$status == "failed"), " comparison tasks failed.")
  }
  invisible(status_table)
}

`%||%` <- function(left, right) if (is.null(left)) right else left

settings <- list(
  smoke = smoke, experiment = experiment, models = models,
  scenarios = unique(vapply(tasks, `[[`, character(1L), "scenario")),
  replications = replications, runs = runs, burn = burn, thin = thin,
  prediction_thin = prediction_thin, n_side = n_side,
  clustering_times = clustering_times, prediction_times = prediction_times,
  n_neighbors = n_neighbors, workers = workers, omp_threads = omp_threads,
  reuse = reuse, started_at = Sys.time()
)
saveRDS(settings, file.path(results_root, "run_settings.rds"))

if (experiment %in% c("clustering", "both")) run_experiment("clustering")
if (experiment %in% c("prediction", "both")) run_experiment("prediction")

cat("Model-comparison simulation jobs completed.\n")
