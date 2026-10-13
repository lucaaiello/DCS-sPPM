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
source(file.path(validation_root, "R", "prediction_experiment_utils.R"))
source(file.path(validation_root, "R", "model_comparison_utils.R"))
source(file.path(validation_root, "R", "europe_model_comparison_utils.R"))

required_packages <- c("Rcpp", "Matrix")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1L), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}
invisible(loadNamespace("Matrix"))

csv_setting <- function(name, default) {
  values <- trimws(strsplit(Sys.getenv(name, default), ",", fixed = TRUE)[[1L]])
  values[nzchar(values)]
}
flag <- function(name, default) {
  tolower(Sys.getenv(name, if (default) "true" else "false")) %in%
    c("true", "1")
}
models <- csv_setting(
  "DCS_REAL_FULL_MODELS",
  "dcs_sppm,ppm_common_shock,k1_common_shock,sppm_no_common_shock"
)
models <- vapply(models, validate_dcs_model_type, character(1L))
runs <- as.integer(Sys.getenv("DCS_REAL_FULL_RUNS", "30000"))
burn <- as.integer(Sys.getenv("DCS_REAL_FULL_BURN", "20000"))
thin <- as.integer(Sys.getenv("DCS_REAL_FULL_THIN", "1"))
prediction_thin <- as.integer(Sys.getenv("DCS_REAL_FULL_PREDICTION_THIN", "1"))
threads <- as.integer(Sys.getenv("DCS_REAL_FULL_THREADS", "8"))
seed <- as.integer(Sys.getenv("DCS_REAL_FULL_SEED", "20260921"))
n_neighbors <- as.integer(Sys.getenv("DCS_REAL_FULL_NEIGHBORS", "10"))
reuse <- flag("DCS_REAL_FULL_REUSE", TRUE)
dry_run <- flag("DCS_REAL_FULL_DRY_RUN", FALSE)
preflight <- flag("DCS_REAL_FULL_PREFLIGHT", FALSE)
results_name <- Sys.getenv(
  "DCS_REAL_FULL_RESULTS_NAME", "production_30000_unique_locations"
)
if (!grepl("^[A-Za-z0-9_.-]+$", results_name)) {
  stop("DCS_REAL_FULL_RESULTS_NAME must be a safe directory name.")
}
if (runs <= burn || burn < 0L || thin < 1L || prediction_thin < 1L ||
    threads < 1L ||
    n_neighbors < 1L) stop("Invalid full-data settings.")

Sys.setenv(OMP_NUM_THREADS = threads)
if (!nzchar(Sys.getenv("DCS_REAL_FULL_DATA_FILE", ""))) {
  Sys.setenv(DCS_REAL_FULL_DATA_FILE = file.path(
    validation_root, "Data", "pm10_unique_locations_2024.rds"
  ))
}
output_root <- file.path(
  validation_results_dir, "Real_data_analysis", "model_comparison", "results", results_name
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(utils::sessionInfo()), file.path(output_root, "SESSION_INFO.txt"))
prepared <- prepare_europe_model_comparison_data(
  validation_root = validation_root,
  validation_results_dir = validation_results_dir,
  holdout_result_name = NULL,
  preparation_name = "production_30000_online_grid",
  n_neighbors = n_neighbors
)
prediction_inputs <- prepare_europe_full_prediction_inputs(
  prepared, n_neighbors = n_neighbors
)
cache_dir <- file.path(
  dirname(validation_root), "tmp", "sourcecpp_real_model_comparison_full"
)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(
  file.path(
    validation_root, "src",
    "sampler_global_beta_real_data_online_prediction.cpp"
  ),
  cacheDir = cache_dir,
  rebuild = FALSE
)
if (!isTRUE(validate_cholmod_runtime_cpp())) {
  stop("Matrix CHOLMOD runtime validation failed.")
}

if (dry_run || preflight) {
  for (model in models) {
    initial <- europe_initial_values(prepared, model)
    if (length(initial$phi) != ncol(prepared$Y) ||
        length(initial$tau2) != ncol(prepared$Y) ||
        any(!is.finite(c(initial$phi, initial$tau2)))) {
      stop("Invalid site-level initialization for model ", model, ".")
    }
  }
  cat(
    "Preflight completed: ", ncol(prepared$Y), " stations, ",
    nrow(prepared$Y), " days, ", nrow(prediction_inputs$grid),
    " prediction locations, and ", n_neighbors, " NNGP neighbors.\n",
    sep = ""
  )
  cat("No MCMC iteration was started.\n")
  quit(save = "no", status = 0L)
}

for (model in models) {
  model_dir <- file.path(output_root, model)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  result_file <- file.path(model_dir, "pilot_fit.rds")
  model_seed <- seed + 100000L * match(model, dcs_model_comparison_specs()$model)
  expected_settings <- list(
    runs = runs,
    burn = burn,
    thin = thin,
    prediction_thin = prediction_thin,
    threads = threads,
    seed = model_seed,
    n_neighbors = n_neighbors,
    source_support = "production_30000_online_grid",
    data_file = prepared$data_file,
    data_md5 = prepared$data_md5,
    training_locations = ncol(prepared$Y),
    prediction_locations = nrow(prediction_inputs$grid),
    model = model
  )
  writeLines(
    c(
      "European PM10 structural-comparison production fit",
      "=================================================",
      paste0("model=", model),
      paste0("runs=", runs),
      paste0("burn=", burn),
      paste0("thin=", thin),
      paste0("prediction_thin=", prediction_thin),
      paste0("threads=", threads),
      paste0("seed=", model_seed),
      paste0("n_neighbors=", n_neighbors),
      paste0("results_name=", results_name),
      "source_support=production_30000_online_grid",
      paste0("training_locations=", ncol(prepared$Y)),
      paste0("data_file=", prepared$data_file),
      paste0("data_md5=", prepared$data_md5),
      paste0("prediction_locations=", nrow(prediction_inputs$grid)),
      "sampler=src/sampler_global_beta_real_data_online_prediction.cpp"
    ),
    file.path(model_dir, "RUN_SETTINGS.txt")
  )
  if (reuse && file.exists(result_file)) {
    cached <- readRDS(result_file)
    if (identical(cached$comparison_settings, expected_settings)) {
      message("Reusing ", result_file)
      next
    }
    stop(
      "Existing result has incompatible settings: ", result_file,
      ". Choose another DCS_REAL_FULL_RESULTS_NAME or set reuse=false."
    )
  }
  initial <- europe_initial_values(prepared, model)
  set.seed(model_seed)
  elapsed <- system.time({
    fit <- dcs_sppm_global_beta_cpp(
      Y = prepared$Y,
      locations = prepared$locations,
      X = prepared$X,
      runs = runs, burn = burn, thin = thin, cohesion = 3L,
      init_beta = initial$beta,
      init_xi2 = initial$xi2,
      init_sigma2 = initial$sigma2,
      init_phi = initial$phi,
      init_tau2 = initial$tau2,
      init_M = initial$M,
      init_kappa = initial$kappa,
      nu = 0.5,
      NN_ind = prepared$NN_training,
      k0 = 1, v0 = 4,
      a_sigma = 3, b_sigma = 0.40,
      a_xi = 2, b_xi = 2,
      a_tau = 2, b_tau = 0.35,
      a_phi = 1, b_phi = 1,
      a_kappa = 3, b_kappa = 2,
      eta_phi = 0.05, eta_kappa = 0.05,
      CC = 3L,
      beta_update_mode = "collapsed",
      init_sigma2_u = initial$sigma2_u,
      a_sigma2_u = 3, b_sigma2_u = 0.60,
      locations_similarity = prepared$similarity_training,
      save_latent = FALSE, save_w = FALSE,
      sigma2_update_mode = "collapsed",
      tau2_update_mode = "pairwise_interwoven",
      joint_variance_initial_scale = 0.03,
      spatial_update_mode = "collapsed",
      spatial_initial_scale = 0.00003,
      spatial_process_model = "nngp",
      progress_file = file.path(model_dir, "progress.txt"),
      show_progress = TRUE,
      prediction_locations = prediction_inputs$locations,
      prediction_locations_similarity = prediction_inputs$similarity,
      prediction_neighbors = prediction_inputs$neighbors,
      prediction_X = prediction_inputs$X,
      prediction_thin = prediction_thin,
      prediction_seed = model_seed + 1L,
      model_type = model
    )
  })
  colnames(fit$beta) <- colnames(prepared$X)
  colnames(fit$sigma2beta) <- colnames(prepared$X)
  fit$model <- model
  fit$elapsed_seconds <- unname(elapsed["elapsed"])
  fit$station_information <- as.data.frame(prepared$station_information)
  fit$analysis_dates <- prepared$analysis_dates
  fit$data_scaling <- prepared$scaling
  fit$prediction$grid <- prediction_inputs$grid
  fit$prediction$nearest_observed_station <-
    prediction_inputs$nearest_observed_station
  fit$comparison_settings <- expected_settings
  fit$pilot_settings <- list(
    runs = runs,
    burn = burn,
    thin = thin,
    neighbors = n_neighbors,
    seed = model_seed,
    requested_openmp_threads = threads,
    joint_variance_initial_scale = 0.03,
    spatial_initial_scale = 0.00003,
    sampler_source = "sampler_global_beta_real_data_online_prediction.cpp",
    prediction_grid_n_lon = 100L,
    prediction_grid_n_lat = 100L,
    prediction_land_locations = nrow(prediction_inputs$grid),
    prediction_thin = prediction_thin,
    source_support = "production_30000_online_grid",
    model = model
  )
  saveRDS(fit, result_file, compress = FALSE)
  utils::write.csv(
    fit$block_timing, file.path(model_dir, "block_timing.csv"),
    row.names = FALSE
  )
  if (!is.null(fit$adaptation_diagnostics)) utils::write.csv(
    fit$adaptation_diagnostics,
    file.path(model_dir, "adaptation_diagnostics.csv"),
    row.names = FALSE
  )
  if (!is.null(fit$regime_block_timing)) utils::write.csv(
    fit$regime_block_timing,
    file.path(model_dir, "regime_block_timing.csv"),
    row.names = FALSE
  )
  if (!is.null(fit$pairwise_adaptation_by_cluster_size)) utils::write.csv(
    fit$pairwise_adaptation_by_cluster_size,
    file.path(model_dir, "pairwise_adaptation_by_cluster_size.csv"),
    row.names = FALSE
  )
  rm(fit)
  invisible(gc(FALSE))
}
cat("Full-data model-comparison fits completed in:", output_root, "\n")
