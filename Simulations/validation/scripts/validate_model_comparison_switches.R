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

cache_dir <- file.path(
  dirname(validation_root), "tmp", "sourcecpp_model_switch_validation"
)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(
  file.path(
    validation_root, "src", "sampler_global_beta_single_w_prediction.cpp"
  ),
  cacheDir = cache_dir,
  rebuild = FALSE
)
Sys.setenv(OMP_NUM_THREADS = "1")

task <- clustering_recovery_tasks(1L)[[
  which(vapply(
    clustering_recovery_tasks(1L),
    function(item) item$scenario == "regular_strips__intermediate",
    logical(1L)
  ))[1L]
]]
generated <- simulate_prediction_experiment_data(
  task, n_side = 6L, n_times = 30L, n_neighbors = 3L
)
prediction <- generated$prediction

run_chain <- function(model_type, similarity = scale(generated$simulation$locations),
                      prediction_similarity = prediction$locations_similarity) {
  run_global_beta_validation_chain(
    simulation = generated$simulation,
    runs = 80L,
    burn = 40L,
    thin = 1L,
    beta_update_mode = "collapsed",
    sigma2_update_mode = "collapsed",
    tau2_update_mode = "pairwise_interwoven",
    spatial_update_mode = "collapsed",
    model_type = model_type,
    similarity_locations = similarity,
    prior_override = pilot_fixed_prior(),
    seed = 20900001L,
    overdispersed_initialization = TRUE,
    save_latent = FALSE,
    save_w = FALSE,
    show_progress = FALSE,
    sampler_extra_args = list(
      prediction_locations = prediction$locations,
      prediction_locations_similarity = prediction_similarity,
      prediction_neighbors = prediction$neighbors,
      prediction_X = prediction$X,
      prediction_thin = 4L,
      prediction_seed = 20950001L,
      prediction_truth_u = prediction$truth_u,
      prediction_truth_w = prediction$truth_w,
      prediction_truth_y = prediction$truth_y,
      prediction_true_cluster = as.integer(prediction$true_cluster),
      prediction_cluster_phi = prediction$cluster_phi,
      prediction_cluster_tau2 = prediction$cluster_tau2
    )
  )
}

default_fit <- run_chain(NULL)
if (default_fit$model_type != "dcs_sppm" ||
    default_fit$partition_prior != "spatial_ppm" ||
    default_fit$predictive_regime_rule != "spatial_similarity" ||
    !isTRUE(default_fit$uses_common_shock)) {
  stop("The omitted model switch does not resolve to the full DCS-sPPM.")
}

ppm_fit <- run_chain("ppm_common_shock")
if (ppm_fit$partition_prior != "ordinary_ppm" ||
    ppm_fit$predictive_regime_rule != "cluster_size_only" ||
    isTRUE(ppm_fit$uses_spatial_partition)) {
  stop("PPM retains a spatial-similarity contribution.")
}

k1_fit <- run_chain("k1_common_shock")
if (!all(k1_fit$k_clust == 1L) || !all(k1_fit$cluster == 1L) ||
    length(k1_fit$phistar[[1L]]) != 1L) {
  stop("K=1 model performed a partition update.")
}

no_common_shock <- run_chain("sppm_no_common_shock")
if (max(abs(no_common_shock$u_mean)) != 0 ||
    max(abs(no_common_shock$prediction$u$mean)) != 0 ||
    !all(is.na(no_common_shock$kappa)) ||
    !all(is.na(no_common_shock$sigma2_u))) {
  stop("No-common-shock model retained a spatial common-shock contribution.")
}

validation <- data.frame(
  check = c(
    "default_resolves_to_full_model",
    "ppm_similarity_disabled",
    "k1_partition_disabled",
    "no_common_shock_zero_u"
  ),
  passed = TRUE,
  stringsAsFactors = FALSE
)
output_dir <- file.path(
  validation_results_dir, "Simulations", "validation", "results"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  validation, file.path(output_dir, "model_switch_validation.csv"),
  row.names = FALSE
)
print(validation)
