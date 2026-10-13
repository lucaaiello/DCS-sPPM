find_global_beta_validation_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    candidate <- file.path(current, "DCS_sPPM_global_beta_sim_validation")
    if (file.exists(file.path(candidate, "README.md"))) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }

    if (file.exists(file.path(current, "README.md")) &&
        file.exists(file.path(current, "R", "project_paths.R")) &&
        file.exists(file.path(current, "src", "sampler_global_beta_single_w.cpp"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }

  stop("Could not locate `DCS_sPPM_global_beta_sim_validation`.")
}

validation_root <- find_global_beta_validation_root()
validation_results_dir <- validation_root
simulation_clustering_results_dir <- file.path(validation_root, "Simulations", "clustering_recovery", "results")
simulation_prediction_results_dir <- file.path(validation_root, "Simulations", "spatial_prediction", "results")
simulation_comparison_results_dir <- file.path(validation_root, "Simulations", "model_comparison", "results")
simulation_benchmark_results_dir <- file.path(validation_root, "Simulations", "gp_nngp_benchmark", "results")
simulation_validation_results_dir <- file.path(validation_root, "Simulations", "validation", "results")
real_data_comparison_results_dir <- file.path(validation_root, "Real_data_analysis", "model_comparison", "results")
real_data_reference_results_dir <- file.path(validation_root, "Real_data_analysis", "reference_inputs", "results")
real_data_diagnostics_results_dir <- file.path(validation_root, "Real_data_analysis", "mixing_diagnostics", "results")
real_data_preparation_results_dir <- file.path(validation_root, "Real_data_analysis", "preparation", "results")
manuscript_assets_dir <- file.path(validation_root, "Manuscript_figures", "results")

