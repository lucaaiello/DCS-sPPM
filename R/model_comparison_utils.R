dcs_model_comparison_specs <- function() {
  data.frame(
    model = c(
      "dcs_sppm", "ppm_common_shock", "k1_common_shock",
      "sppm_no_common_shock"
    ),
    label = c(
      "DCS-sPPM", "PPM + common shock", "K=1 + common shock",
      "sPPM without common shock"
    ),
    uses_partition = c(TRUE, TRUE, FALSE, TRUE),
    uses_spatial_partition = c(TRUE, FALSE, FALSE, TRUE),
    uses_common_shock = c(TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
}

validate_dcs_model_type <- function(model) {
  choices <- dcs_model_comparison_specs()$model
  if (length(model) != 1L || is.na(model) || !model %in% choices) {
    stop("Unknown model type: ", paste(model, collapse = ", "))
  }
  model
}

dcs_model_label <- function(model) {
  specs <- dcs_model_comparison_specs()
  specs$label[match(validate_dcs_model_type(model), specs$model)]
}

dcs_model_uses <- function(model, component) {
  specs <- dcs_model_comparison_specs()
  if (!component %in% c(
    "uses_partition", "uses_spatial_partition", "uses_common_shock"
  )) {
    stop("Unknown model component: ", component)
  }
  specs[[component]][match(validate_dcs_model_type(model), specs$model)]
}

comparison_tasks <- function(replications = 100L, models = NULL) {
  if (is.null(models)) models <- dcs_model_comparison_specs()$model
  models <- vapply(models, validate_dcs_model_type, character(1L))
  base_tasks <- clustering_recovery_tasks(replications)
  tasks <- vector("list", length(base_tasks) * length(models))
  index <- 0L
  for (base_task in base_tasks) {
    for (model in models) {
      index <- index + 1L
      task <- base_task
      task$model <- model
      task$model_label <- dcs_model_label(model)
      task$display_name <- paste(model, base_task$display_name, sep = "__")
      task$comparison_index <- index
      tasks[[index]] <- task
    }
  }
  tasks
}

comparison_chain_seed <- function(task, seed_base) {
  model_index <- match(task$model, dcs_model_comparison_specs()$model)
  as.integer(
    seed_base + 1000000L * model_index +
      10000L * task$scenario_index + task$replicate
  )
}

comparison_partition_metrics <- function(fit, simulation, task) {
  common <- data.frame(
    model = task$model,
    model_label = task$model_label,
    scenario = task$scenario,
    geography = task$geography,
    separation = task$separation,
    replicate = task$replicate,
    stringsAsFactors = FALSE
  )
  if (!dcs_model_uses(task$model, "uses_partition")) {
    return(cbind(
      common,
      adjusted_rand_index = NA_real_,
      variation_of_information = NA_real_,
      pairwise_sensitivity = NA_real_,
      pairwise_specificity = NA_real_,
      coclustering_brier = NA_real_,
      exact_partition_recovery = NA,
      true_clusters = length(unique(simulation$true$cluster)),
      posterior_mean_clusters = 1,
      posterior_mode_clusters = 1L,
      spatial_acceptance = fit$acc_spatial,
      phi_acceptance = fit$accphi,
      elapsed_seconds = fit$elapsed_seconds
    ))
  }
  result <- clustering_recovery_partition_metrics(fit, simulation, task)$metrics
  result$model <- task$model
  result$model_label <- task$model_label
  result[, c("model", "model_label", setdiff(names(result), c("model", "model_label")))]
}

comparison_parameter_metrics <- function(fit, simulation, task, experiment) {
  if (dcs_model_uses(task$model, "uses_common_shock")) {
    metrics <- full_parameter_metrics(
      fit, simulation, task$scenario, experiment, task$replicate
    )
  } else {
    metrics <- pilot_parameter_metrics(
      fit, simulation, task$scenario, task$replicate,
      include_spatial = FALSE
    )
    metrics$family <- experiment
    metrics <- metrics[, c(
      "scenario", "family", setdiff(names(metrics), c("scenario", "family"))
    )]
  }
  metrics$model <- task$model
  metrics$model_label <- task$model_label
  metrics$applicable <- TRUE
  metrics[, c(
    "model", "model_label", "applicable",
    setdiff(names(metrics), c("model", "model_label", "applicable"))
  )]
}

run_clustering_comparison_task <- function(
    task, results_dir, runs, burn, thin, n_side, n_times, n_neighbors,
    fixed_prior, reuse = TRUE, sampler_source = "sampler_global_beta_single_w.cpp"
) {
  validate_dcs_model_type(task$model)
  task_dir <- file.path(
    results_dir, task$model, task$scenario,
    sprintf("replicate_%03d", task$replicate)
  )
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  result_file <- file.path(task_dir, "metrics.rds")
  chain_seed <- comparison_chain_seed(task, 20700000L)
  settings <- list(
    version = 1L, model = task$model, sampler_source = sampler_source,
    scenario = task$scenario, geography = task$geography,
    separation = task$separation, replicate = task$replicate,
    n_side = n_side, n = n_side^2, n_times = n_times,
    n_neighbors = n_neighbors, runs = runs, burn = burn, thin = thin,
    data_seed_invariant_across_models = TRUE, chain_seed = chain_seed
  )
  if (reuse && file.exists(result_file)) {
    cached <- tryCatch(readRDS(result_file), error = function(error) NULL)
    if (!is.null(cached) && identical(cached$settings, settings)) {
      return(invisible(list(status = "reused", file = result_file)))
    }
  }

  simulation <- simulate_clustering_recovery_data(
    task, n_side = n_side, n_times = n_times, n_neighbors = n_neighbors
  )
  fit <- run_global_beta_validation_chain(
    simulation = simulation,
    runs = runs,
    burn = burn,
    thin = thin,
    beta_update_mode = "collapsed",
    sigma2_update_mode = "collapsed",
    tau2_update_mode = "pairwise_interwoven",
    spatial_update_mode = "collapsed",
    spatial_process_model = "nngp",
    model_type = task$model,
    prior_override = fixed_prior,
    seed = chain_seed,
    overdispersed_initialization = TRUE,
    save_latent = FALSE,
    save_w = FALSE,
    show_progress = FALSE
  )
  partition <- comparison_partition_metrics(fit, simulation, task)
  parameters <- comparison_parameter_metrics(
    fit, simulation, task, "model_comparison_clustering"
  )
  representative <- NULL
  posterior_similarity <- NULL
  if (dcs_model_uses(task$model, "uses_partition")) {
    posterior_similarity <- posterior_similarity_matrix(fit$cluster)
    representative <- dahl_partition(
      fit$cluster, posterior_similarity
    )$partition
  }
  result <- list(
    settings = settings,
    truth = simulation$true,
    partition = partition,
    parameters = parameters,
    representative_partition = representative,
    posterior_similarity = posterior_similarity,
    mcmc = list(
      elapsed_seconds = fit$elapsed_seconds,
      retained_draws = nrow(fit$cluster),
      k_mean = mean(fit$k_clust),
      k_ess = if (length(unique(fit$k_clust)) > 1L) {
        effective_sample_size(fit$k_clust)
      } else {
        NA_real_
      },
      phi_acceptance = fit$accphi,
      spatial_acceptance = fit$acc_spatial,
      structural_validation = list(
        model_type = fit$model_type,
        uses_common_shock = fit$uses_common_shock,
        uses_partition = fit$uses_partition,
        uses_spatial_partition = fit$uses_spatial_partition,
        all_k_equal_one = all(fit$k_clust == 1L),
        max_abs_u_mean = max(abs(fit$u_mean)),
        spatial_chains_all_na = all(is.na(fit$kappa)) &&
          all(is.na(fit$sigma2_u))
      )
    )
  )
  atomic_save_rds(result, result_file)
  invisible(list(status = "completed", file = result_file, result = result))
}

comparison_prediction_metrics <- function(fit, generated, task) {
  processes <- c("w", "y")
  if (dcs_model_uses(task$model, "uses_common_shock")) {
    processes <- c("u", processes)
  }
  rows <- do.call(rbind, lapply(processes, function(process) {
    prediction_process_metrics(
      fit$prediction[[process]],
      generated$prediction[[paste0("truth_", process)]],
      generated$prediction$location_type,
      process
    )
  }))
  rows$model <- task$model
  rows$model_label <- task$model_label
  rows$scenario <- task$scenario
  rows$geography <- task$geography
  rows$separation <- task$separation
  rows$replicate <- task$replicate
  rows$prediction_seconds <- fit$prediction$seconds
  rows$elapsed_seconds <- fit$elapsed_seconds
  rows[, c(
    "model", "model_label", "scenario", "geography", "separation",
    "replicate", setdiff(names(rows), c(
      "model", "model_label", "scenario", "geography", "separation",
      "replicate"
    ))
  )]
}

comparison_regime_prediction_metrics <- function(fit, generated, task) {
  if (!dcs_model_uses(task$model, "uses_partition")) return(NULL)
  metrics <- prediction_regime_metrics(fit, generated, task)
  metrics$model <- task$model
  metrics$model_label <- task$model_label
  metrics[, c("model", "model_label", setdiff(names(metrics), c("model", "model_label")))]
}

run_prediction_comparison_task <- function(
    task, results_dir, runs, burn, thin, prediction_thin, n_side, n_times,
    n_neighbors, fixed_prior, reuse = TRUE,
    sampler_source = "sampler_global_beta_single_w_prediction.cpp"
) {
  validate_dcs_model_type(task$model)
  task_dir <- file.path(
    results_dir, task$model, task$scenario,
    sprintf("replicate_%03d", task$replicate)
  )
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  result_file <- file.path(task_dir, "metrics.rds")
  chain_seed <- comparison_chain_seed(task, 20800000L)
  prediction_seed <- chain_seed + 500000L
  settings <- list(
    version = 1L, model = task$model, sampler_source = sampler_source,
    scenario = task$scenario, geography = task$geography,
    separation = task$separation, replicate = task$replicate,
    n_side = n_side, n = n_side^2, n_times = n_times,
    n_neighbors = n_neighbors, runs = runs, burn = burn, thin = thin,
    prediction_thin = prediction_thin,
    data_seed_invariant_across_models = TRUE,
    chain_seed = chain_seed, prediction_seed = prediction_seed
  )
  if (reuse && file.exists(result_file)) {
    cached <- tryCatch(readRDS(result_file), error = function(error) NULL)
    if (!is.null(cached) && identical(cached$settings, settings)) {
      return(invisible(list(status = "reused", file = result_file)))
    }
  }

  generated <- simulate_prediction_experiment_data(
    task, n_side, n_times, n_neighbors
  )
  prediction <- generated$prediction
  fit <- run_global_beta_validation_chain(
    simulation = generated$simulation,
    runs = runs,
    burn = burn,
    thin = thin,
    beta_update_mode = "collapsed",
    sigma2_update_mode = "collapsed",
    tau2_update_mode = "pairwise_interwoven",
    spatial_update_mode = "collapsed",
    spatial_process_model = "nngp",
    model_type = task$model,
    prior_override = fixed_prior,
    seed = chain_seed,
    overdispersed_initialization = TRUE,
    save_latent = FALSE,
    save_w = FALSE,
    show_progress = FALSE,
    sampler_extra_args = list(
      prediction_locations = prediction$locations,
      prediction_locations_similarity = prediction$locations_similarity,
      prediction_neighbors = prediction$neighbors,
      prediction_X = prediction$X,
      prediction_thin = as.integer(prediction_thin),
      prediction_seed = as.integer(prediction_seed),
      prediction_truth_u = prediction$truth_u,
      prediction_truth_w = prediction$truth_w,
      prediction_truth_y = prediction$truth_y,
      prediction_true_cluster = as.integer(prediction$true_cluster),
      prediction_cluster_phi = prediction$cluster_phi,
      prediction_cluster_tau2 = prediction$cluster_tau2
    )
  )
  result <- list(
    settings = settings,
    prediction = comparison_prediction_metrics(fit, generated, task),
    regime_prediction = comparison_regime_prediction_metrics(
      fit, generated, task
    ),
    partition = comparison_partition_metrics(
      fit, generated$simulation, task
    ),
    parameters = comparison_parameter_metrics(
      fit, generated$simulation, task, "model_comparison_prediction"
    ),
    runtime = list(
      elapsed_seconds = fit$elapsed_seconds,
      prediction_seconds = fit$prediction$seconds,
      retained_prediction_draws = fit$prediction$retained_draws,
      prediction_memory_mb = fit$prediction$working_memory_mb,
      structural_validation = list(
        model_type = fit$model_type,
        uses_common_shock = fit$uses_common_shock,
        uses_partition = fit$uses_partition,
        uses_spatial_partition = fit$uses_spatial_partition,
        all_k_equal_one = all(fit$k_clust == 1L),
        max_abs_u_mean = max(abs(fit$u_mean)),
        max_abs_predicted_u_mean = max(abs(fit$prediction$u$mean)),
        spatial_chains_all_na = all(is.na(fit$kappa)) &&
          all(is.na(fit$sigma2_u))
      )
    )
  )
  atomic_save_rds(result, result_file)
  invisible(list(status = "completed", file = result_file, result = result))
}

collect_model_comparison_results <- function(results_dir) {
  files <- list.files(
    results_dir, pattern = "^metrics[.]rds$", recursive = TRUE,
    full.names = TRUE
  )
  if (!length(files)) return(list())
  lapply(files, readRDS)
}

summarize_metric_long <- function(data, id, metric_names) {
  if (is.null(data) || !nrow(data)) return(data.frame())
  rows <- lapply(metric_names, function(metric) {
    if (!metric %in% names(data)) return(NULL)
    groups <- interaction(data[, id, drop = FALSE], drop = TRUE, lex.order = TRUE)
    split_values <- split(data[[metric]], groups)
    split_rows <- split(data[, id, drop = FALSE], groups)
    do.call(rbind, Map(function(values, keys) {
      values <- values[is.finite(values)]
      if (!length(values)) return(NULL)
      data.frame(
        keys[1L, , drop = FALSE], metric = metric,
        mean = mean(values), sd = stats::sd(values),
        q025 = unname(stats::quantile(values, 0.025)),
        q975 = unname(stats::quantile(values, 0.975)),
        n_replications = length(values), check.names = FALSE
      )
    }, split_values, split_rows))
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
