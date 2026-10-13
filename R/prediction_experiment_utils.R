prediction_evaluation_locations <- function() {
  interior_axis <- seq(0.15, 4.85, length.out = 10L)
  interior <- as.matrix(expand.grid(
    x = interior_axis,
    y = interior_axis,
    KEEP.OUT.ATTRS = FALSE
  ))
  edge_axis <- seq(0.15, 4.85, length.out = 10L)
  exterior <- rbind(
    cbind(x = edge_axis, y = -0.35),
    cbind(x = edge_axis, y = 5.35),
    cbind(x = -0.35, y = edge_axis),
    cbind(x = 5.35, y = edge_axis)
  )
  locations <- rbind(interior, exterior)
  colnames(locations) <- c("x", "y")
  list(
    locations = locations,
    type = c(rep("interpolation", nrow(interior)), rep("extrapolation", nrow(exterior)))
  )
}

prediction_regions <- function(
    geography,
    replicate,
    training_locations,
    prediction_locations,
    seed_base = 20310000L
) {
  training_truth <- clustering_recovery_truth(
    geography,
    replicate,
    training_locations,
    seed_base
  )
  if (geography == "random_allocation") {
    set.seed(seed_base + 30000L + as.integer(replicate))
    prediction_region <- sample(
      rep(1:3, length.out = nrow(prediction_locations)),
      replace = FALSE
    )
  } else {
    training_score <- if (geography == "regular_strips") {
      training_locations[, 1]
    } else {
      training_locations[, 1] -
        0.85 * sin(2 * pi * training_locations[, 2] / 5) -
        0.20 * sin(4 * pi * training_locations[, 2] / 5)
    }
    prediction_score <- if (geography == "regular_strips") {
      prediction_locations[, 1]
    } else {
      prediction_locations[, 1] -
        0.85 * sin(2 * pi * prediction_locations[, 2] / 5) -
        0.20 * sin(4 * pi * prediction_locations[, 2] / 5)
    }
    ordered <- sort(training_score)
    region_size <- length(ordered) / 3L
    breaks <- c(
      mean(ordered[c(region_size, region_size + 1L)]),
      mean(ordered[c(2L * region_size, 2L * region_size + 1L)])
    )
    prediction_region <- as.integer(cut(
      prediction_score,
      breaks = c(-Inf, breaks, Inf),
      labels = FALSE,
      include.lowest = TRUE
    ))
  }

  list(
    training_cluster = training_truth$true_cluster,
    prediction_cluster = training_truth$region_to_parameter[prediction_region],
    training_region = training_truth$region,
    prediction_region = as.integer(prediction_region),
    region_to_parameter = training_truth$region_to_parameter
  )
}

nearest_observed_neighbors <- function(
    observed_locations,
    prediction_locations,
    n_neighbors = 10L
) {
  observed_locations <- as.matrix(observed_locations)
  prediction_locations <- as.matrix(prediction_locations)
  n_neighbors <- min(as.integer(n_neighbors), nrow(observed_locations))
  neighbors <- matrix(0L, nrow(prediction_locations), n_neighbors)
  for (prediction in seq_len(nrow(prediction_locations))) {
    distances <- sqrt(rowSums(
      (observed_locations - matrix(
        prediction_locations[prediction, ],
        nrow(observed_locations),
        2L,
        byrow = TRUE
      ))^2
    ))
    neighbors[prediction, ] <- order(distances)[seq_len(n_neighbors)]
  }
  neighbors
}

split_time_major_design <- function(X, n_total, selected_sites, n_times) {
  rows <- unlist(lapply(seq_len(n_times), function(time) {
    (time - 1L) * n_total + selected_sites
  }), use.names = FALSE)
  X[rows, , drop = FALSE]
}

simulate_prediction_experiment_data <- function(
    task,
    n_side = 12L,
    n_times = 100L,
    n_neighbors = 10L,
    seed_base = 20500000L
) {
  training_locations <- clustering_recovery_locations(n_side)
  evaluation <- prediction_evaluation_locations()
  prediction_locations <- evaluation$locations
  geography_truth <- prediction_regions(
    task$geography,
    task$replicate,
    training_locations,
    prediction_locations
  )
  representatives <- clustering_recovery_separations()[[task$separation]]
  full_locations <- rbind(training_locations, prediction_locations)
  full_cluster <- c(
    geography_truth$training_cluster,
    geography_truth$prediction_cluster
  )
  geography_index <- match(
    task$geography,
    clustering_recovery_geographies()$geography
  )
  data_seed <- seed_base + 100000L * geography_index + task$replicate
  full <- simulate_dcs_sppm_global_beta(
    n_side = n_side,
    n_times = n_times,
    n_neighbors = n_neighbors,
    cluster_phi = representatives$phi,
    cluster_tau2 = representatives$tau2,
    spatial_dgp = "nngp",
    seed = data_seed,
    locations = full_locations,
    true_cluster = full_cluster
  )

  n_training <- nrow(training_locations)
  n_prediction <- nrow(prediction_locations)
  n_total <- n_training + n_prediction
  training_sites <- seq_len(n_training)
  prediction_sites <- n_training + seq_len(n_prediction)
  training_design <- split_time_major_design(
    full$X,
    n_total,
    training_sites,
    n_times
  )
  prediction_design <- split_time_major_design(
    full$X,
    n_total,
    prediction_sites,
    n_times
  )

  simulation <- list(
    Y = full$Y[, training_sites, drop = FALSE],
    X = training_design,
    locations = training_locations,
    NN_ind = full$NN_ind[training_sites, , drop = FALSE],
    true = list(
      beta = full$true$beta,
      sigma2 = full$true$sigma2,
      sigma2_u = full$true$sigma2_u,
      kappa = full$true$kappa,
      nu = full$true$nu,
      spatial_dgp = full$true$spatial_dgp,
      cluster = geography_truth$training_cluster,
      phi = representatives$phi[geography_truth$training_cluster],
      tau2 = representatives$tau2[geography_truth$training_cluster],
      cluster_phi = representatives$phi,
      cluster_tau2 = representatives$tau2,
      w = full$true$w[, training_sites, drop = FALSE],
      u = full$true$u[, training_sites, drop = FALSE],
      eta = full$true$eta[, training_sites, drop = FALSE],
      fixed_effect = full$true$fixed_effect[, training_sites, drop = FALSE]
    ),
    dimensions = c(n = n_training, T = n_times, p = ncol(training_design)),
    seed = data_seed,
    design_order = full$design_order
  )
  location_center <- colMeans(training_locations)
  location_scale <- apply(training_locations, 2L, stats::sd)
  prediction_similarity <- sweep(
    sweep(prediction_locations, 2L, location_center, "-"),
    2L,
    location_scale,
    "/"
  )

  list(
    simulation = simulation,
    prediction = list(
      locations = prediction_locations,
      location_type = evaluation$type,
      locations_similarity = prediction_similarity,
      neighbors = nearest_observed_neighbors(
        training_locations,
        prediction_locations,
        n_neighbors
      ),
      X = prediction_design,
      truth_u = full$true$u[, prediction_sites, drop = FALSE],
      truth_w = full$true$w[, prediction_sites, drop = FALSE],
      truth_y = full$Y[, prediction_sites, drop = FALSE],
      true_cluster = geography_truth$prediction_cluster,
      cluster_phi = representatives$phi,
      cluster_tau2 = representatives$tau2
    )
  )
}

prediction_process_metrics <- function(summary, truth, location_type, process) {
  strata <- c("all", "interpolation", "extrapolation")
  do.call(rbind, lapply(strata, function(stratum) {
    selected <- if (stratum == "all") {
      seq_along(location_type)
    } else {
      which(location_type == stratum)
    }
    error <- summary$mean[, selected, drop = FALSE] -
      truth[, selected, drop = FALSE]
    data.frame(
      process = process,
      stratum = stratum,
      n_locations = length(selected),
      n_times = nrow(truth),
      bias = mean(error),
      rmse = sqrt(mean(summary$squared_error[, selected, drop = FALSE])),
      mae = mean(summary$absolute_error[, selected, drop = FALSE]),
      coverage = mean(summary$coverage[, selected, drop = FALSE]),
      interval_width = mean(summary$interval_width[, selected, drop = FALSE]),
      crps = mean(summary$crps[, selected, drop = FALSE]),
      mean_posterior_sd = mean(summary$sd[, selected, drop = FALSE]),
      stringsAsFactors = FALSE
    )
  }))
}

prediction_experiment_metrics <- function(fit, generated, task) {
  prediction <- generated$prediction
  rows <- rbind(
    prediction_process_metrics(
      fit$prediction$u,
      prediction$truth_u,
      prediction$location_type,
      "u"
    ),
    prediction_process_metrics(
      fit$prediction$w,
      prediction$truth_w,
      prediction$location_type,
      "w"
    ),
    prediction_process_metrics(
      fit$prediction$y,
      prediction$truth_y,
      prediction$location_type,
      "y"
    )
  )
  rows$scenario <- task$scenario
  rows$geography <- task$geography
  rows$separation <- task$separation
  rows$replicate <- task$replicate
  rows$prediction_seconds <- fit$prediction$seconds
  rows$elapsed_seconds <- fit$elapsed_seconds
  rows$prediction_draws <- fit$prediction$retained_draws
  rows$prediction_memory_mb <- fit$prediction$working_memory_mb
  rows[, c(
    "scenario", "geography", "separation", "replicate",
    setdiff(names(rows), c("scenario", "geography", "separation", "replicate"))
  )]
}

prediction_regime_metrics <- function(fit, generated, task) {
  type <- generated$prediction$location_type
  probability <- fit$prediction$regime$true_class_probability
  modal <- fit$prediction$regime$modal_class
  truth <- generated$prediction$true_cluster
  strata <- c("all", "interpolation", "extrapolation")
  do.call(rbind, lapply(strata, function(stratum) {
    selected <- if (stratum == "all") seq_along(type) else which(type == stratum)
    data.frame(
      scenario = task$scenario,
      geography = task$geography,
      separation = task$separation,
      replicate = task$replicate,
      stratum = stratum,
      classification_accuracy = mean(modal[selected] == truth[selected]),
      mean_true_class_probability = mean(probability[selected]),
      stringsAsFactors = FALSE
    )
  }))
}

run_prediction_experiment_task <- function(
    task,
    results_dir,
    runs,
    burn,
    thin,
    prediction_thin,
    n_side,
    n_times,
    n_neighbors,
    fixed_prior,
    deadline,
    reuse,
    sampler_source,
    chain_progress_file = ""
) {
  task_dir <- file.path(
    results_dir,
    task$scenario,
    sprintf("replicate_%03d", task$replicate)
  )
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  metrics_file <- file.path(task_dir, "metrics.rds")
  failure_file <- file.path(task_dir, "failure.rds")
  chain_seed <- 20600000L + 10000L * task$scenario_index + task$replicate
  prediction_seed <- chain_seed + 500000L
  expected_settings <- list(
    version = 1L,
    sampler_source = sampler_source,
    scenario = task$scenario,
    geography = task$geography,
    separation = task$separation,
    replicate = task$replicate,
    n_side = n_side,
    n = n_side^2,
    n_times = n_times,
    n_neighbors = n_neighbors,
    runs = runs,
    burn = burn,
    thin = thin,
    prediction_thin = prediction_thin,
    chain_seed = chain_seed,
    prediction_seed = prediction_seed
  )

  if (reuse && file.exists(metrics_file)) {
    cached <- tryCatch(readRDS(metrics_file), error = function(error) NULL)
    if (!is.null(cached) && identical(cached$settings, expected_settings)) {
      return(data.frame(
        scenario = task$scenario,
        replicate = task$replicate,
        status = "reused",
        stringsAsFactors = FALSE
      ))
    }
  }
  if (file.exists(metrics_file)) unlink(metrics_file)
  if (as.numeric(Sys.time()) >= deadline) {
    return(data.frame(
      scenario = task$scenario,
      replicate = task$replicate,
      status = "deadline",
      stringsAsFactors = FALSE
    ))
  }

  tryCatch({
    generated <- simulate_prediction_experiment_data(
      task,
      n_side,
      n_times,
      n_neighbors
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
      prior_override = fixed_prior,
      seed = chain_seed,
      overdispersed_initialization = TRUE,
      save_latent = FALSE,
      save_w = FALSE,
      progress_file = chain_progress_file,
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
    prediction_metrics <- prediction_experiment_metrics(fit, generated, task)
    regime_metrics <- prediction_regime_metrics(fit, generated, task)
    partition <- full_partition_metrics(
      fit,
      generated$simulation,
      task$scenario,
      "prediction_experiment",
      task$replicate
    )
    parameters <- full_parameter_metrics(
      fit,
      generated$simulation,
      task$scenario,
      "prediction_experiment",
      task$replicate
    )
    compact <- list(
      settings = expected_settings,
      prediction = prediction_metrics,
      regime_prediction = regime_metrics,
      partition = partition,
      parameters = parameters,
      runtime = list(
        elapsed_seconds = fit$elapsed_seconds,
        prediction_seconds = fit$prediction$seconds,
        prediction_draws = fit$prediction$retained_draws,
        prediction_memory_mb = fit$prediction$working_memory_mb,
        openmp_threads = fit$openmp_threads
      )
    )
    atomic_save_rds(compact, metrics_file)
    if (file.exists(failure_file)) unlink(failure_file)
    rm(fit, generated)
    invisible(gc(FALSE))
    data.frame(
      scenario = task$scenario,
      replicate = task$replicate,
      status = "completed",
      stringsAsFactors = FALSE
    )
  }, error = function(error) {
    atomic_save_rds(
      list(settings = expected_settings, message = conditionMessage(error)),
      failure_file
    )
    data.frame(
      scenario = task$scenario,
      replicate = task$replicate,
      status = "failed",
      stringsAsFactors = FALSE
    )
  })
}

run_prediction_experiment_worker <- function(
    tasks,
    worker,
    status_file,
    chain_progress_file,
    results_dir,
    runs,
    burn,
    thin,
    prediction_thin,
    n_side,
    n_times,
    n_neighbors,
    fixed_prior,
    deadline,
    reuse,
    sampler_source
) {
  total <- length(tasks)
  completed <- 0L
  failures <- 0L
  statuses <- vector("list", total)
  for (index in seq_along(tasks)) {
    task <- tasks[[index]]
    write_clustering_worker_status(
      status_file,
      worker,
      total,
      completed,
      "running",
      task$display_name,
      failures
    )
    writeLines("0", chain_progress_file)
    statuses[[index]] <- run_prediction_experiment_task(
      task = task,
      results_dir = results_dir,
      runs = runs,
      burn = burn,
      thin = thin,
      prediction_thin = prediction_thin,
      n_side = n_side,
      n_times = n_times,
      n_neighbors = n_neighbors,
      fixed_prior = fixed_prior,
      deadline = deadline,
      reuse = reuse,
      sampler_source = sampler_source,
      chain_progress_file = chain_progress_file
    )
    completed <- completed + 1L
    failures <- failures + as.integer(statuses[[index]]$status == "failed")
    writeLines("0", chain_progress_file)
    write_clustering_worker_status(
      status_file,
      worker,
      total,
      completed,
      if (completed == total) "done" else "running",
      task$display_name,
      failures
    )
  }
  do.call(rbind, statuses)
}
