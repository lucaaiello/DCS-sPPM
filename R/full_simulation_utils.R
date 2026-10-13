full_simulation_scenarios <- function() {
  scenarios <- list()
  scenario_index <- 0L

  add_scenario <- function(id, family, replications, simulation_args) {
    scenario_index <<- scenario_index + 1L
    scenarios[[scenario_index]] <<- list(
      scenario_index = scenario_index,
      id = id,
      family = family,
      replications = as.integer(replications),
      simulation_args = simulation_args
    )
  }

  for (sigma2_u in c(0.15, 0.35, 0.75)) {
    for (kappa in c(0.5, 1.5, 3.0)) {
      add_scenario(
        sprintf(
          "main_sigma2u_%03d_kappa_%03d",
          round(100 * sigma2_u),
          round(100 * kappa)
        ),
        "primary_spatial_grid",
        30L,
        list(sigma2_u = sigma2_u, kappa = kappa)
      )
    }
  }

  add_scenario(
    "stress_weak_regimes",
    "stress_test",
    10L,
    list(
      cluster_phi = c(0.35, 0.55, 0.75),
      cluster_tau2 = c(0.25, 0.35, 0.50)
    )
  )
  add_scenario(
    "stress_high_measurement_noise",
    "stress_test",
    10L,
    list(sigma2 = 0.25)
  )
  add_scenario(
    "stress_exact_gp_dgp",
    "stress_test",
    10L,
    list(spatial_dgp = "exact_gp")
  )
  add_scenario(
    "stress_larger_spatial_grid",
    "stress_test",
    10L,
    list(n_side = 8L)
  )
  scenarios
}

scalability_simulation_scenarios <- function() {
  list(
    list(
      scenario_index = 14L,
      id = "scalability_n100_t120",
      family = "scalability",
      replications = 30L,
      simulation_args = list(n_side = 10L, n_times = 120L)
    ),
    list(
      scenario_index = 15L,
      id = "scalability_n036_t365",
      family = "scalability",
      replications = 30L,
      simulation_args = list(n_side = 6L, n_times = 365L)
    ),
    list(
      scenario_index = 16L,
      id = "scalability_n144_t365",
      family = "scalability",
      replications = 30L,
      simulation_args = list(n_side = 12L, n_times = 365L)
    )
  )
}

full_simulation_manifest <- function(scenarios = full_simulation_scenarios()) {
  rows <- lapply(seq_along(scenarios), function(index) {
    scenario <- scenarios[[index]]
    arguments <- modifyList(
      list(
        n_side = 6L,
        n_times = 120L,
        n_neighbors = 8L,
        sigma2 = 0.08,
        sigma2_u = 0.35,
        kappa = 1.5,
        cluster_phi = c(0.15, 0.55, 0.88),
        cluster_tau2 = c(0.12, 0.35, 0.75),
        spatial_dgp = "nngp"
      ),
      scenario$simulation_args
    )
    data.frame(
      scenario_index = scenario$scenario_index,
      scenario = scenario$id,
      family = scenario$family,
      replications = scenario$replications,
      n_side = arguments$n_side,
      n = arguments$n_side^2,
      n_times = arguments$n_times,
      n_neighbors = arguments$n_neighbors,
      sigma2 = arguments$sigma2,
      sigma2_u = arguments$sigma2_u,
      kappa = arguments$kappa,
      cluster_phi = paste(arguments$cluster_phi, collapse = ";"),
      cluster_tau2 = paste(arguments$cluster_tau2, collapse = ";"),
      spatial_dgp = arguments$spatial_dgp,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

full_simulation_tasks <- function(scenarios = full_simulation_scenarios()) {
  tasks <- list()
  task_index <- 0L
  for (replicate in seq_len(max(vapply(
    scenarios,
    function(scenario) scenario$replications,
    integer(1)
  )))) {
    for (scenario_index in seq_along(scenarios)) {
      scenario <- scenarios[[scenario_index]]
      if (replicate <= scenario$replications) {
        task_index <- task_index + 1L
        tasks[[task_index]] <- list(
          scenario_index = scenario$scenario_index,
          scenario = scenario,
          replicate = replicate
        )
      }
    }
  }
  tasks
}

full_simulation_worker_preflight <- function(fixed_prior, seed = 20270000L) {
  simulation <- simulate_dcs_sppm_global_beta(
    n_side = 3L,
    n_times = 8L,
    n_neighbors = 3L,
    seed = seed
  )
  fit <- run_global_beta_validation_chain(
    simulation = simulation,
    runs = 20L,
    burn = 10L,
    beta_update_mode = "collapsed",
    sigma2_update_mode = "collapsed",
    tau2_update_mode = "pairwise_interwoven",
    spatial_update_mode = "collapsed",
    prior_override = fixed_prior,
    seed = seed + 1L,
    overdispersed_initialization = TRUE,
    save_latent = FALSE,
    save_w = FALSE
  )
  if (!all(is.finite(fit$sigma2_u)) || !all(is.finite(fit$kappa))) {
    stop("Sparse-sampler worker preflight returned non-finite draws.")
  }
  TRUE
}

variation_of_information <- function(reference, estimate) {
  contingency <- table(reference, estimate)
  probabilities <- contingency / sum(contingency)
  reference_probabilities <- rowSums(probabilities)
  estimate_probabilities <- colSums(probabilities)
  entropy_reference <- -sum(
    reference_probabilities * log(reference_probabilities)
  )
  entropy_estimate <- -sum(
    estimate_probabilities * log(estimate_probabilities)
  )
  positive <- probabilities > 0
  expected <- outer(reference_probabilities, estimate_probabilities)
  mutual_information <- sum(
    probabilities[positive] * log(probabilities[positive] / expected[positive])
  )
  entropy_reference + entropy_estimate - 2 * mutual_information
}

pairwise_partition_accuracy <- function(reference, estimate) {
  true_pairs <- outer(reference, reference, FUN = "==")
  estimated_pairs <- outer(estimate, estimate, FUN = "==")
  selected <- upper.tri(true_pairs)
  true_pairs <- true_pairs[selected]
  estimated_pairs <- estimated_pairs[selected]
  sensitivity <- if (any(true_pairs)) {
    mean(estimated_pairs[true_pairs])
  } else {
    NA_real_
  }
  specificity <- if (any(!true_pairs)) {
    mean(!estimated_pairs[!true_pairs])
  } else {
    NA_real_
  }
  c(sensitivity = sensitivity, specificity = specificity)
}

full_parameter_metrics <- function(
    fit,
    simulation,
    scenario,
    family,
    replicate
) {
  metrics <- pilot_parameter_metrics(
    fit,
    simulation,
    scenario,
    replicate
  )
  microergodic_draws <- fit$sigma2_u / fit$kappa
  microergodic_truth <-
    simulation$true$sigma2_u / simulation$true$kappa
  interval <- stats::quantile(
    microergodic_draws,
    c(0.025, 0.975),
    names = FALSE
  )
  microergodic <- data.frame(
    scenario = scenario,
    replicate = replicate,
    parameter_group = "spatial",
    parameter = "sigma2_u_over_kappa",
    truth = microergodic_truth,
    posterior_mean = mean(microergodic_draws),
    posterior_sd = stats::sd(microergodic_draws),
    q025 = interval[1],
    q975 = interval[2],
    covered = microergodic_truth >= interval[1] &&
      microergodic_truth <= interval[2],
    interval_width = diff(interval),
    ess = effective_sample_size(microergodic_draws),
    lag1_autocorrelation = stats::acf(
      microergodic_draws,
      lag.max = 1L,
      plot = FALSE
    )$acf[2],
    stringsAsFactors = FALSE
  )
  metrics <- rbind(metrics, microergodic)
  metrics$family <- family
  metrics[, c("scenario", "family", setdiff(names(metrics), c("scenario", "family")))]
}

full_partition_metrics <- function(
    fit,
    simulation,
    scenario,
    family,
    replicate
) {
  similarity <- posterior_similarity_matrix(fit$cluster)
  representative <- dahl_partition(fit$cluster, similarity)$partition
  ari <- adjusted_rand_index(simulation$true$cluster, representative)
  pairwise <- pairwise_partition_accuracy(
    simulation$true$cluster,
    representative
  )
  data.frame(
    scenario = scenario,
    family = family,
    replicate = replicate,
    adjusted_rand_index = ari,
    variation_of_information = variation_of_information(
      simulation$true$cluster,
      representative
    ),
    pairwise_sensitivity = unname(pairwise["sensitivity"]),
    pairwise_specificity = unname(pairwise["specificity"]),
    exact_partition_recovery = isTRUE(all.equal(ari, 1)),
    true_clusters = length(unique(simulation$true$cluster)),
    posterior_mean_clusters = mean(fit$k_clust),
    posterior_mode_clusters = as.integer(names(which.max(table(fit$k_clust)))),
    spatial_acceptance = fit$acc_spatial,
    elapsed_seconds = fit$elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

summarize_full_parameters <- function(metrics) {
  summary <- summarize_pilot_parameters(metrics)
  family_lookup <- unique(metrics[, c("scenario", "family")])
  summary <- merge(summary, family_lookup, by = "scenario", sort = FALSE)
  summary$coverage_mcse <- sqrt(
    summary$coverage * (1 - summary$coverage) / summary$replications
  )
  summary[, c("scenario", "family", setdiff(names(summary), c("scenario", "family")))]
}

summarize_full_partitions <- function(metrics) {
  groups <- split(metrics, metrics$scenario)
  summary <- do.call(rbind, lapply(groups, function(group) {
    data.frame(
      scenario = group$scenario[1],
      family = group$family[1],
      replications = nrow(group),
      mean_ari = mean(group$adjusted_rand_index),
      median_ari = stats::median(group$adjusted_rand_index),
      mean_variation_of_information = mean(group$variation_of_information),
      mean_pairwise_sensitivity = mean(group$pairwise_sensitivity),
      mean_pairwise_specificity = mean(group$pairwise_specificity),
      exact_recovery_rate = mean(group$exact_partition_recovery),
      mean_posterior_clusters = mean(group$posterior_mean_clusters),
      correct_modal_cluster_rate = mean(
        group$posterior_mode_clusters == group$true_clusters
      ),
      mean_spatial_acceptance = mean(group$spatial_acceptance),
      mean_elapsed_seconds = mean(group$elapsed_seconds),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  summary[order(summary$scenario), ]
}

write_full_simulation_summaries <- function(results_dir) {
  metric_files <- list.files(
    results_dir,
    pattern = "^metrics\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(metric_files)) {
    stop("No completed full-simulation checkpoints were found.")
  }
  results <- lapply(metric_files, readRDS)
  parameter_metrics <- do.call(rbind, lapply(results, `[[`, "parameters"))
  partition_metrics <- do.call(rbind, lapply(results, `[[`, "partition"))
  parameter_summary <- summarize_full_parameters(parameter_metrics)
  partition_summary <- summarize_full_partitions(partition_metrics)

  key_parameters <- c("sigma2", "sigma2_u", "kappa")
  key_metrics <- parameter_metrics[
    parameter_metrics$parameter %in% key_parameters,
  ]
  ess_flags <- key_metrics[key_metrics$ess < 100, ]
  failure_files <- list.files(
    results_dir,
    pattern = "^failure\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  failures <- if (length(failure_files)) {
    do.call(rbind, lapply(failure_files, readRDS))
  } else {
    data.frame(
      scenario = character(),
      family = character(),
      replicate = integer(),
      message = character()
    )
  }

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
    ess_flags,
    file.path(results_dir, "ess_flags.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    failures,
    file.path(results_dir, "failures.csv"),
    row.names = FALSE
  )

  list(
    parameter_metrics = parameter_metrics,
    parameter_summary = parameter_summary,
    partition_metrics = partition_metrics,
    partition_summary = partition_summary,
    ess_flags = ess_flags,
    failures = failures
  )
}
