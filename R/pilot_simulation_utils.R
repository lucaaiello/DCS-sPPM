pilot_fixed_prior <- function() {
  list(
    a_sigma = 3,
    b_sigma = 0.16,
    a_xi = 2,
    b_xi = 2,
    a_tau = 2,
    b_tau = 0.35,
    a_phi = 1,
    b_phi = 1,
    a_kappa = 3,
    b_kappa = 2,
    a_sigma2_u = 3,
    b_sigma2_u = 0.70
  )
}

pilot_scenarios <- function() {
  data.frame(
    scenario = c("low_spatial", "baseline", "high_spatial"),
    sigma2_u = c(0.15, 0.35, 0.75),
    stringsAsFactors = FALSE
  )
}

atomic_save_rds <- function(object, file) {
  temporary_file <- paste0(file, ".tmp")
  saveRDS(object, temporary_file)
  if (file.exists(file)) {
    unlink(file)
  }
  if (!file.rename(temporary_file, file)) {
    stop("Could not finalize `", file, "`.")
  }
  invisible(file)
}

pilot_parameter_metrics <- function(
    fit, simulation, scenario, replicate, include_spatial = TRUE
) {
  summarize_draws <- function(draws, truth, parameter, parameter_group) {
    interval <- stats::quantile(draws, c(0.025, 0.975), names = FALSE)
    data.frame(
      scenario = scenario,
      replicate = replicate,
      parameter_group = parameter_group,
      parameter = parameter,
      truth = truth,
      posterior_mean = mean(draws),
      posterior_sd = stats::sd(draws),
      q025 = interval[1],
      q975 = interval[2],
      covered = truth >= interval[1] && truth <= interval[2],
      interval_width = diff(interval),
      ess = effective_sample_size(draws),
      lag1_autocorrelation = stats::acf(
        draws,
        lag.max = 1L,
        plot = FALSE
      )$acf[2],
      stringsAsFactors = FALSE
    )
  }

  rows <- lapply(seq_along(simulation$true$beta), function(index) {
    summarize_draws(
      fit$beta[, index],
      simulation$true$beta[index],
      paste0("beta_", names(simulation$true$beta)[index]),
      "global"
    )
  })
  rows[[length(rows) + 1L]] <- summarize_draws(
    fit$sigma2[, 1], simulation$true$sigma2, "sigma2", "global"
  )
  if (include_spatial) {
    rows[[length(rows) + 1L]] <- summarize_draws(
      fit$sigma2_u, simulation$true$sigma2_u, "sigma2_u", "spatial"
    )
    rows[[length(rows) + 1L]] <- summarize_draws(
      fit$kappa, simulation$true$kappa, "kappa", "spatial"
    )
  }

  for (regime in sort(unique(simulation$true$cluster))) {
    sites <- which(simulation$true$cluster == regime)
    rows[[length(rows) + 1L]] <- summarize_draws(
      rowMeans(fit$phi[, sites, drop = FALSE]),
      simulation$true$cluster_phi[regime],
      paste0("phi_regime_", regime),
      "regime"
    )
    rows[[length(rows) + 1L]] <- summarize_draws(
      rowMeans(fit$tau2[, sites, drop = FALSE]),
      simulation$true$cluster_tau2[regime],
      paste0("tau2_regime_", regime),
      "regime"
    )
  }
  do.call(rbind, rows)
}

pilot_partition_metrics <- function(fit, simulation, scenario, replicate) {
  similarity <- posterior_similarity_matrix(fit$cluster)
  representative <- dahl_partition(fit$cluster, similarity)$partition
  ari <- adjusted_rand_index(simulation$true$cluster, representative)
  data.frame(
    scenario = scenario,
    replicate = replicate,
    adjusted_rand_index = ari,
    exact_partition_recovery = isTRUE(all.equal(ari, 1)),
    true_clusters = length(unique(simulation$true$cluster)),
    posterior_mean_clusters = mean(fit$k_clust),
    posterior_mode_clusters = as.integer(names(which.max(table(fit$k_clust)))),
    spatial_acceptance = fit$acc_spatial,
    elapsed_seconds = fit$elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

summarize_pilot_parameters <- function(metrics) {
  groups <- split(metrics, list(metrics$scenario, metrics$parameter), drop = TRUE)
  rows <- lapply(groups, function(group) {
    errors <- group$posterior_mean - group$truth
    data.frame(
      scenario = group$scenario[1],
      parameter_group = group$parameter_group[1],
      parameter = group$parameter[1],
      truth = group$truth[1],
      replications = nrow(group),
      mean_estimate = mean(group$posterior_mean),
      bias = mean(errors),
      relative_bias = mean(errors) / abs(group$truth[1]),
      rmse = sqrt(mean(errors^2)),
      empirical_sd = stats::sd(group$posterior_mean),
      mean_posterior_sd = mean(group$posterior_sd),
      coverage = mean(group$covered),
      mean_interval_width = mean(group$interval_width),
      median_ess = stats::median(group$ess),
      mean_lag1_autocorrelation = mean(group$lag1_autocorrelation),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  summary[order(summary$scenario, summary$parameter_group, summary$parameter), ]
}

summarize_pilot_partitions <- function(metrics) {
  groups <- split(metrics, metrics$scenario)
  rows <- lapply(groups, function(group) {
    data.frame(
      scenario = group$scenario[1],
      replications = nrow(group),
      mean_ari = mean(group$adjusted_rand_index),
      median_ari = stats::median(group$adjusted_rand_index),
      exact_recovery_rate = mean(group$exact_partition_recovery),
      mean_posterior_clusters = mean(group$posterior_mean_clusters),
      correct_modal_cluster_rate = mean(
        group$posterior_mode_clusters == group$true_clusters
      ),
      mean_spatial_acceptance = mean(group$spatial_acceptance),
      mean_elapsed_seconds = mean(group$elapsed_seconds),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  summary[order(summary$scenario), ]
}
