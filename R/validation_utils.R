combination_two <- function(x) x * (x - 1) / 2

adjusted_rand_index <- function(reference, estimate) {
  contingency <- table(reference, estimate)
  n <- sum(contingency)
  if (n < 2L) {
    return(NA_real_)
  }

  sum_cells <- sum(combination_two(contingency))
  sum_rows <- sum(combination_two(rowSums(contingency)))
  sum_columns <- sum(combination_two(colSums(contingency)))
  expected <- sum_rows * sum_columns / combination_two(n)
  maximum <- 0.5 * (sum_rows + sum_columns)
  denominator <- maximum - expected
  if (abs(denominator) < .Machine$double.eps) {
    return(as.numeric(identical(as.integer(reference), as.integer(estimate))))
  }
  (sum_cells - expected) / denominator
}

posterior_similarity_matrix <- function(cluster_draws) {
  cluster_draws <- as.matrix(cluster_draws)
  n <- ncol(cluster_draws)
  similarity <- matrix(0, nrow = n, ncol = n)
  for (draw in seq_len(nrow(cluster_draws))) {
    labels <- cluster_draws[draw, ]
    similarity <- similarity + outer(labels, labels, FUN = "==")
  }
  similarity / nrow(cluster_draws)
}

dahl_partition <- function(cluster_draws, similarity = NULL) {
  cluster_draws <- as.matrix(cluster_draws)
  if (is.null(similarity)) {
    similarity <- posterior_similarity_matrix(cluster_draws)
  }

  losses <- vapply(seq_len(nrow(cluster_draws)), function(draw) {
    allocation <- cluster_draws[draw, ]
    association <- outer(allocation, allocation, FUN = "==")
    sum((association - similarity)^2)
  }, numeric(1))

  list(
    partition = as.integer(cluster_draws[which.min(losses), ]),
    draw = which.min(losses),
    loss = min(losses),
    losses = losses
  )
}

align_partition_to_truth <- function(estimate, truth) {
  estimate <- as.integer(factor(estimate))
  truth <- as.integer(factor(truth))
  estimated_levels <- sort(unique(estimate))
  aligned <- integer(length(estimate))
  mapping <- integer(length(estimated_levels))

  for (index in seq_along(estimated_levels)) {
    members <- estimate == estimated_levels[index]
    counts <- table(factor(truth[members], levels = sort(unique(truth))))
    mapping[index] <- as.integer(names(counts)[which.max(counts)])
    aligned[members] <- mapping[index]
  }

  list(aligned = aligned, mapping = mapping)
}

effective_sample_size <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 4L || stats::var(x) == 0) {
    return(length(x))
  }
  if (requireNamespace("coda", quietly = TRUE)) {
    return(as.numeric(coda::effectiveSize(x)))
  }

  autocorrelation <- stats::acf(x, plot = FALSE, lag.max = min(1000L, length(x) - 1L))$acf[-1]
  positive <- autocorrelation[autocorrelation > 0]
  length(x) / (1 + 2 * sum(positive))
}

matern_correlation <- function(distance, kappa, nu) {
  correlation <- numeric(length(distance))
  zero <- distance == 0
  correlation[zero] <- 1
  if (all(zero)) {
    return(correlation)
  }

  positive_distance <- distance[!zero]
  if (abs(nu - 0.5) < 1e-12) {
    correlation[!zero] <- exp(-positive_distance / kappa)
  } else {
    argument <- sqrt(2 * nu) * positive_distance / kappa
    correlation[!zero] <-
      2^(1 - nu) / gamma(nu) * argument^nu * besselK(argument, nu)
  }
  correlation
}

profile_true_common_shock_kappa <- function(simulation, fit, grid = NULL) {
  truth <- simulation$true
  locations <- as.matrix(simulation$locations)
  neighbor_index <- as.matrix(simulation$NN_ind)
  common_shock <- as.matrix(truth$u)
  n <- ncol(common_shock)
  n_times <- nrow(common_shock)

  if (is.null(grid)) {
    grid <- seq(0.4 * truth$kappa, 2.4 * truth$kappa, length.out = 241L)
  }

  evaluate <- function(kappa) {
    log_conditional_variance <- 0
    quadratic_form <- 0

    for (site in seq_len(n)) {
      neighbors <- neighbor_index[site, ]
      neighbors <- neighbors[neighbors > 0L]

      if (!length(neighbors)) {
        coefficients <- numeric(0)
        conditional_variance <- 1
      } else {
        neighbor_locations <- locations[neighbors, , drop = FALSE]
        neighbor_distances <- as.matrix(stats::dist(neighbor_locations))
        neighbor_correlation <- matrix(
          matern_correlation(as.vector(neighbor_distances), kappa, truth$nu),
          nrow = length(neighbors),
          ncol = length(neighbors)
        )
        cross_distances <- sqrt(rowSums(
          (neighbor_locations - matrix(
            locations[site, ], nrow = length(neighbors), ncol = 2L, byrow = TRUE
          ))^2
        ))
        cross_correlation <- matern_correlation(cross_distances, kappa, truth$nu)
        coefficients <- solve(neighbor_correlation, cross_correlation)
        conditional_variance <- 1 - sum(cross_correlation * coefficients)
      }

      innovations <- common_shock[2:n_times, site]
      if (length(neighbors)) {
        innovations <- innovations - as.vector(
          common_shock[2:n_times, neighbors, drop = FALSE] %*% coefficients
        )
      }
      log_conditional_variance <-
        log_conditional_variance + log(conditional_variance)
      quadratic_form <-
        quadratic_form + sum(innovations^2) / conditional_variance
    }

    n_innovations <- n * (n_times - 1L)
    log_kappa_prior <-
      (fit$prior$a_kappa - 1) * log(kappa) - fit$prior$b_kappa * kappa
    conditional_log_posterior <-
      log_kappa_prior -
      0.5 * (n_times - 1L) * log_conditional_variance -
      0.5 * n_innovations * log(truth$sigma2_u) -
      0.5 * quadratic_form / truth$sigma2_u
    integrated_log_posterior <-
      log_kappa_prior -
      0.5 * (n_times - 1L) * log_conditional_variance -
      (fit$prior$a_sigma2_u + n_innovations / 2) *
        log(fit$prior$b_sigma2_u + quadratic_form / 2)
    posterior_sigma2_u_mean <-
      (fit$prior$b_sigma2_u + quadratic_form / 2) /
      (fit$prior$a_sigma2_u + n_innovations / 2 - 1)

    c(
      conditional_log_posterior = conditional_log_posterior,
      integrated_log_posterior = integrated_log_posterior,
      posterior_sigma2_u_mean = posterior_sigma2_u_mean,
      quadratic_form = quadratic_form
    )
  }

  profile <- as.data.frame(t(vapply(grid, evaluate, numeric(4))))
  profile$kappa <- grid
  profile$conditional_relative <-
    profile$conditional_log_posterior - max(profile$conditional_log_posterior)
  profile$integrated_relative <-
    profile$integrated_log_posterior - max(profile$integrated_log_posterior)
  profile
}

plot_kappa_profile <- function(profile, truth, posterior_draws, file) {
  grDevices::pdf(file, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(
    profile$kappa,
    profile$integrated_relative,
    type = "l",
    lwd = 2,
    col = "#0072B2",
    xlab = expression(kappa),
    ylab = "Relative log posterior",
    ylim = c(-50, 0),
    main = "Oracle profile using the true common shocks"
  )
  graphics::lines(
    profile$kappa,
    profile$conditional_relative,
    lwd = 2,
    lty = 2,
    col = "#D55E00"
  )
  graphics::abline(v = truth, col = "black", lty = 3)
  graphics::abline(v = stats::median(posterior_draws), col = "#009E73", lty = 3)
  graphics::legend(
    "bottomleft",
    legend = c(
      "sigma2_u integrated out",
      "conditional on true sigma2_u",
      "true kappa",
      "MCMC posterior median"
    ),
    col = c("#0072B2", "#D55E00", "black", "#009E73"),
    lty = c(1, 2, 3, 3),
    lwd = c(2, 2, 1, 1),
    bty = "n"
  )
}

summarize_scalar_draws <- function(draws, truth, parameter) {
  interval <- stats::quantile(draws, c(0.025, 0.5, 0.975), names = FALSE)
  data.frame(
    parameter = parameter,
    truth = truth,
    mean = mean(draws),
    median = interval[2],
    q025 = interval[1],
    q975 = interval[3],
    covered = truth >= interval[1] && truth <= interval[3],
    ess = effective_sample_size(draws),
    row.names = NULL
  )
}

global_parameter_summary <- function(fit, simulation) {
  beta_draws <- as.matrix(fit$beta)
  beta_names <- names(simulation$true$beta)
  if (is.null(beta_names)) {
    beta_names <- paste0("beta_", seq_len(ncol(beta_draws)))
  }

  rows <- lapply(seq_len(ncol(beta_draws)), function(index) {
    summarize_scalar_draws(
      beta_draws[, index],
      simulation$true$beta[index],
      beta_names[index]
    )
  })
  rows[[length(rows) + 1L]] <- summarize_scalar_draws(
    fit$sigma2[, 1], simulation$true$sigma2, "sigma2"
  )
  rows[[length(rows) + 1L]] <- summarize_scalar_draws(
    fit$sigma2_u, simulation$true$sigma2_u, "sigma2_u"
  )
  rows[[length(rows) + 1L]] <- summarize_scalar_draws(
    fit$kappa, simulation$true$kappa, "kappa"
  )
  do.call(rbind, rows)
}

cluster_parameter_summary <- function(fit, representative_partition, truth) {
  estimated_levels <- sort(unique(representative_partition))
  rows <- lapply(estimated_levels, function(cluster) {
    sites <- which(representative_partition == cluster)
    matched_true <- as.integer(names(which.max(table(truth$cluster[sites]))))
    phi_draw <- rowMeans(fit$phi[, sites, drop = FALSE])
    tau2_draw <- rowMeans(fit$tau2[, sites, drop = FALSE])
    phi_interval <- stats::quantile(phi_draw, c(0.025, 0.5, 0.975))
    tau2_interval <- stats::quantile(tau2_draw, c(0.025, 0.5, 0.975))

    data.frame(
      estimated_cluster = cluster,
      matched_true_cluster = matched_true,
      size = length(sites),
      phi_true = truth$cluster_phi[matched_true],
      phi_median = phi_interval[2],
      phi_q025 = phi_interval[1],
      phi_q975 = phi_interval[3],
      tau2_true = truth$cluster_tau2[matched_true],
      tau2_median = tau2_interval[2],
      tau2_q025 = tau2_interval[1],
      tau2_q975 = tau2_interval[3],
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}

write_trace_diagnostics <- function(fit, simulation, file) {
  truth <- simulation$true
  anchor_sites <- vapply(sort(unique(truth$cluster)), function(cluster) {
    which(truth$cluster == cluster)[1]
  }, integer(1))

  grDevices::pdf(file, width = 11, height = 8.5)
  old_par <- graphics::par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), cex.main = 0.9)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  for (index in seq_len(ncol(fit$beta))) {
    graphics::plot(
      fit$beta[, index], type = "l", xlab = "Retained draw", ylab = "",
      main = names(truth$beta)[index]
    )
    graphics::abline(h = truth$beta[index], col = "firebrick", lty = 2)
  }
  scalar_series <- list(
    sigma2 = list(draw = fit$sigma2[, 1], truth = truth$sigma2),
    sigma2_u = list(draw = fit$sigma2_u, truth = truth$sigma2_u),
    kappa = list(draw = fit$kappa, truth = truth$kappa),
    K = list(draw = fit$k_clust, truth = length(unique(truth$cluster)))
  )
  for (name in names(scalar_series)) {
    item <- scalar_series[[name]]
    graphics::plot(item$draw, type = "l", xlab = "Retained draw", ylab = "", main = name)
    graphics::abline(h = item$truth, col = "firebrick", lty = 2)
  }

  plot_anchor <- function(site) {
    graphics::plot(
      fit$phi[, site], type = "l", xlab = "Retained draw", ylab = "",
      main = paste0("phi, anchor site ", site)
    )
    graphics::abline(h = truth$phi[site], col = "firebrick", lty = 2)
    graphics::plot(
      fit$tau2[, site], type = "l", xlab = "Retained draw", ylab = "",
      main = paste0("tau2, anchor site ", site)
    )
    graphics::abline(h = truth$tau2[site], col = "firebrick", lty = 2)
  }

  plot_anchor(anchor_sites[1])
  graphics::par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), cex.main = 0.9)
  for (site in anchor_sites[-1]) {
    plot_anchor(site)
  }
  invisible(anchor_sites)
}

plot_partition_validation <- function(simulation, estimated_partition, file) {
  locations <- simulation$locations
  alignment <- align_partition_to_truth(estimated_partition, simulation$true$cluster)
  colors <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00")

  grDevices::pdf(file, width = 10, height = 4.5)
  old_par <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::plot(
    locations[, 1], locations[, 2], pch = 21, cex = 1.8,
    bg = colors[simulation$true$cluster], xlab = "x", ylab = "y",
    main = "True partition", asp = 1
  )
  graphics::grid()
  graphics::plot(
    locations[, 1], locations[, 2], pch = 21, cex = 1.8,
    bg = colors[alignment$aligned], xlab = "x", ylab = "y",
    main = "Representative posterior partition", asp = 1
  )
  graphics::grid()
  invisible(alignment)
}

plot_similarity_matrix <- function(similarity, truth, file) {
  ordering <- order(truth)
  grDevices::pdf(file, width = 6, height = 5.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::image(
    similarity[ordering, ordering],
    col = grDevices::hcl.colors(40, "YlOrRd", rev = TRUE),
    axes = FALSE,
    xlab = "Sites ordered by true regime",
    ylab = "Sites ordered by true regime",
    main = "Posterior similarity matrix"
  )
  boundaries <- cumsum(as.numeric(table(truth))) / length(truth)
  graphics::abline(v = boundaries[-length(boundaries)], lty = 2)
  graphics::abline(h = boundaries[-length(boundaries)], lty = 2)
}

plot_latent_reconstruction <- function(fit, simulation, file) {
  if (is.null(fit$w) || is.null(fit$u)) {
    return(invisible(NULL))
  }

  n <- simulation$dimensions["n"]
  n_times <- simulation$dimensions["T"]
  w_mean <- matrix(colMeans(fit$w), nrow = n_times, ncol = n, byrow = TRUE)
  u_mean <- matrix(colMeans(fit$u), nrow = n_times, ncol = n, byrow = TRUE)
  sites <- vapply(sort(unique(simulation$true$cluster)), function(cluster) {
    which(simulation$true$cluster == cluster)[1]
  }, integer(1))

  grDevices::pdf(file, width = 11, height = 7)
  old_par <- graphics::par(mfrow = c(2, length(sites)), mar = c(3, 3, 2, 1))
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  for (site in sites) {
    graphics::matplot(
      cbind(simulation$true$w[, site], w_mean[, site]), type = "l", lty = c(1, 2),
      col = c("black", "#0072B2"), xlab = "Time", ylab = "w",
      main = paste0("Site ", site, ", regime ", simulation$true$cluster[site])
    )
  }
  for (site in sites) {
    graphics::matplot(
      cbind(simulation$true$u[, site], u_mean[, site]), type = "l", lty = c(1, 2),
      col = c("black", "#D55E00"), xlab = "Time", ylab = "u",
      main = paste0("Common shock, site ", site)
    )
  }
  invisible(list(w_mean = w_mean, u_mean = u_mean, sites = sites))
}
