make_previous_neighbor_matrix <- function(locations, n_neighbors = 8L) {
  locations <- as.matrix(locations)
  n <- nrow(locations)
  n_neighbors <- min(as.integer(n_neighbors), n - 1L)

  if (n < 2L || ncol(locations) != 2L || n_neighbors < 1L) {
    stop("Require an n x 2 location matrix, n >= 2, and at least one neighbor.")
  }

  neighbor_index <- matrix(0L, nrow = n, ncol = n_neighbors)
  for (site in 2:n) {
    candidates <- seq_len(site - 1L)
    distances <- sqrt(rowSums(
      (locations[candidates, , drop = FALSE] -
         matrix(locations[site, ], nrow = length(candidates), ncol = 2L,
                byrow = TRUE))^2
    ))
    selected <- candidates[order(distances)][seq_len(min(n_neighbors, site - 1L))]
    neighbor_index[site, seq_along(selected)] <- selected
  }
  neighbor_index
}

simulate_nngp_common_shocks <- function(
    locations,
    neighbor_index,
    n_times,
    sigma2_u,
    kappa,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  locations <- as.matrix(locations)
  neighbor_index <- as.matrix(neighbor_index)
  n <- nrow(locations)
  common_shock <- matrix(0, nrow = n_times, ncol = n)
  conditional_coefficients <- vector("list", n)
  conditional_variances <- numeric(n)

  for (site in seq_len(n)) {
    neighbors <- neighbor_index[site, ]
    neighbors <- neighbors[neighbors > 0L]

    if (!length(neighbors)) {
      conditional_coefficients[[site]] <- numeric(0)
      conditional_variances[site] <- sigma2_u
      next
    }

    neighbor_distances <- as.matrix(dist(locations[neighbors, , drop = FALSE]))
    neighbor_covariance <- sigma2_u * exp(-neighbor_distances / kappa)
    cross_distances <- sqrt(rowSums(
      (locations[neighbors, , drop = FALSE] -
         matrix(locations[site, ], nrow = length(neighbors), ncol = 2L,
                byrow = TRUE))^2
    ))
    cross_covariance <- sigma2_u * exp(-cross_distances / kappa)
    coefficients <- solve(neighbor_covariance, cross_covariance)
    conditional_variance <- sigma2_u - sum(cross_covariance * coefficients)

    if (!is.finite(conditional_variance) || conditional_variance <= 0) {
      stop("Invalid NNGP conditional variance while simulating site ", site, ".")
    }

    conditional_coefficients[[site]] <- as.numeric(coefficients)
    conditional_variances[site] <- conditional_variance
  }

  for (time in 2:n_times) {
    for (site in seq_len(n)) {
      neighbors <- neighbor_index[site, ]
      neighbors <- neighbors[neighbors > 0L]
      conditional_mean <- if (length(neighbors)) {
        sum(conditional_coefficients[[site]] * common_shock[time, neighbors])
      } else {
        0
      }
      common_shock[time, site] <- stats::rnorm(
        1L,
        conditional_mean,
        sqrt(conditional_variances[site])
      )
    }
  }

  list(
    u = common_shock,
    coefficients = conditional_coefficients,
    conditional_variances = conditional_variances
  )
}

simulate_exact_gp_common_shocks <- function(
    locations,
    n_times,
    sigma2_u,
    kappa,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  locations <- as.matrix(locations)
  n <- nrow(locations)
  covariance <- sigma2_u * exp(
    -as.matrix(stats::dist(locations)) / kappa
  )
  covariance <- covariance + diag(.Machine$double.eps * sigma2_u, n)
  upper_cholesky <- chol(covariance)
  common_shock <- matrix(0, nrow = n_times, ncol = n)
  common_shock[2:n_times, ] <-
    matrix(stats::rnorm((n_times - 1L) * n), nrow = n_times - 1L) %*%
    upper_cholesky

  list(u = common_shock, covariance = covariance)
}

standardize_covariate <- function(x) {
  scale_value <- stats::sd(x)
  if (!is.finite(scale_value) || scale_value <= 0) {
    stop("Cannot standardize a constant or non-finite covariate.")
  }
  as.numeric((x - mean(x)) / scale_value)
}

simulate_dcs_sppm_global_beta <- function(
    n_side = 6L,
    n_times = 120L,
    n_neighbors = 8L,
    beta = c(1.0, 1.25, -0.8),
    sigma2 = 0.08,
    sigma2_u = 0.35,
    kappa = 1.5,
    nu = 0.5,
    cluster_phi = c(0.15, 0.55, 0.88),
    cluster_tau2 = c(0.12, 0.35, 0.75),
    spatial_dgp = c("nngp", "exact_gp"),
    seed = 20260910L,
    locations = NULL,
    true_cluster = NULL
) {
  set.seed(seed)
  spatial_dgp <- match.arg(spatial_dgp)

  if (nu != 0.5) {
    stop("The first simulator currently uses the exponential Matern case nu = 0.5.")
  }
  if (length(cluster_phi) != 3L || length(cluster_tau2) != 3L) {
    stop("The first validation scenario requires exactly three regimes.")
  }

  if (is.null(locations)) {
    coordinates <- expand.grid(
      x = seq(0, 5, length.out = n_side),
      y = seq(0, 5, length.out = n_side),
      KEEP.OUT.ATTRS = FALSE
    )
    coordinates <- coordinates[
      order(coordinates$x, coordinates$y),
      ,
      drop = FALSE
    ]
    rownames(coordinates) <- NULL
    locations <- as.matrix(coordinates)
  } else {
    locations <- as.matrix(locations)
    if (ncol(locations) != 2L || nrow(locations) < 3L ||
        !all(is.finite(locations))) {
      stop("`locations` must be a finite matrix with two columns and at least three rows.")
    }
    colnames(locations) <- c("x", "y")
  }
  n <- nrow(locations)

  if (is.null(true_cluster)) {
    x_breaks <- stats::quantile(locations[, "x"], c(1 / 3, 2 / 3))
    true_cluster <- as.integer(cut(
      locations[, "x"],
      breaks = c(-Inf, x_breaks, Inf),
      labels = FALSE,
      include.lowest = TRUE
    ))
  } else {
    true_cluster <- as.integer(true_cluster)
    if (length(true_cluster) != n || anyNA(true_cluster) ||
        !identical(sort(unique(true_cluster)), 1:3)) {
      stop("`true_cluster` must have one entry per site and use labels 1, 2, and 3.")
    }
  }
  phi <- cluster_phi[true_cluster]
  tau2 <- cluster_tau2[true_cluster]

  neighbor_index <- make_previous_neighbor_matrix(locations, n_neighbors)
  spatial_process <- if (spatial_dgp == "nngp") {
    simulate_nngp_common_shocks(
      locations = locations,
      neighbor_index = neighbor_index,
      n_times = n_times,
      sigma2_u = sigma2_u,
      kappa = kappa
    )
  } else {
    simulate_exact_gp_common_shocks(
      locations = locations,
      n_times = n_times,
      sigma2_u = sigma2_u,
      kappa = kappa
    )
  }

  time_index <- rep(seq_len(n_times), each = n)
  site_index <- rep(seq_len(n), times = n_times)
  x_coordinate <- standardize_covariate(locations[, "x"])
  y_coordinate <- standardize_covariate(locations[, "y"])

  raw_covariate_1 <-
    sin(2 * pi * time_index / 24) +
    0.35 * x_coordinate[site_index] +
    stats::rnorm(n * n_times, sd = 0.55)
  raw_covariate_2 <-
    cos(2 * pi * time_index / 17) * y_coordinate[site_index] +
    stats::rnorm(n * n_times, sd = 0.65)

  X <- cbind(
    intercept = 1,
    space_time_1 = standardize_covariate(raw_covariate_1),
    space_time_2 = standardize_covariate(raw_covariate_2)
  )
  if (ncol(X) != length(beta)) {
    stop("The global beta vector must have one entry per design column.")
  }

  eta <- matrix(0, nrow = n_times, ncol = n)
  w <- matrix(0, nrow = n_times, ncol = n)
  w[1, ] <- stats::rnorm(n, sd = sqrt(tau2 / (1 - phi^2)))
  for (time in 2:n_times) {
    eta[time, ] <- stats::rnorm(n, sd = sqrt(tau2))
    w[time, ] <-
      phi * w[time - 1L, ] + spatial_process$u[time, ] + eta[time, ]
  }

  fixed_effect <- matrix(X %*% beta, nrow = n_times, ncol = n, byrow = TRUE)
  epsilon <- matrix(
    stats::rnorm(n * n_times, sd = sqrt(sigma2)),
    nrow = n_times,
    ncol = n
  )
  Y <- fixed_effect + w + epsilon

  colnames(Y) <- paste0("site_", seq_len(n))
  colnames(w) <- colnames(Y)
  colnames(spatial_process$u) <- colnames(Y)

  list(
    Y = Y,
    X = X,
    locations = locations,
    NN_ind = neighbor_index,
    true = list(
      beta = setNames(beta, colnames(X)),
      sigma2 = sigma2,
      sigma2_u = sigma2_u,
      kappa = kappa,
      nu = nu,
      spatial_dgp = spatial_dgp,
      cluster = true_cluster,
      phi = phi,
      tau2 = tau2,
      cluster_phi = cluster_phi,
      cluster_tau2 = cluster_tau2,
      w = w,
      u = spatial_process$u,
      eta = eta,
      fixed_effect = fixed_effect
    ),
    dimensions = c(n = n, T = n_times, p = ncol(X)),
    seed = seed,
    design_order = "time-major: row t*n+i contains x_t(s_i)"
  )
}
