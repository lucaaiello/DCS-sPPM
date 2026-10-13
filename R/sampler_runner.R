global_beta_validation_prior <- function(simulation) {
  list(
    a_sigma = 3,
    b_sigma = 2 * simulation$true$sigma2,
    a_xi = 2,
    b_xi = 2,
    a_tau = 2,
    b_tau = 0.35,
    a_phi = 1,
    b_phi = 1,
    a_kappa = 3,
    b_kappa = 2,
    a_sigma2_u = 3,
    b_sigma2_u = 2 * simulation$true$sigma2_u
  )
}

run_global_beta_validation_chain <- function(
    simulation,
    runs,
    burn,
    thin = 1L,
    beta_update_mode = "conditional",
    seed = 20260911L,
    overdispersed_initialization = FALSE,
    save_latent = FALSE,
    save_w = FALSE,
    fixed_u = NULL,
    sigma2_update_mode = "conditional",
    eta_sigma2 = 0.15,
    tau2_update_mode = "conditional",
    eta_tau2 = 0.15,
    joint_variance_initial_scale = 0.15,
    spatial_update_mode = "conditional",
    spatial_initial_scale = 0.08,
    spatial_process_model = "nngp",
    model_type = "dcs_sppm",
    similarity_locations = NULL,
    prior_override = NULL,
    debug_checks = FALSE,
    progress_file = "",
    show_progress = TRUE,
    sampler_extra_args = list()
) {
  if (!exists("dcs_sppm_global_beta_cpp", mode = "function")) {
    stop("Compile `src/sampler_global_beta_single_w.cpp` before running a chain.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The `Matrix` package is required by the sparse CHOLMOD sampler.")
  }
  if (runs <= burn || burn < 0L || thin <= 0L) {
    stop("Require runs > burn >= 0 and thin > 0.")
  }

  set.seed(seed)
  n <- unname(simulation$dimensions["n"])
  p <- unname(simulation$dimensions["p"])
  prior <- global_beta_validation_prior(simulation)
  if (!is.null(prior_override)) {
    unknown_prior_names <- setdiff(names(prior_override), names(prior))
    if (length(unknown_prior_names)) {
      stop(
        "Unknown prior override(s): ",
        paste(unknown_prior_names, collapse = ", ")
      )
    }
    prior[names(prior_override)] <- prior_override
  }

  if (overdispersed_initialization) {
    init_beta <- stats::rnorm(p, 0, 0.75)
    init_sigma2 <- exp(stats::rnorm(1L, log(simulation$true$sigma2), 0.6))
    init_phi <- rep(stats::runif(1L, 0.1, 0.9), n)
    init_tau2 <- rep(exp(stats::rnorm(1L, log(0.35), 0.5)), n)
    init_kappa <- exp(stats::rnorm(1L, log(simulation$true$kappa), 0.4))
    init_sigma2_u <- exp(
      stats::rnorm(1L, log(simulation$true$sigma2_u), 0.5)
    )
  } else {
    init_beta <- rep(0, p)
    init_sigma2 <- simulation$true$sigma2
    init_phi <- rep(0.5, n)
    init_tau2 <- rep(0.35, n)
    init_kappa <- simulation$true$kappa
    init_sigma2_u <- simulation$true$sigma2_u
  }
  if (is.null(similarity_locations)) {
    similarity_locations <- scale(simulation$locations)
  }

  sampler_arguments <- list(
      Y = simulation$Y,
      locations = simulation$locations,
      X = simulation$X,
      runs = as.integer(runs),
      burn = as.integer(burn),
      thin = as.integer(thin),
      cohesion = 3L,
      init_beta = init_beta,
      init_xi2 = rep(2, p),
      init_sigma2 = init_sigma2,
      init_phi = init_phi,
      init_tau2 = init_tau2,
      init_M = 0.5,
      init_kappa = init_kappa,
      nu = simulation$true$nu,
      NN_ind = simulation$NN_ind,
      k0 = 1,
      v0 = 4,
      a_sigma = prior$a_sigma,
      b_sigma = prior$b_sigma,
      a_xi = prior$a_xi,
      b_xi = prior$b_xi,
      a_tau = prior$a_tau,
      b_tau = prior$b_tau,
      a_phi = prior$a_phi,
      b_phi = prior$b_phi,
      a_kappa = prior$a_kappa,
      b_kappa = prior$b_kappa,
      eta_phi = 0.05,
      eta_kappa = 0.05,
      beta_update_mode = beta_update_mode,
      init_sigma2_u = init_sigma2_u,
      a_sigma2_u = prior$a_sigma2_u,
      b_sigma2_u = prior$b_sigma2_u,
      debug_checks = debug_checks,
      locations_similarity = similarity_locations,
      save_latent = save_latent,
      save_w = save_w,
      fixed_u = fixed_u,
      sigma2_update_mode = sigma2_update_mode,
      eta_sigma2 = eta_sigma2,
      tau2_update_mode = tau2_update_mode,
      eta_tau2 = eta_tau2,
      joint_variance_initial_scale = joint_variance_initial_scale,
      spatial_update_mode = spatial_update_mode,
      spatial_initial_scale = spatial_initial_scale,
      spatial_process_model = spatial_process_model
  )
  if (!is.null(model_type)) {
    sampler_arguments$model_type <- model_type
  }
  sampler_formals <- names(formals(dcs_sppm_global_beta_cpp))
  progress_requested <- nzchar(progress_file) || !isTRUE(show_progress)
  progress_supported <- all(
    c("progress_file", "show_progress") %in% sampler_formals
  )
  if (progress_requested && !progress_supported) {
    stop("The compiled sampler does not support the requested progress dashboard.")
  }
  if (progress_supported) {
    sampler_arguments$progress_file <- progress_file
    sampler_arguments$show_progress <- show_progress
  }
  if (!is.list(sampler_extra_args)) {
    stop("`sampler_extra_args` must be a named list.")
  }
  if (length(sampler_extra_args)) {
    if (is.null(names(sampler_extra_args)) ||
        any(!nzchar(names(sampler_extra_args)))) {
      stop("Every entry of `sampler_extra_args` must be named.")
    }
    duplicated_arguments <- intersect(
      names(sampler_extra_args),
      names(sampler_arguments)
    )
    if (length(duplicated_arguments)) {
      stop(
        "Duplicated sampler argument(s): ",
        paste(duplicated_arguments, collapse = ", ")
      )
    }
    sampler_arguments <- c(sampler_arguments, sampler_extra_args)
  }

  elapsed <- system.time({
    fit <- do.call(dcs_sppm_global_beta_cpp, sampler_arguments)
  })

  colnames(fit$beta) <- names(simulation$true$beta)
  fit$elapsed_seconds <- unname(elapsed["elapsed"])
  fit$mcmc_settings <- list(
    runs = as.integer(runs),
    burn = as.integer(burn),
    thin = as.integer(thin),
    beta_update_mode = beta_update_mode,
    sigma2_update_mode = sigma2_update_mode,
    tau2_update_mode = tau2_update_mode,
    spatial_update_mode = spatial_update_mode,
    spatial_process_model = spatial_process_model,
    model_type = if (is.null(model_type)) "cpp_default" else model_type,
    prior_override = prior_override,
    seed = as.integer(seed),
    overdispersed_initialization = overdispersed_initialization,
    condition_on_fixed_u = !is.null(fixed_u)
  )
  fit$prior <- prior
  fit
}
