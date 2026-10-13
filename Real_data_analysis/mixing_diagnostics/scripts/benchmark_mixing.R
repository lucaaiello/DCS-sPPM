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
source(file.path(validation_root, "R", "validation_utils.R"))

required_packages <- c("Rcpp", "RcppArmadillo", "Matrix", "data.table")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install the required package(s): ", paste(missing_packages, collapse = ", "))
}

runs <- as.integer(Sys.getenv("DCS_MIXING_BENCHMARK_RUNS", "500"))
burn <- as.integer(Sys.getenv("DCS_MIXING_BENCHMARK_BURN", "250"))
thin <- as.integer(Sys.getenv("DCS_MIXING_BENCHMARK_THIN", "1"))
threads <- as.integer(Sys.getenv(
  "DCS_MIXING_BENCHMARK_THREADS",
  as.character(min(6L, parallel::detectCores(logical = FALSE)))
))
seed <- as.integer(Sys.getenv("DCS_MIXING_BENCHMARK_SEED", "20260923"))
if (anyNA(c(runs, burn, thin, threads, seed)) || runs <= burn || burn < 0L ||
    thin < 1L || threads < 1L) {
  stop("Require runs > burn >= 0, thin >= 1, and threads >= 1.")
}
Sys.setenv(OMP_NUM_THREADS = threads)

benchmark_name <- Sys.getenv(
  "DCS_MIXING_BENCHMARK_NAME",
  sprintf("mixing_benchmark_%d_confounding_diagnostics", runs)
)
if (!grepl("^[A-Za-z0-9_.-]+$", benchmark_name) ||
    benchmark_name %in% c(".", "..")) {
  stop("`DCS_MIXING_BENCHMARK_NAME` must be a safe directory name.")
}
benchmark_dir <- file.path(
  validation_results_dir,
  "Real_data_analysis", "mixing_diagnostics", "results", benchmark_name
)
dir.create(benchmark_dir, recursive = TRUE, showWarnings = FALSE)

preparation_file <- file.path(validation_root, "Data", "production_10000_single_chain_preparation.rds")
if (!file.exists(preparation_file)) {
  stop("Frozen diagnostic preparation input not found: ", preparation_file)
}
preparation <- readRDS(preparation_file)
station_information <- data.table::as.data.table(
  preparation$station_information
)
analysis_dates <- as.Date(preparation$analysis_dates)
data_scaling <- preparation$data_scaling
initial_values <- preparation$initial_values
if (is.null(initial_values)) stop("Missing frozen diagnostic starting values.")
names(initial_values) <- paste0("init_", names(initial_values))
rm(preparation)
invisible(gc())

pm10 <- data.table::as.data.table(readRDS(
  file.path(validation_root, "Data", "pm10_combined_2024.rds")
))
pm10[, Date := as.Date(as.numeric(unclass(Date)), origin = "1970-01-01")]
pm10 <- pm10[
  Samplingpoint %in% station_information$Samplingpoint & Date %in% analysis_dates
]
Y_raw <- matrix(
  NA_real_,
  nrow = length(analysis_dates),
  ncol = nrow(station_information),
  dimnames = list(
    as.character(analysis_dates),
    station_information$Samplingpoint
  )
)
date_index <- match(pm10$Date, analysis_dates)
station_index <- match(pm10$Samplingpoint, station_information$Samplingpoint)
if (anyDuplicated(cbind(date_index, station_index))) {
  stop("Duplicate station-date observations remain in the benchmark data.")
}
Y_raw[cbind(date_index, station_index)] <- pm10$Value
impute_short_gaps <- function(values) {
  observed <- which(is.finite(values))
  if (length(observed) < 2L) {
    stop("A retained station has fewer than two finite observations.")
  }
  stats::approx(
    x = observed,
    y = values[observed],
    xout = seq_along(values),
    method = "linear",
    rule = 2
  )$y
}
Y_raw <- apply(Y_raw, 2L, impute_short_gaps)
Y <- (Y_raw - data_scaling$response_center) / data_scaling$response_scale
rm(pm10, Y_raw)
invisible(gc())

weather <- data.table::as.data.table(readRDS(
  file.path(validation_root, "Data", "era5_combined_2024.rds")
))
weather[, Date := as.Date("2024-01-01") + Date]
weather[, Longitude := as.numeric(Longitude)]
weather[, Latitude := as.numeric(Latitude)]
weather_cells <- unique(station_information[, .(
  weather_longitude,
  weather_latitude,
  weather_cell
)])
weather <- weather[
  weather_cells,
  on = c(
    "Longitude" = "weather_longitude",
    "Latitude" = "weather_latitude"
  ),
  nomatch = 0L,
  allow.cartesian = TRUE
][Date %in% analysis_dates]
weather_matrix <- function(variable) {
  result <- matrix(
    NA_real_,
    nrow = length(analysis_dates),
    ncol = max(weather_cells$weather_cell)
  )
  result[cbind(
    match(weather$Date, analysis_dates),
    weather$weather_cell
  )] <- weather[[variable]]
  result[, station_information$weather_cell, drop = FALSE]
}
temperature <- weather_matrix("Temperature")
precipitation <- weather_matrix("Precipitation")
wind_speed <- weather_matrix("Wind")
if (!all(is.finite(temperature)) || !all(is.finite(precipitation)) ||
    !all(is.finite(wind_speed))) {
  stop("Weather matching produced missing benchmark covariates.")
}
time_major_vector <- function(values) as.vector(t(values))
X <- cbind(
  intercept = 1,
  temperature = time_major_vector(
    (temperature - data_scaling$temperature_center) /
      data_scaling$temperature_scale
  ),
  precipitation = time_major_vector(
    (precipitation - data_scaling$precipitation_center) /
      data_scaling$precipitation_scale
  ),
  wind_speed = time_major_vector(
    (wind_speed - data_scaling$wind_speed_center) /
      data_scaling$wind_speed_scale
  )
)
rm(weather, temperature, precipitation, wind_speed)
invisible(gc())

latitude_reference <- mean(station_information$Latitude) * pi / 180
locations <- cbind(
  x_1000km = station_information$Longitude *
    cos(latitude_reference) * 111.32 / 1000,
  y_1000km = station_information$Latitude * 110.57 / 1000
)
NN_ind <- make_previous_neighbor_matrix(locations, 10L)

sampler_arguments <- c(
  list(
    Y = Y,
    locations = locations,
    X = X,
    runs = runs,
    burn = burn,
    thin = thin,
    cohesion = 3L
  ),
  initial_values,
  list(
    nu = 0.50,
    NN_ind = NN_ind,
    k0 = 1,
    v0 = 4,
    a_sigma = 3,
    b_sigma = 0.40,
    a_xi = 2,
    b_xi = 2,
    a_tau = 2,
    b_tau = 0.35,
    a_phi = 1,
    b_phi = 1,
    a_kappa = 3,
    b_kappa = 2,
    eta_phi = 0.05,
    eta_kappa = 0.05,
    CC = 3L,
    beta_update_mode = "collapsed",
    a_sigma2_u = 3,
    b_sigma2_u = 0.60,
    locations_similarity = scale(locations),
    save_latent = FALSE,
    save_w = FALSE,
    sigma2_update_mode = "collapsed",
    tau2_update_mode = "pairwise_interwoven",
    joint_variance_initial_scale = 0.03,
    spatial_update_mode = "collapsed",
    spatial_initial_scale = 0.00003,
    spatial_process_model = "nngp",
    show_progress = TRUE
  )
)

diagnostic_series <- function(fit) {
  values <- cbind(
    fit$beta,
    sigma2 = as.numeric(fit$sigma2),
    kappa = as.numeric(fit$kappa),
    sigma2_u = as.numeric(fit$sigma2_u),
    Mdp = as.numeric(fit$Mdp),
    K = as.numeric(fit$k_clust)
  )
  colnames(values)[seq_len(ncol(fit$beta))] <- paste0(
    "beta_", colnames(X)
  )
  values
}

summarize_diagnostics <- function(fit, kernel) {
  values <- diagnostic_series(fit)
  rows <- lapply(seq_len(ncol(values)), function(column) {
    draws <- values[, column]
    ess <- effective_sample_size(draws)
    data.frame(
      kernel = kernel,
      parameter = colnames(values)[column],
      posterior_mean = mean(draws),
      posterior_sd = stats::sd(draws),
      lag1_autocorrelation = if (length(draws) > 1L) {
        stats::cor(draws[-length(draws)], draws[-1L])
      } else {
        NA_real_
      },
      ess = ess,
      mcse = stats::sd(draws) / sqrt(ess),
      ess_per_second = ess / fit$elapsed_seconds,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_kernel <- function(fit, kernel) {
  cluster_change <- if (nrow(fit$cluster) > 1L) {
    fit$cluster[-1L, , drop = FALSE] !=
      fit$cluster[-nrow(fit$cluster), , drop = FALSE]
  } else {
    matrix(FALSE, nrow = 0L, ncol = ncol(fit$cluster))
  }
  data.frame(
    kernel = kernel,
    elapsed_seconds = fit$elapsed_seconds,
    seconds_per_iteration = fit$elapsed_seconds / runs,
    mean_K = mean(fit$k_clust),
    sd_K = stats::sd(fit$k_clust),
    proportion_iterations_partition_changed = if (nrow(cluster_change)) {
      mean(rowSums(cluster_change) > 0L)
    } else {
      NA_real_
    },
    mean_proportion_sites_changed = if (nrow(cluster_change)) {
      mean(cluster_change)
    } else {
      NA_real_
    },
    phi_acceptance = fit$accphi,
    pairwise_variance_acceptance = fit$acc_pairwise_variance,
    spatial_acceptance = fit$acc_spatial,
    split_acceptance = if (!is.null(fit$acc_split)) fit$acc_split else NA_real_,
    merge_acceptance = if (!is.null(fit$acc_merge)) fit$acc_merge else NA_real_,
    split_proposals = if (!is.null(fit$n_split_prop)) fit$n_split_prop else 0L,
    merge_proposals = if (!is.null(fit$n_merge_prop)) fit$n_merge_prop else 0L,
    stringsAsFactors = FALSE
  )
}

summarize_latent_projections <- function(fit, kernel) {
  projection_sources <- list(
    w = fit$w_design_projection,
    u_transition = fit$u_transition_design_projection
  )
  if (any(vapply(projection_sources, is.null, logical(1L)))) {
    stop("The diagnostic sampler did not return latent projections.")
  }
  rows <- list()
  row_index <- 1L
  for (source_name in names(projection_sources)) {
    projection <- projection_sources[[source_name]]
    for (column in seq_len(ncol(projection))) {
      draws <- projection[, column]
      ess <- effective_sample_size(draws)
      rows[[row_index]] <- data.frame(
        kernel = kernel,
        latent_component = source_name,
        covariate = colnames(X)[column],
        posterior_mean = mean(draws),
        posterior_sd = stats::sd(draws),
        lag1_autocorrelation = stats::cor(
          draws[-length(draws)],
          draws[-1L]
        ),
        ess = ess,
        mcse = stats::sd(draws) / sqrt(ess),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  do.call(rbind, rows)
}

summarize_beta_latent_correlations <- function(fit, kernel) {
  projection_sources <- list(
    w = fit$w_design_projection,
    u_transition = fit$u_transition_design_projection
  )
  rows <- list()
  row_index <- 1L
  for (source_name in names(projection_sources)) {
    projection <- projection_sources[[source_name]]
    for (beta_column in seq_len(ncol(fit$beta))) {
      for (latent_column in seq_len(ncol(projection))) {
        rows[[row_index]] <- data.frame(
          kernel = kernel,
          beta = colnames(X)[beta_column],
          latent_component = source_name,
          latent_covariate = colnames(X)[latent_column],
          posterior_correlation = stats::cor(
            fit$beta[, beta_column],
            projection[, latent_column]
          ),
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L
      }
    }
  }
  do.call(rbind, rows)
}

run_kernel <- function(kernel, source_file) {
  kernel_dir <- file.path(benchmark_dir, kernel)
  dir.create(kernel_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(
    dirname(validation_root),
    "tmp",
    paste0("sourcecpp_", kernel)
  )
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(
    file.path(validation_root, "src", source_file),
    cacheDir = cache_dir,
    rebuild = FALSE,
    showOutput = TRUE
  )
  sampler <- dcs_sppm_global_beta_cpp
  arguments <- sampler_arguments
  arguments$progress_file <- file.path(kernel_dir, "progress.txt")
  set.seed(seed)
  elapsed <- system.time({
    fit <- do.call(sampler, arguments)
  })
  fit$elapsed_seconds <- unname(elapsed["elapsed"])
  fit$benchmark_settings <- list(
    runs = runs,
    burn = burn,
    thin = thin,
    threads = threads,
    seed = seed,
    source_file = source_file
  )
  diagnostics <- summarize_diagnostics(fit, kernel)
  summary <- summarize_kernel(fit, kernel)
  latent_projections <- summarize_latent_projections(fit, kernel)
  beta_latent_correlations <- summarize_beta_latent_correlations(fit, kernel)
  saveRDS(fit, file.path(kernel_dir, "fit.rds"), compress = FALSE)
  utils::write.csv(
    fit$block_timing,
    file.path(kernel_dir, "block_timing.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    diagnostics,
    file.path(kernel_dir, "diagnostics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summary,
    file.path(kernel_dir, "kernel_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    latent_projections,
    file.path(kernel_dir, "latent_projection_diagnostics.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    beta_latent_correlations,
    file.path(kernel_dir, "beta_latent_correlations.csv"),
    row.names = FALSE
  )
  rm(fit)
  invisible(gc())
  list(
    summary = summary,
    diagnostics = diagnostics,
    latent_projections = latent_projections,
    beta_latent_correlations = beta_latent_correlations
  )
}

original_results_name <- Sys.getenv(
  "DCS_MIXING_ORIGINAL_RESULTS_NAME",
  "mixing_benchmark_500"
)
original_results_dir <- file.path(
  validation_results_dir,
  "Real_data_analysis", "mixing_diagnostics", "results", original_results_name,
  "original"
)
original_summary_file <- file.path(
  original_results_dir,
  "kernel_summary.csv"
)
original_diagnostics_file <- file.path(
  original_results_dir,
  "diagnostics.csv"
)
if (!file.exists(original_summary_file) ||
    !file.exists(original_diagnostics_file)) {
  stop(
    "Completed original-kernel benchmark summaries were not found under: ",
    original_results_dir
  )
}

kernel_results <- list(
  original = list(
    summary = utils::read.csv(
      original_summary_file,
      stringsAsFactors = FALSE
    ),
    diagnostics = utils::read.csv(
      original_diagnostics_file,
      stringsAsFactors = FALSE
    )
  ),
  confounding_diagnostics = run_kernel(
    "confounding_diagnostics",
    "sampler_global_beta_real_data_2.cpp"
  )
)

kernel_summary <- do.call(rbind, lapply(kernel_results, `[[`, "summary"))
diagnostic_summary <- do.call(
  rbind,
  lapply(kernel_results, `[[`, "diagnostics")
)
latent_projection_summary <-
  kernel_results$confounding_diagnostics$latent_projections
beta_latent_correlation_summary <-
  kernel_results$confounding_diagnostics$beta_latent_correlations
utils::write.csv(
  kernel_summary,
  file.path(benchmark_dir, "kernel_comparison.csv"),
  row.names = FALSE
)
utils::write.csv(
  diagnostic_summary,
  file.path(benchmark_dir, "diagnostic_comparison.csv"),
  row.names = FALSE
)
utils::write.csv(
  latent_projection_summary,
  file.path(benchmark_dir, "latent_projection_diagnostics.csv"),
  row.names = FALSE
)
utils::write.csv(
  beta_latent_correlation_summary,
  file.path(benchmark_dir, "beta_latent_correlations.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    kernel_summary = kernel_summary,
    diagnostic_summary = diagnostic_summary,
    latent_projection_summary = latent_projection_summary,
    beta_latent_correlation_summary = beta_latent_correlation_summary,
    settings = list(
      runs = runs,
      burn = burn,
      thin = thin,
      threads = threads,
      seed = seed
    )
  ),
  file.path(benchmark_dir, "benchmark_summary.rds")
)

print(kernel_summary)
print(diagnostic_summary)
message("Benchmark results written to: ", benchmark_dir)
