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
source(file.path(validation_root, "R", "model_comparison_utils.R"))
source(file.path(validation_root, "R", "validation_utils.R"))
if (!requireNamespace("salso", quietly = TRUE)) stop("Install the salso package.")

results_name <- Sys.getenv(
  "DCS_REAL_FULL_RESULTS_NAME", "production_30000_unique_locations"
)
root <- file.path(
  validation_results_dir, "Real_data_analysis", "model_comparison", "results", results_name
)
dir.create(root, recursive = TRUE, showWarnings = FALSE)
files <- list.files(
  root, pattern = "^(pilot_fit|fit)[.]rds$", recursive = TRUE,
  full.names = TRUE
)
fits <- lapply(files, readRDS)
if (!length(fits)) stop("No full-data fits found.")
if (!any(vapply(fits, function(fit) identical(fit$model, "dcs_sppm"), logical(1L)))) {
  stop("The comparison requires a DCS-sPPM fit in the same result directory; historical fits are not substituted.")
}
if (nzchar(Sys.getenv("DCS_REAL_FULL_DATA_FILE", ""))) {
  if (!setequal(vapply(fits, function(fit) fit$model, character(1L)),
      c("dcs_sppm", "ppm_common_shock", "k1_common_shock", "sppm_no_common_shock"))) {
    stop("The unique-location comparison requires all four new fits.")
  }
  reference_inventory <- fits[[1L]]$station_information$Samplingpoint
  if (!all(vapply(fits, function(fit) identical(
      fit$station_information$Samplingpoint, reference_inventory), logical(1L)))) {
    stop("Comparison fits have different station inventories.")
  }
  if (!all(vapply(fits, function(fit) identical(fit$comparison_settings$data_md5,
      unname(tools::md5sum(Sys.getenv("DCS_REAL_FULL_DATA_FILE")))), logical(1L)))) {
    stop("Comparison fits do not match the selected dataset checksum.")
  }
}

vi_runs <- as.integer(Sys.getenv("DCS_COMPARISON_VI_RUNS", "32"))
vi_cores <- as.integer(Sys.getenv("DCS_COMPARISON_VI_CORES", "8"))
set.seed(as.integer(Sys.getenv("DCS_COMPARISON_VI_SEED", "20261002")))
`%||%` <- function(left, right) if (is.null(left)) right else left
summary_rows <- list()
size_rows <- list()
parameter_rows <- list()
runtime_rows <- list()
prediction_rows <- list()
chain_rows <- list()
diagnostic_rows <- list()

append_chain_summary <- function(rows, fit, model, quantity, values, labels = NULL) {
  values <- as.matrix(values)
  if (is.null(labels)) labels <- seq_len(ncol(values))
  for (column in seq_len(ncol(values))) {
    chain <- values[, column]
    if (!any(is.finite(chain))) next
    rows[[length(rows) + 1L]] <- data.frame(
      model = model,
      model_label = dcs_model_label(model),
      quantity = quantity,
      component = as.character(labels[column]),
      mean = mean(chain),
      sd = stats::sd(chain),
      q025 = unname(stats::quantile(chain, 0.025)),
      median = stats::median(chain),
      q975 = unname(stats::quantile(chain, 0.975)),
      ess = effective_sample_size(chain),
      stringsAsFactors = FALSE
    )
  }
  rows
}

for (index in seq_along(fits)) {
  fit <- fits[[index]]
  model <- validate_dcs_model_type(fit$model %||% "dcs_sppm")
  uses_partition <- dcs_model_uses(model, "uses_partition")
  if (uses_partition) {
    vi_partition <- as.integer(salso::salso(
      fit$cluster, loss = salso::VI(), nRuns = vi_runs,
      nCores = vi_cores
    ))
    vi_partition <- match(vi_partition, unique(vi_partition))
    cluster_sizes <- as.integer(table(vi_partition))
    vi_k <- length(cluster_sizes)
    singleton_count <- sum(cluster_sizes == 1L)
    station_information <- fit$station_information
    latitude_reference <- mean(station_information$Latitude) * pi / 180
    coordinates <- cbind(
      station_information$Longitude * cos(latitude_reference) * 111.32,
      station_information$Latitude * 110.57
    )
    within_distances <- unlist(lapply(split(seq_len(nrow(coordinates)), vi_partition), function(members) {
      if (length(members) < 2L) return(numeric())
      as.numeric(stats::dist(coordinates[members, , drop = FALSE]))
    }), use.names = FALSE)
    mean_within_distance <- mean(within_distances)
    median_within_distance <- stats::median(within_distances)
    size_rows[[length(size_rows) + 1L]] <- data.frame(
      model = model, cluster = seq_along(cluster_sizes), size = cluster_sizes
    )
  } else {
    vi_k <- 1L
    singleton_count <- 0L
    cluster_sizes <- ncol(fit$cluster)
    mean_within_distance <- NA_real_
    median_within_distance <- NA_real_
    size_rows[[length(size_rows) + 1L]] <- data.frame(
      model = model, cluster = 1L, size = cluster_sizes
    )
  }
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    model = model,
    model_label = dcs_model_label(model),
    posterior_mean_K = mean(fit$k_clust),
    posterior_median_K = stats::median(fit$k_clust),
    K_q025 = unname(stats::quantile(fit$k_clust, 0.025)),
    K_q975 = unname(stats::quantile(fit$k_clust, 0.975)),
    VI_optimal_K = vi_k,
    VI_singletons = singleton_count,
    mean_VI_cluster_size = mean(cluster_sizes),
    median_VI_cluster_size = stats::median(cluster_sizes),
    mean_within_cluster_distance_km = mean_within_distance,
    median_within_cluster_distance_km = median_within_distance,
    stringsAsFactors = FALSE
  )
  parameter_rows[[length(parameter_rows) + 1L]] <- data.frame(
    model = model,
    phi_mean = mean(unlist(fit$phistar)),
    phi_sd = stats::sd(unlist(fit$phistar)),
    tau2_mean = mean(unlist(fit$tau2star)),
    tau2_sd = stats::sd(unlist(fit$tau2star)),
    stringsAsFactors = FALSE
  )
  runtime_rows[[length(runtime_rows) + 1L]] <- data.frame(
    model = model,
    model_label = dcs_model_label(model),
    elapsed_seconds = fit$elapsed_seconds %||% NA_real_,
    prediction_seconds = if (!is.null(fit$prediction)) {
      fit$prediction$seconds %||% NA_real_
    } else {
      NA_real_
    },
    retained_draws = nrow(fit$cluster),
    prediction_draws = if (!is.null(fit$prediction)) {
      fit$prediction$retained_draws %||% NA_integer_
    } else {
      NA_integer_
    },
    stringsAsFactors = FALSE
  )
  if (!is.null(fit$prediction)) {
    for (process in intersect(c("u", "w", "y"), names(fit$prediction))) {
      process_summary <- fit$prediction[[process]]
      prediction_rows[[length(prediction_rows) + 1L]] <- data.frame(
        model = model,
        model_label = dcs_model_label(model),
        process = process,
        mean_surface_mean = mean(process_summary$mean),
        mean_surface_sd = stats::sd(as.numeric(process_summary$mean)),
        mean_posterior_sd = mean(process_summary$sd),
        median_posterior_sd = stats::median(process_summary$sd),
        structurally_zero = model == "sppm_no_common_shock" && process == "u",
        stringsAsFactors = FALSE
      )
    }
  }
  chain_rows <- append_chain_summary(
    chain_rows, fit, model, "beta", fit$beta,
    colnames(fit$beta) %||% seq_len(ncol(fit$beta))
  )
  chain_rows <- append_chain_summary(
    chain_rows, fit, model, "regression_prior_variance", fit$sigma2beta,
    colnames(fit$sigma2beta) %||% seq_len(ncol(fit$sigma2beta))
  )
  chain_rows <- append_chain_summary(
    chain_rows, fit, model, "measurement_variance", fit$sigma2, "sigma2"
  )
  if (uses_partition) {
    chain_rows <- append_chain_summary(
      chain_rows, fit, model, "partition_mass", fit$Mdp, "M"
    )
  }
  if (dcs_model_uses(model, "uses_common_shock")) {
    chain_rows <- append_chain_summary(
      chain_rows, fit, model, "spatial_range", fit$kappa, "kappa"
    )
    chain_rows <- append_chain_summary(
      chain_rows, fit, model, "common_shock_variance", fit$sigma2_u,
      "sigma2_u"
    )
  }
  diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
    model = model,
    model_label = dcs_model_label(model),
    phi_acceptance = fit$accphi %||% NA_real_,
    spatial_acceptance = if (dcs_model_uses(model, "uses_common_shock")) {
      fit$acc_spatial %||% fit$acckappa %||% NA_real_
    } else {
      NA_real_
    },
    mean_K = mean(fit$k_clust),
    K_ess = if (length(unique(fit$k_clust)) > 1L) {
      effective_sample_size(as.numeric(fit$k_clust))
    } else {
      NA_real_
    },
    loglik_ess = effective_sample_size(as.numeric(fit$loglik)),
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, summary_rows)
cluster_sizes <- do.call(rbind, size_rows)
parameter_table <- do.call(rbind, parameter_rows)
runtime_table <- do.call(rbind, runtime_rows)
prediction_table <- if (length(prediction_rows)) {
  do.call(rbind, prediction_rows)
} else {
  data.frame()
}
chain_table <- do.call(rbind, chain_rows)
diagnostic_table <- do.call(rbind, diagnostic_rows)
utils::write.csv(
  summary_table, file.path(root, "structural_comparison.csv"), row.names = FALSE
)
utils::write.csv(
  cluster_sizes, file.path(root, "vi_cluster_sizes.csv"), row.names = FALSE
)
utils::write.csv(
  parameter_table, file.path(root, "regime_parameter_comparison.csv"),
  row.names = FALSE
)
utils::write.csv(
  runtime_table, file.path(root, "runtime_comparison.csv"), row.names = FALSE
)
utils::write.csv(
  prediction_table, file.path(root, "prediction_surface_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  chain_table, file.path(root, "posterior_parameter_comparison.csv"),
  row.names = FALSE
)
utils::write.csv(
  diagnostic_table, file.path(root, "mcmc_diagnostics.csv"),
  row.names = FALSE
)

reference_index <- which(vapply(
  fits, function(fit) identical(fit$model %||% "dcs_sppm", "dcs_sppm"),
  logical(1L)
))[1L]
prediction_differences <- list()
if (length(reference_index) && !is.na(reference_index)) {
  reference_prediction <- fits[[reference_index]]$prediction
  for (index in seq_along(fits)) {
    model <- validate_dcs_model_type(fits[[index]]$model %||% "dcs_sppm")
    if (model == "dcs_sppm" || is.null(fits[[index]]$prediction)) next
    for (process in intersect(c("u", "w", "y"), names(reference_prediction))) {
      reference_mean <- as.numeric(reference_prediction[[process]]$mean)
      candidate_mean <- as.numeric(fits[[index]]$prediction[[process]]$mean)
      reference_sd <- as.numeric(reference_prediction[[process]]$sd)
      candidate_sd <- as.numeric(fits[[index]]$prediction[[process]]$sd)
      prediction_differences[[length(prediction_differences) + 1L]] <- data.frame(
        model = model,
        model_label = dcs_model_label(model),
        process = process,
        mean_surface_correlation = if (stats::sd(candidate_mean) > 0) {
          stats::cor(reference_mean, candidate_mean)
        } else {
          NA_real_
        },
        mean_surface_rmse = sqrt(mean((candidate_mean - reference_mean)^2)),
        mean_surface_mae = mean(abs(candidate_mean - reference_mean)),
        posterior_sd_rmse = sqrt(mean((candidate_sd - reference_sd)^2)),
        mean_sd_ratio = mean(candidate_sd) / mean(reference_sd),
        stringsAsFactors = FALSE
      )
    }
  }
}
utils::write.csv(
  if (length(prediction_differences)) {
    do.call(rbind, prediction_differences)
  } else {
    data.frame()
  },
  file.path(root, "prediction_surface_difference_vs_dcs.csv"),
  row.names = FALSE
)
print(summary_table)
