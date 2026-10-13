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

if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Install the Rcpp package.")
cache_dir <- file.path(dirname(validation_root), "tmp", "sourcecpp_predictive_regime")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(
  file.path(validation_root, "src", "predictive_regime_postprocess.cpp"),
  cacheDir = cache_dir,
  rebuild = TRUE,
  showOutput = FALSE
)

allocation <- rbind(c(1L, 1L, 2L, 2L), c(1L, 1L, 2L, 2L))
phi_star <- list(c(0.2, 0.8), c(0.25, 0.75))
tau2_star <- list(c(0.4, 1.2), c(0.5, 1.1))
mass <- c(1, 1)
observed_similarity <- rbind(
  c(-1.2, -1.0), c(-0.8, -1.1), c(0.9, 1.0), c(1.1, 0.8)
)
prediction_similarity <- rbind(c(-1, -1), c(1, 1))
vi_partition <- c(1L, 1L, 2L, 2L)
arguments <- list(
  allocation = allocation,
  phi_star = phi_star,
  tau2_star = tau2_star,
  mass = mass,
  observed_similarity = observed_similarity,
  prediction_similarity = prediction_similarity,
  vi_partition = vi_partition,
  draw_indices = 1:2,
  threads = 1L
)

full <- do.call(
  predictive_regime_summary_cpp,
  c(arguments, list(model_type = "dcs_sppm"))
)
no_shock <- do.call(
  predictive_regime_summary_cpp,
  c(arguments, list(model_type = "sppm_no_common_shock"))
)
ordinary <- do.call(
  predictive_regime_summary_cpp,
  c(arguments, list(model_type = "ppm_common_shock"))
)

stopifnot(
  isTRUE(all.equal(full$regime_probability, no_shock$regime_probability)),
  !isTRUE(all.equal(
    full$regime_probability[1L, ], full$regime_probability[2L, ]
  )),
  isTRUE(all.equal(
    ordinary$regime_probability[1L, ],
    ordinary$regime_probability[2L, ]
  ))
)

k1_arguments <- arguments
k1_arguments$allocation <- matrix(1L, nrow = 2L, ncol = 4L)
k1_arguments$phi_star <- list(0.4, 0.5)
k1_arguments$tau2_star <- list(0.8, 0.9)
k1_arguments$vi_partition <- rep.int(1L, 4L)
k1 <- do.call(
  predictive_regime_summary_cpp,
  c(k1_arguments, list(model_type = "k1_common_shock"))
)
stopifnot(
  all(k1$maximum_probability == 1),
  all(k1$normalized_entropy == 0),
  all(k1$new_regime_probability == 0),
  all(k1$dominant_regime == 1L)
)

diagnostics <- data.frame(
  assertion = c(
    "spatial partition independent of common-shock switch",
    "spatial PPM probabilities vary with prediction coordinates",
    "ordinary PPM probabilities ignore prediction coordinates",
    "fixed K=1 has certainty one and entropy zero"
  ),
  passed = TRUE
)
output_dir <- file.path(
  validation_root, "Simulations", "validation", "results"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  diagnostics,
  file.path(output_dir, "predictive_regime_switch_validation.csv"),
  row.names = FALSE
)
print(diagnostics)
