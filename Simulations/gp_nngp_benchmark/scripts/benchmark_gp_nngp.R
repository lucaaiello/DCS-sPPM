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

if (!requireNamespace("Rcpp", quietly = TRUE)) stop("The `Rcpp` package is required.")
cache_dir <- file.path(dirname(validation_root), "tmp", "sourcecpp_gp_nngp_benchmark")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(
  file.path(validation_root, "src", "benchmark_gp_nngp.cpp"),
  cacheDir = cache_dir,
  rebuild = FALSE
)

sample_sizes <- c(100L, 500L, 1000L, 2000L, 3000L, 4000L)
replications <- 5L
n_neighbors <- 10L
results_dir <- simulation_benchmark_results_dir
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
results_file <- file.path(results_dir, "timings.csv")
results <- if (file.exists(results_file)) {
  utils::read.csv(results_file, stringsAsFactors = FALSE)
} else {
  data.frame()
}

for (sample_size in sample_sizes) {
  for (replicate in seq_len(replications)) {
    set.seed(20300000L + sample_size + replicate)
    locations <- matrix(stats::runif(2L * sample_size), ncol = 2L)
    field <- stats::rnorm(sample_size)
    neighbors <- make_previous_neighbor_matrix(locations, n_neighbors)
    for (method in c("exact_gp", "nngp")) {
      already_done <- nrow(results) && any(
        results$n == sample_size & results$replicate == replicate &
          results$method == method
      )
      if (already_done) next
      message(
        "Benchmarking ", method, ", n=", sample_size,
        ", replicate=", replicate, "."
      )
      timing <- benchmark_gp_nngp_cpp(
        locations,
        neighbors,
        field,
        method,
        variance = 1,
        range = 0.2,
        nu = 0.5
      )
      row <- data.frame(
        method = method,
        n = sample_size,
        replicate = replicate,
        n_neighbors = if (method == "nngp") n_neighbors else NA_integer_,
        setup_seconds = timing$setup_seconds,
        evaluation_seconds = timing$evaluation_seconds,
        total_seconds = timing$total_seconds,
        storage_mb = timing$storage_mb,
        log_density = timing$log_density,
        stringsAsFactors = FALSE
      )
      results <- rbind(results, row)
      utils::write.csv(results, results_file, row.names = FALSE)
    }
  }
}

groups <- split(results, list(results$method, results$n), drop = TRUE)
summary <- do.call(rbind, lapply(groups, function(group) {
  data.frame(
    method = group$method[1],
    n = group$n[1],
    replications = nrow(group),
    median_setup_seconds = stats::median(group$setup_seconds),
    median_evaluation_seconds = stats::median(group$evaluation_seconds),
    median_total_seconds = stats::median(group$total_seconds),
    median_storage_mb = stats::median(group$storage_mb),
    stringsAsFactors = FALSE
  )
}))
exact <- summary[summary$method == "exact_gp", c("n", "median_total_seconds")]
names(exact)[2] <- "exact_total_seconds"
nngp <- summary[summary$method == "nngp", , drop = FALSE]
comparison <- merge(nngp, exact, by = "n", sort = TRUE)
comparison$speedup_over_exact <-
  comparison$exact_total_seconds / comparison$median_total_seconds
utils::write.csv(summary, file.path(results_dir, "summary.csv"), row.names = FALSE)
utils::write.csv(comparison, file.path(results_dir, "comparison.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
print(comparison, row.names = FALSE)
