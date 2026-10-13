script_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_argument)) stop("Run this script with Rscript.")
project <- dirname(normalizePath(sub("^--file=", "", script_argument[[1L]]), winslash = "/"))
while (!file.exists(file.path(project, "R", "project_paths.R"))) {
  parent <- dirname(project)
  if (identical(parent, project)) stop("Could not locate the project.")
  project <- parent
}
if (!requireNamespace("data.table", quietly = TRUE)) stop("Install data.table.")
source(file.path(project, "R", "europe_model_comparison_utils.R"))
for (reference in c("production_30000_online_grid", "production_10000_single_chain")) {
  output <- file.path(project, "Data", paste0(reference, "_preparation.rds"))
  if (file.exists(output)) stop("Refusing to overwrite frozen preparation input: ", output)
  input <- file.path(project, "Real_data_analysis", "reference_inputs", "results", reference, "pilot_fit.rds")
  message("Extracting preparation metadata from ", reference)
  fit <- readRDS(input)
  initial <- NULL
  if (reference == "production_10000_single_chain") {
    initial <- list(
      beta = as.numeric(tail(fit$beta, 1L)), xi2 = as.numeric(tail(fit$sigma2beta, 1L)),
      sigma2 = as.numeric(tail(fit$sigma2, 1L)), phi = as.numeric(tail(fit$phi, 1L)),
      tau2 = as.numeric(tail(fit$tau2, 1L)), M = as.numeric(tail(fit$Mdp, 1L)),
      kappa = as.numeric(tail(fit$kappa, 1L)), sigma2_u = as.numeric(tail(fit$sigma2_u, 1L))
    )
  }
  inputs <- list(
    schema_version = 1L,
    station_information = as.data.frame(fit$station_information),
    analysis_dates = as.Date(fit$analysis_dates), data_scaling = fit$data_scaling,
    study_bbox = fit$pilot_settings$study_bbox,
    prediction = if (is.null(fit$prediction)) NULL else list(
      grid = fit$prediction$grid,
      nearest_observed_station = fit$prediction$nearest_observed_station
    ),
    initial_values = initial,
    provenance = list(reference = reference, source_md5 = unname(tools::md5sum(input)))
  )
  stopifnot(identical(inputs$station_information, as.data.frame(fit$station_information)),
    identical(inputs$analysis_dates, as.Date(fit$analysis_dates)),
    identical(inputs$data_scaling, fit$data_scaling),
    identical(inputs$prediction$grid, fit$prediction$grid))
  rm(fit)
  invisible(gc())
  saveRDS(inputs, output, compress = "xz")
  message("Created ", output)
}
