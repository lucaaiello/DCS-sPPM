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

csv_setting <- function(name, default) {
  values <- trimws(strsplit(Sys.getenv(name, default), ",", fixed = TRUE)[[1L]])
  values[nzchar(values)]
}
flag <- function(name, default) {
  tolower(Sys.getenv(name, if (default) "true" else "false")) %in%
    c("true", "1")
}
models <- csv_setting(
  "DCS_REAL_FULL_MODELS",
  "dcs_sppm,ppm_common_shock,k1_common_shock,sppm_no_common_shock"
)
models <- vapply(models, validate_dcs_model_type, character(1L))
results_name <- Sys.getenv(
  "DCS_REAL_FULL_RESULTS_NAME", "production_30000_unique_locations"
)
render_outputs <- flag("DCS_REAL_FULL_RENDER", TRUE)
render_daily_maps <- flag("DCS_REAL_FULL_RENDER_DAILY_MAPS", FALSE)
summarize_outputs <- flag("DCS_REAL_FULL_SUMMARIZE", TRUE)
dry_run <- flag("DCS_REAL_FULL_DRY_RUN", FALSE)

rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
run_child <- function(script, environment = character()) {
  script_path <- file.path(validation_root, script)
  if (dry_run) {
    cat(
      "DRY RUN: ", paste(environment, collapse = " "),
      if (length(environment)) " " else "", rscript, " ", script_path, "\n",
      sep = ""
    )
    return(invisible(0L))
  }
  environment_names <- sub("=.*$", "", environment)
  environment_values <- sub("^[^=]*=", "", environment)
  previous_values <- Sys.getenv(environment_names, unset = NA_character_)
  if (length(environment_names)) {
    do.call(Sys.setenv, as.list(stats::setNames(
      environment_values, environment_names
    )))
  }
  on.exit({
    present <- !is.na(previous_values)
    if (any(present)) {
      do.call(Sys.setenv, as.list(stats::setNames(
        previous_values[present], environment_names[present]
      )))
    }
    if (any(!present)) Sys.unsetenv(environment_names[!present])
  }, add = TRUE)
  status <- system2(
    rscript,
    args = shQuote(script_path)
  )
  if (!identical(status, 0L)) {
    stop(script, " failed with exit status ", status, ".")
  }
}

for (model in models) {
  cat("\n===== Full European analysis: ", model, " =====\n", sep = "")
  run_child(
    "Real_data_analysis/model_comparison/scripts/run_full_analysis.R",
    paste0("DCS_REAL_FULL_MODELS=", model)
  )

  if (render_outputs) {
    result_dir <- file.path(
      validation_results_dir, "Real_data_analysis", "model_comparison", "results",
      results_name, model
    )
    resolved_result_dir <- if (dry_run) {
      normalizePath(dirname(result_dir), winslash = "/", mustWork = FALSE) |>
        file.path(basename(result_dir))
    } else {
      normalizePath(result_dir, winslash = "/", mustWork = TRUE)
    }
    render_environment <- c(
      paste0("DCS_GRID_RESULT_DIR=", resolved_result_dir),
      paste0(
        "DCS_RENDER_PREDICTIVE_MAPS=",
        if (render_daily_maps) "true" else "false"
      ),
      "DCS_REUSE_VI=true"
    )
    run_child(
      "Real_data_analysis/maps/scripts/render_prediction_and_clustering.R",
      render_environment
    )
    run_child(
      "Real_data_analysis/maps/scripts/generate_predictive_regime_maps.R",
      c(render_environment, "DCS_REUSE_REGIME_SUMMARY=true")
    )
  }
}

if (summarize_outputs) {
  cat("\n===== Summarizing full European model comparison =====\n")
  run_child("Real_data_analysis/model_comparison/scripts/summarize_full_analysis.R")
}

cat(
  "\nFull reproducibility workflow completed under:\n",
  file.path(
    validation_results_dir, "Real_data_analysis", "model_comparison", "results", results_name
  ),
  "\n",
  sep = ""
)
