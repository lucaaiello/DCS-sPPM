script_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_argument)) stop("Run this script with Rscript.")
project <- dirname(normalizePath(sub("^--file=", "", script_argument[[1L]]), winslash = "/"))
while (!file.exists(file.path(project, "R", "project_paths.R"))) {
  parent <- dirname(project)
  if (identical(parent, project)) stop("Could not locate the project.")
  project <- parent
}
source(file.path(project, "R", "simulate_dcs_sppm.R"))
source(file.path(project, "R", "prediction_experiment_utils.R"))
source(file.path(project, "R", "model_comparison_utils.R"))
source(file.path(project, "R", "europe_model_comparison_utils.R"))
Sys.setenv(DCS_REAL_FULL_DATA_FILE = file.path(project, "Data", "pm10_unique_locations_2024.rds"))
prepared <- prepare_europe_model_comparison_data(project, project,
  holdout_result_name = NULL, preparation_name = "production_30000_online_grid")
prediction <- prepare_europe_full_prediction_inputs(prepared)
stopifnot(ncol(prepared$Y) == 2989L, nrow(prepared$Y) == 366L,
  nrow(prediction$grid) == 5297L, all(is.finite(prepared$Y)),
  all(is.finite(prepared$X)), all(is.finite(prediction$X)),
  identical(prediction$grid, as.data.frame(prepared$preparation$prediction$grid)),
  is.null(prepared$production_fit), is.null(prepared$preparation$beta),
  is.null(prepared$preparation$cluster))
for (model in c("dcs_sppm", "ppm_common_shock", "k1_common_shock", "sppm_no_common_shock")) {
  initial <- europe_initial_values(prepared, model)
  stopifnot(identical(initial$beta, rep(0, ncol(prepared$X))),
    length(initial$phi) == ncol(prepared$Y), all(initial$phi == 0.5))
}
message("Preparation validated: 2,989 stations, 366 days, 5,297 prediction sites; no historical fit was read.")
