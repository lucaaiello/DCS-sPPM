argument <- grep("^--file=", commandArgs(), value = TRUE)
project <- dirname(normalizePath(sub("^--file=", "", argument[[1L]]), winslash = "/"))
while (!file.exists(file.path(project, "R", "project_paths.R"))) project <- dirname(project)
source(file.path(project, "R", "manuscript_preparation_helpers.R"))
settings <- preparation_arguments(list(pm10 = "", weather = "", "check-only" = "false", output = file.path(project, "Data", "preparation", "results", "reconstructed")))
if (!nzchar(settings$pm10) || !nzchar(settings$weather)) stop("Supply --pm10 and --weather prepared RDS files.")
manifest_dir <- file.path(project, "Data", "preparation", "manifests")
frozen <- dget(file.path(manifest_dir, "analysis_settings.R"))
stations <- read.csv(file.path(manifest_dir, "selected_stations.csv"), stringsAsFactors = FALSE)
grid <- read.csv(file.path(manifest_dir, "prediction_grid.csv"))
if (nrow(stations) != 2989L || nrow(grid) != 5297L || anyDuplicated(stations$Samplingpoint) || anyDuplicated(stations[, c("Longitude", "Latitude")])) stop("Invalid final cohort manifest.")
dates <- as.Date(frozen$analysis_dates)
pm10 <- as.data.frame(readRDS(settings$pm10))
if (!all(c("Samplingpoint", "Date", "Value") %in% names(pm10))) stop("PM10 input must have Samplingpoint, Date and Value.")
pm10$Date <- preparation_day(pm10$Date)
pm10 <- pm10[pm10$Samplingpoint %in% stations$Samplingpoint & pm10$Date %in% dates, ]
if (anyDuplicated(pm10[, c("Samplingpoint", "Date")])) stop("Duplicate station/date observations.")
observed <- matrix(NA_real_, length(dates), nrow(stations), dimnames = list(as.character(dates), stations$Samplingpoint))
observed[cbind(match(pm10$Date, dates), match(pm10$Samplingpoint, stations$Samplingpoint))] <- pm10$Value
imputed <- observed
report <- data.frame(Samplingpoint = stations$Samplingpoint, observed_match = FALSE, imputed_match = FALSE)
for (index in seq_len(ncol(observed))) {
  valid <- which(is.finite(observed[, index]))
  if (length(valid) >= 2L) imputed[, index] <- stats::approx(valid, observed[valid, index], seq_along(dates), rule = 2)$y
  report$observed_match[index] <- identical(preparation_fingerprint(observed[, index]), stations$observed_fingerprint[index])
  report$imputed_match[index] <- identical(preparation_fingerprint(imputed[, index]), stations$imputed_fingerprint[index])
}
weather <- as.data.frame(readRDS(settings$weather))
if (!all(c("Longitude", "Latitude", "Date", "Temperature", "Precipitation", "Wind") %in% names(weather))) stop("Invalid weather schema.")
calendar_day <- preparation_day(weather$Date)
weather$Date <- as.integer(calendar_day - as.Date("2024-01-01"))
weather_checks <- read.csv(file.path(manifest_dir, "weather_fingerprints.csv"), stringsAsFactors = FALSE)
weather_key <- paste(weather$Longitude, weather$Latitude)
cell_match <- match(weather_key, paste(weather_checks$Longitude, weather_checks$Latitude))
cell_rows <- split(which(!is.na(cell_match)), cell_match[!is.na(cell_match)])
weather_report <- weather_checks[, c("Longitude", "Latitude")]
for (variable in c("Temperature", "Precipitation", "Wind")) {
  weather_report[[variable]] <- vapply(seq_len(nrow(weather_checks)), function(index) {
    rows <- cell_rows[[as.character(index)]]
    rows <- rows[order(weather$Date[rows])]
    identical(as.character(calendar_day[rows]), frozen$analysis_dates) &&
      identical(preparation_fingerprint(weather[[variable]][rows]), weather_checks[[variable]][index])
  }, logical(1L))
}
dir.create(settings$output, recursive = TRUE, showWarnings = FALSE)
write.csv(report, file.path(settings$output, "pm10_reconstruction_check.csv"), row.names = FALSE)
write.csv(weather_report, file.path(settings$output, "weather_reconstruction_check.csv"), row.names = FALSE)
if (!all(report$observed_match & report$imputed_match) || !all(as.matrix(weather_report[, c("Temperature", "Precipitation", "Wind")]))) {
  stop("Reconstruction differs from manuscript inputs. Inspect the reports; no production inputs were written.")
}
if (tolower(settings[["check-only"]]) == "true") {
  message("Validation passed; candidate observations, imputation and matched weather reproduce the manuscript fingerprints. No RDS files written.")
  quit(save = "no", status = 0L)
}
unique_input <- list(observations = pm10, station_information = stations,
  analysis_dates = dates, Y_observed = observed, Y_imputed = imputed,
  imputed_mask = is.na(observed), settings = list(minimum_coverage = 0.90, source = "frozen manuscript manifest"))
preparation <- list(schema_version = 1L, station_information = stations,
  analysis_dates = dates, data_scaling = frozen$data_scaling, study_bbox = frozen$study_bbox,
  prediction = list(grid = grid), initial_values = NULL,
  provenance = list(source = "manifest-validated reconstruction; no historical fit"))
preparation_write_rds(unique_input, file.path(settings$output, "pm10_unique_locations_2024.rds"))
preparation_write_rds(preparation, file.path(settings$output, "production_30000_online_grid_preparation.rds"))
preparation_write_rds(weather, file.path(settings$output, "era5_combined_2024.rds"))
message("Validated reconstruction of 2,989 stations, 366 days, 5,297 prediction sites and matched meteorology. Production inputs remain untouched.")
