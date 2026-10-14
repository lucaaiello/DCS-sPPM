argument <- grep("^--file=", commandArgs(), value = TRUE)
project <- dirname(normalizePath(sub("^--file=", "", argument[[1L]]), winslash = "/"))
while (!file.exists(file.path(project, "R", "project_paths.R"))) project <- dirname(project)
source(file.path(project, "R", "manuscript_preparation_helpers.R"))
settings <- preparation_arguments(list(output = file.path(project, "Data", "preparation", "manifests")))
dir.create(settings$output, recursive = TRUE, showWarnings = FALSE)
paths <- file.path(settings$output, c("analysis_settings.R", "selected_stations.csv", "prediction_grid.csv", "weather_fingerprints.csv"))
if (any(file.exists(paths))) stop("Refusing to overwrite an existing manuscript manifest.")
reference <- readRDS(file.path(project, "Data", "production_30000_online_grid_preparation.rds"))
selected <- readRDS(file.path(project, "Data", "pm10_unique_locations_2024.rds"))
stations <- as.data.frame(selected$station_information)
inventory <- as.data.frame(reference$station_information)
matching <- vapply(seq_len(nrow(stations)), function(index) {
  same <- which(inventory$Samplingpoint == stations$Samplingpoint[index])
  if (length(same) == 1L) return(same)
  which.min((inventory$Longitude - stations$Longitude[index])^2 + (inventory$Latitude - stations$Latitude[index])^2)
}, integer(1L))
weather_columns <- c("weather_longitude", "weather_latitude", "weather_cell")
manifest <- stations[, c("Samplingpoint", "Longitude", "Latitude", "station_key", "selection_reason")]
manifest[weather_columns] <- inventory[matching, weather_columns]
manifest$observed_fingerprint <- vapply(seq_len(nrow(stations)), function(index) preparation_fingerprint(selected$Y_observed[, index]), character(1L))
manifest$imputed_fingerprint <- vapply(seq_len(nrow(stations)), function(index) preparation_fingerprint(selected$Y_imputed[, index]), character(1L))
frozen <- list(schema_version = 1L, analysis_dates = as.character(selected$analysis_dates),
  data_scaling = reference$data_scaling, study_bbox = reference$study_bbox,
  fingerprint_format = "MD5 of one %.12g value per line, NA for missing values; UTC daily order")
dput(frozen, paths[1L], control = "all")
write.csv(manifest, paths[2L], row.names = FALSE, na = "")
write.csv(reference$prediction$grid, paths[3L], row.names = FALSE)
weather <- readRDS(file.path(project, "Data", "era5_combined_2024.rds"))
weather$calendar_day <- preparation_day(weather$Date)
cells <- unique(manifest[weather_columns])
weather_key <- paste(weather$Longitude, weather$Latitude)
cell_match <- match(weather_key, paste(cells$weather_longitude, cells$weather_latitude))
cell_rows <- split(which(!is.na(cell_match)), cell_match[!is.na(cell_match)])
checks <- lapply(seq_len(nrow(cells)), function(index) {
  cell <- cells[index, ]
  rows <- cell_rows[[as.character(index)]]
  rows <- rows[order(weather$calendar_day[rows])]
  if (!identical(as.character(weather$calendar_day[rows]), frozen$analysis_dates)) stop("Weather dates do not match the analysed sample.")
  data.frame(Longitude = cell$weather_longitude, Latitude = cell$weather_latitude,
    Temperature = preparation_fingerprint(weather$Temperature[rows]),
    Precipitation = preparation_fingerprint(weather$Precipitation[rows]),
    Wind = preparation_fingerprint(weather$Wind[rows]))
})
write.csv(do.call(rbind, checks), paths[4L], row.names = FALSE)
message("Exported station decisions, grid, scaling constants and numeric fingerprints; no concentrations or posterior arrays.")
