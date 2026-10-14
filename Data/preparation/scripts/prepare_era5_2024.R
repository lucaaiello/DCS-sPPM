argument <- grep("^--file=", commandArgs(), value = TRUE)
project <- dirname(normalizePath(sub("^--file=", "", argument[[1L]]), winslash = "/"))
while (!file.exists(file.path(project, "R", "project_paths.R"))) project <- dirname(project)
source(file.path(project, "R", "manuscript_preparation_helpers.R"))
settings <- preparation_arguments(list(input = "", output = file.path(project, "Data", "preparation", "results", "era5_combined_2024.rds")))
if (!requireNamespace("ncdf4", quietly = TRUE)) stop("Install ncdf4.")
if (!nzchar(settings$input)) stop("Supply --input containing the four quarterly NetCDF directories.")
if (file.exists(settings$output)) stop("Refusing to overwrite: ", settings$output)
quarters <- c("01_02_03_2024_euro", "04_05_06_2024_euro", "07_08_09_2024_euro", "10_11_12_2024_euro")
read_variable <- function(directory, variable) {
  paths <- list.files(directory, pattern = "[.]nc$", full.names = TRUE)
  candidates <- list()
  for (path in paths) {
    handle <- ncdf4::nc_open(path)
    tryCatch({
      if (variable %in% names(handle$var)) {
        dimensions <- vapply(handle$var[[variable]]$dim, function(dimension) dimension$name, character(1L))
        if (length(dimensions) != 3L || !all(c("longitude", "latitude", "valid_time") %in% dimensions)) stop("Unsupported NetCDF dimensions in ", path)
        values <- ncdf4::ncvar_get(handle, variable, collapse_degen = FALSE)
        values <- aperm(values, match(c("longitude", "latitude", "valid_time"), dimensions))
        units <- ncdf4::ncatt_get(handle, "valid_time", "units")$value
        pieces <- strsplit(units, " since ", fixed = TRUE)[[1L]]
        multiplier <- if (length(pieces) == 2L) switch(pieces[1L], seconds = 1, hours = 3600, days = 86400, NA_real_) else NA_real_
        if (!is.finite(multiplier)) stop("Unsupported time units: ", units)
        origin <- as.POSIXct(pieces[2L], tz = "UTC")
        if (is.na(origin)) stop("Invalid NetCDF time origin.")
        candidates[[length(candidates) + 1L]] <- list(values = values,
          longitude = ncdf4::ncvar_get(handle, "longitude"), latitude = ncdf4::ncvar_get(handle, "latitude"),
          dates = as.Date(origin + multiplier * ncdf4::ncvar_get(handle, "valid_time"), tz = "UTC"))
      }
    }, finally = ncdf4::nc_close(handle))
  }
  if (length(candidates) != 1L) stop("Expected exactly one ", variable, " array in ", directory)
  candidates[[1L]]
}
blocks <- lapply(quarters, function(quarter) {
  directory <- file.path(settings$input, quarter)
  variables <- lapply(c("t2m", "tp", "u10", "v10"), function(variable) read_variable(directory, variable))
  reference <- variables[[1L]]
  for (variable in variables[-1L]) {
    if (!identical(variable$longitude, reference$longitude) || !identical(variable$latitude, reference$latitude) || !identical(variable$dates, reference$dates)) stop("NetCDF variables have different coordinates or dates.")
  }
  result <- expand.grid(Longitude = reference$longitude, Latitude = reference$latitude,
    Date = as.integer(reference$dates - as.Date("2024-01-01")))
  result$Temperature <- as.vector(reference$values) - 273.15
  result$Precipitation <- as.vector(variables[[2L]]$values)
  result$Wind <- as.vector(sqrt(variables[[3L]]$values^2 + variables[[4L]]$values^2))
  result
})
weather <- do.call(rbind, blocks)
weather <- weather[order(weather$Longitude, weather$Latitude, weather$Date), ]
rownames(weather) <- NULL
if (!identical(sort(unique(weather$Date)), 0:365) || anyDuplicated(weather[, c("Longitude", "Latitude", "Date")])) stop("Weather must contain every 2024 UTC day exactly once per grid cell.")
if (any(!is.finite(as.matrix(weather)))) stop("Non-finite meteorological values.")
preparation_write_rds(weather, settings$output)
message("Prepared all four quarters. Wind is the magnitude of daily-mean u10/v10 components; precipitation retains the supplied tp daily statistic, with no daily-total conversion.")
