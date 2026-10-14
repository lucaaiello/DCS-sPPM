argument <- grep("^--file=", commandArgs(), value = TRUE)
project <- dirname(normalizePath(sub("^--file=", "", argument[[1L]]), winslash = "/"))
while (!file.exists(file.path(project, "R", "project_paths.R"))) project <- dirname(project)
source(file.path(project, "R", "manuscript_preparation_helpers.R"))
settings <- preparation_arguments(list(hourly = "", daily = "", metadata = "",
  output = file.path(project, "Data", "preparation", "results", "pm10_combined_2024.rds")))
packages <- c("arrow", "dplyr", "data.table")
missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Install: ", paste(missing, collapse = ", "))
if (!nzchar(settings$hourly) || !nzchar(settings$daily) || !nzchar(settings$metadata)) stop("Supply --hourly parquet directory, --daily CSV and --metadata CSV.")
if (file.exists(settings$output)) stop("Refusing to overwrite: ", settings$output)
start <- as.POSIXct("2024-01-01", tz = "UTC")
end <- as.POSIXct("2025-01-01", tz = "UTC")
hourly <- arrow::open_dataset(settings$hourly, format = "parquet")
required <- c("Samplingpoint", "Start", "Value", "Validity")
if (!all(required %in% names(hourly))) stop("Hourly input is missing required columns.")
hourly <- dplyr::filter(hourly, Start >= start, Start < end, Value > 0, Validity == 1)
hourly <- dplyr::mutate(hourly, Date = as.Date(Start))
hourly <- dplyr::group_by(hourly, Samplingpoint, Date)
hourly <- dplyr::collect(dplyr::summarise(hourly, Value = mean(Value, na.rm = TRUE), valid_hours = dplyr::n(), .groups = "drop"))
hourly <- data.table::as.data.table(hourly)
hourly[, Date := preparation_day(Date)]
daily <- data.table::fread(settings$daily)
if (!all(required %in% names(daily))) stop("Daily input is missing required columns.")
daily[, Start := as.POSIXct(Start, tz = "UTC")]
daily <- daily[Start >= start & Start < end & Value > 0 & Value < 500 & Validity == 1,
  .(Samplingpoint, Date = as.Date(Start, tz = "UTC"), Value)]
metadata <- data.table::fread(settings$metadata)
if (!all(c("Sampling.Point.Id", "Longitude", "Latitude") %in% names(metadata))) stop("Metadata must use the supplied Sampling.Point.Id schema.")
metadata <- unique(metadata[, .(Sampling.Point.Id, Longitude, Latitude)], by = "Sampling.Point.Id")
combined <- data.table::rbindlist(list(hourly, daily), fill = TRUE)
combined[, SamplingPoint_clean := sub("^[A-Z]{2}/", "", Samplingpoint)]
combined[, source_order := seq_len(.N)]
data.table::setorderv(combined, c("SamplingPoint_clean", "Date", "Value", "source_order"))
combined <- combined[, .SD[1L], by = .(SamplingPoint_clean, Date)]
combined[, source_order := NULL]
combined <- merge(combined, metadata, by.x = "SamplingPoint_clean", by.y = "Sampling.Point.Id", all.x = TRUE, sort = FALSE)
preparation_write_rds(combined, settings$output)
message("Prepared using the supplied source-script rules: positive values, Validity=1, daily values below 500, hourly means and minimum across products. Verification and an 18-hour threshold are not added.")
message("Candidate output: ", settings$output, ". Validate against the manuscript manifest before use.")
