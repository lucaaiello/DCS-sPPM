if (!requireNamespace("data.table", quietly = TRUE)) stop("Install data.table first.")
library(data.table)
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_argument)) stop("Run this standalone script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]), winslash = "/")
project_dir <- dirname(script_path)
while (!file.exists(file.path(project_dir, "R", "project_paths.R"))) {
  parent <- dirname(project_dir)
  if (identical(parent, project_dir)) stop("Could not locate the manuscript project.")
  project_dir <- parent
}
repository_dir <- dirname(project_dir)
source(file.path(project_dir, "R", "unique_location_selection_helpers.R"))
preparation_dir <- file.path(project_dir, "Real_data_analysis", "preparation", "results", "unique_location_preparation")
output_dir <- file.path(preparation_dir, "detailed_audit")
audit_path <- file.path(preparation_dir, "sampling_point_selection_audit.csv")
if (!file.exists(audit_path)) stop("Run script 54 to produce the selection audit first.")
stations <- fread(audit_path, na.strings = c("", "NA"))
normalize_id <- function(values) {
  toupper(gsub(":", "_", sub("^[A-Z]{2}/", "", trimws(as.character(values))), fixed = TRUE))
}
collapse_values <- function(values) paste(sort(unique(na.omit(values))), collapse = "; ")
start_date <- as.Date("2024-01-01")
end_date <- as.Date("2024-12-31")
calendar <- seq.Date(start_date, end_date, by = "day")
message("Checking measurement procedures and their activity dates.")
metadata <- pm10_activity_metadata(fread(file.path(repository_dir, "Data", "eea_pm10", "metadata", "DataExtract.csv"),
  encoding = "UTF-8", check.names = FALSE), start_date, end_date)
if (!"metadata_lookup_id" %in% names(stations)) stations[, metadata_lookup_id := normalized_id]
metadata <- metadata[normalized_id %in% stations$metadata_lookup_id]
activity_summary <- metadata[, .(
  metadata_records = .N,
  compatible_records = sum(activity_compatible_with_2024),
  documented_compatible_records = sum(activity_compatible_with_2024 & process_start_documented),
  compatible_method_ranks = collapse_values(as.character(activity_method_rank[activity_compatible_with_2024])),
  compatible_processes = collapse_values(`Process Id`[activity_compatible_with_2024]),
  metadata_years_available = collapse_values(as.character(Year)),
  compatible_station_names = collapse_values(`Air Quality Station Name`[activity_compatible_with_2024])
), by = normalized_id]
setnames(activity_summary, "normalized_id", "metadata_lookup_id")
review <- merge(stations, activity_summary, by = "metadata_lookup_id", all.x = TRUE)
review[, competing_series := station_series > 1L | original_coordinate_series > 1L]

message("Comparing the existing daily series at shared coordinates and within station codes.")
observations <- as.data.table(readRDS(file.path(project_dir, "Data", "pm10_combined_2024.rds")))
observations[, Date := as.Date(as.numeric(unclass(Date)), origin = "1970-01-01")]
observations <- observations[Samplingpoint %in% stations$Samplingpoint & Date %in% calendar,
  .(Samplingpoint, Date, Value)]
if (anyDuplicated(observations[, .(Samplingpoint, Date)])) stop("Duplicate sampling-point/date observations.")
setkey(observations, Samplingpoint)
coordinate_pairs <- stations[, {
  if (.N > 1L) {
    pairs <- combn(sort(Samplingpoint), 2L)
    .(first_series = pairs[1L, ], second_series = pairs[2L, ])
  } else NULL
}, by = .(Longitude, Latitude)][, .(first_series, second_series)]
station_pairs <- stations[, {
  if (.N > 1L) {
    pairs <- combn(sort(Samplingpoint), 2L)
    .(first_series = pairs[1L, ], second_series = pairs[2L, ])
  } else NULL
}, by = station_key][, .(first_series, second_series)]
pairs <- unique(rbindlist(list(coordinate_pairs, station_pairs)))
pair_results <- rbindlist(lapply(seq_len(nrow(pairs)), function(pair_index) {
  first_id <- pairs$first_series[pair_index]
  second_id <- pairs$second_series[pair_index]
  first <- observations[.(first_id), .(Date, first_value = Value)]
  second <- observations[.(second_id), .(Date, second_value = Value)]
  paired <- merge(first, second, by = "Date")
  paired <- paired[is.finite(first_value) & is.finite(second_value)]
  differences <- paired$first_value - paired$second_value
  first_station <- stations[Samplingpoint == first_id]
  second_station <- stations[Samplingpoint == second_id]
  latitude_midpoint <- (first_station$Latitude + second_station$Latitude) / 2 * pi / 180
  separation_m <- 6371000 * sqrt(
    ((first_station$Latitude - second_station$Latitude) * pi / 180)^2 +
    (cos(latitude_midpoint) * (first_station$Longitude - second_station$Longitude) * pi / 180)^2)
  data.table(
    first_series = first_id, second_series = second_id,
    first_station_code = first_station$station_code, second_station_code = second_station$station_code,
    coordinate_separation_m_approx = separation_m,
    overlapping_days = nrow(paired),
    different_days = sum(abs(differences) > 1e-10),
    identical_over_observed_overlap = nrow(paired) > 0L && all(abs(differences) <= 1e-10),
    identical_observation_dates = identical(sort(first$Date[is.finite(first$first_value)]),
      sort(second$Date[is.finite(second$second_value)])),
    signed_mean_difference = if (length(differences)) mean(differences) else NA_real_,
    mean_absolute_difference = if (length(differences)) mean(abs(differences)) else NA_real_,
    maximum_absolute_difference = if (length(differences)) max(abs(differences)) else NA_real_,
    correlation = if (nrow(paired) > 2L && sd(paired$first_value) > 0 && sd(paired$second_value) > 0) {
      cor(paired$first_value, paired$second_value)
    } else NA_real_,
    first_selected = first_station$selected_location, second_selected = second_station$selected_location,
    coverage_difference_first_minus_second = first_station$observed_days - second_station$observed_days
  )
}))

message("Checking source-cache agreement and recording unmatched observations.")
source_paths <- list.files(file.path(repository_dir, "Data", "eea_pm10", "2024", "daily_chunks"),
  pattern = "[.]rds$", full.names = TRUE)
if (!length(source_paths)) stop("Country-level daily caches were not found.")
sources <- rbindlist(lapply(source_paths, function(path) {
  source_data <- as.data.table(readRDS(path))
  if (!nrow(source_data)) return(NULL)
  source_data[, .(normalized_id = normalize_id(Samplingpoint), Date = as.Date(date),
    source = as.character(source), source_value = as.numeric(value), cache_file = basename(path))]
}), fill = TRUE)
sources <- unique(sources[normalized_id %in% stations$normalized_id & Date %in% calendar])
observations[, normalized_id := normalize_id(Samplingpoint)]
source_comparison <- merge(observations, sources, by = c("normalized_id", "Date"),
  all.x = TRUE, allow.cartesian = TRUE)
source_comparison[, match_value := is.finite(Value) & is.finite(source_value) &
  abs(Value - source_value) <= 1e-7 * pmax(1, abs(Value))]
source_day_audit <- source_comparison[is.finite(Value), .(
  Value = Value[1L],
  any_cache_match = any(match_value),
  native_daily_match = any(match_value & source == "day", na.rm = TRUE),
  hourly_match = any(match_value & source == "hour", na.rm = TRUE),
  available_cache_sources = collapse_values(source),
  closest_cache_difference = if (any(is.finite(source_value))) min(abs(Value - source_value), na.rm = TRUE) else NA_real_
), by = .(Samplingpoint, Date)]
provenance_summary <- source_day_audit[, .(
  observed_days_checked = .N,
  unmatched_days = sum(!any_cache_match),
  days_without_cache_record = sum(is.na(closest_cache_difference)),
  native_daily_matching_days_checked = sum(native_daily_match),
  hourly_matching_days_checked = sum(hourly_match)
), by = Samplingpoint]
review <- merge(review, provenance_summary, by = "Samplingpoint", all.x = TRUE)
review[, audit_assessment := fifelse(
  station_coordinate_conflict, "Unresolved: check relocation, inlet coordinates or coordinate correction",
  fifelse(!metadata_matched, "Unresolved: recover the exact sampling-point metadata",
    fifelse(is.na(compatible_records) | compatible_records == 0L,
      "Unresolved: no metadata procedure compatible with 2024",
      fifelse(competing_series & (method_rank_conflict | source_rank == 3L |
        is.na(documented_compatible_records) | documented_compatible_records == 0L),
        "Unresolved: check competing methods, activity dates and provenance",
        fifelse(competing_series, "Provisionally retain the ranked choice; activity-compatible metadata exist",
          "Retain: no competing eligible series; metadata/provenance flags remain documented"))))
)]
if ("decision_applied" %in% names(review)) {
  review[decision_applied == TRUE & metadata_matched == TRUE & compatible_records > 0L,
    audit_assessment := "Resolved by documented selection decision; see resolution_evidence"]
  review[decision_applied == TRUE & metadata_matched == FALSE & station_identity_verified == TRUE,
    audit_assessment := "Retain with documented station/source evidence; exact instrument metadata remain unspecified"]
}
review[, human_decision := ""]
review[, decision_evidence := ""]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(review, file.path(output_dir, "all_series_evidence.csv"))
fwrite(review[competing_series == TRUE], file.path(output_dir, "selection_review_queue.csv"))
fwrite(review[station_coordinate_conflict == TRUE], file.path(output_dir, "coordinate_conflict_evidence.csv"))
metadata_export <- merge(metadata, stations[, .(normalized_id = metadata_lookup_id, Samplingpoint, selected_location)],
  by = "normalized_id", allow.cartesian = TRUE)
fwrite(metadata_export, file.path(output_dir, "measurement_activity_records.csv"))
fwrite(pair_results, file.path(output_dir, "paired_series_comparison.csv"))
fwrite(pair_results[identical_over_observed_overlap == TRUE], file.path(output_dir, "identical_series_pairs.csv"))
fwrite(source_day_audit[any_cache_match == FALSE], file.path(output_dir, "unmatched_daily_values.csv"))
fwrite(provenance_summary, file.path(output_dir, "source_provenance_summary.csv"))
summary_lines <- c(
  "Unique-location selection: detailed evidence audit",
  paste("Series reviewed:", nrow(review)),
  paste("Series involved in station/coordinate competition:", sum(review$competing_series)),
  paste("Series with unresolved audit assessment:", sum(startsWith(review$audit_assessment, "Unresolved"))),
  paste("Compared pairs:", nrow(pair_results)),
  paste("Pairs identical over observed overlap:", sum(pair_results$identical_over_observed_overlap)),
  paste("Finite daily observations unmatched to existing caches:", sum(!source_day_audit$any_cache_match)),
  "",
  "Start with selection_review_queue.csv and coordinate_conflict_evidence.csv.",
  "Fill human_decision with retain/change/unresolved and decision_evidence with the justification.",
  "Activity-compatible means documented dates do not exclude operation in 2024; it does not prove deployment.",
  "Identical overlap alone does not establish duplicate reporting: compare date coverage and instrument metadata.",
  "Source-cache mismatches are evidence to inspect, not proof of incorrect concentrations.",
  "The audit does not change series selection, data, sampler settings or fitted results.",
  paste("Audit outputs:", output_dir)
)
writeLines(summary_lines, file.path(output_dir, "audit_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
